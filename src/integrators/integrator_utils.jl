@inline function loopheader!(integrator::SDEIntegrator)
  # Apply right after iterators / callbacks

  # Accept or reject the step
  if integrator.iter > 0
    if ((integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive) && !integrator.force_stepfail
      integrator.success_iter += 1
      apply_step!(integrator)
    elseif integrator.opts.adaptive && !integrator.accept_step
      if integrator.isout
        integrator.dtnew = integrator.dt*integrator.opts.qmin
      elseif !integrator.force_stepfail
        integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
      end
      choose_algorithm!(integrator,integrator.cache)
      fix_dtnew_at_bounds!(integrator)
      modify_dtnew_for_tstops!(integrator)
      reject_step!(integrator.W,integrator.dtnew)
      integrator.dt = integrator.dtnew
      integrator.sqdt = sqrt(abs(integrator.dt))
    end
  end

  integrator.iter += 1
  integrator.force_stepfail = false
end

@inline function fix_dtnew_at_bounds!(integrator)
  integrator.dtnew = integrator.tdir*min(abs(integrator.opts.dtmax),abs(integrator.dtnew))
  integrator.dtnew = integrator.tdir*max(abs(integrator.dtnew),abs(integrator.opts.dtmin))
end

@inline function modify_dt_for_tstops!(integrator)
  tstops = integrator.opts.tstops
  @fastmath if !isempty(tstops)
    if integrator.opts.adaptive
      if integrator.tdir > 0
        integrator.dt = min(abs(integrator.dt),abs(top(tstops)-integrator.t)) # step! to the end
      else
        integrator.dt = -min(abs(integrator.dt),abs(top(tstops)-integrator.t))
      end
    elseif integrator.dtcache == zero(integrator.t) && integrator.dtchangeable # Use integrator.opts.tstops
      integrator.dt = integrator.tdir*abs(top(tstops)-integrator.t)
    elseif integrator.dtchangeable && !integrator.force_stepfail
      # always try to step! with dtcache, but lower if a tstops
      integrator.dt = @fastmath integrator.tdir*min(abs(integrator.dtcache),abs(top(tstops)-integrator.t)) # step! to the end
    end
  end
end

@inline function modify_dtnew_for_tstops!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    if integrator.tdir > 0
      integrator.dt = min(abs(integrator.dtnew),abs(top(tstops)-integrator.t)) # step! to the end
    else
      integrator.dt = -min(abs(integrator.dtnew),abs(top(tstops)-integrator.t))
    end
  end
end

last_step_failed(integrator::SDEIntegrator) =
  integrator.last_stepfail && !integrator.opts.adaptive

@inline function savevalues!(integrator::SDEIntegrator,force_save=false)::Tuple{Bool,Bool}
  saved, savedexactly = false, false
  !integrator.opts.save_on && return saved, savedexactly
  while !isempty(integrator.opts.saveat) && integrator.tdir*top(integrator.opts.saveat) <= integrator.tdir*integrator.t # Perform saveat
    integrator.saveiter += 1; saved = true
    curt = pop!(integrator.opts.saveat)
    if curt!=integrator.t # If <t, interpolate
      Θ = (curt - integrator.tprev)/integrator.dt
      val = sde_interpolant(Θ,integrator,integrator.opts.save_idxs,Val{0}) # out of place, but force copy later
      save_val = val
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,save_val,Val{false})
      if typeof(integrator.alg) <: StochasticDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      savedexactly = true
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      if integrator.opts.save_idxs === nothing
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      else
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
      end
      if typeof(integrator.alg) <: Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    end
  end
  if force_save || integrator.opts.save_everystep
    integrator.saveiter += 1; saved, savedexactly = true, true
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if typeof(integrator.alg) <: Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
  end
  return saved, savedexactly
end

@inline function loopfooter!(integrator::SDEIntegrator)
  ttmp = integrator.t + integrator.dt
  if integrator.force_stepfail
    if integrator.opts.adaptive
      integrator.dtnew = integrator.dt/integrator.opts.failfactor
    elseif integrator.last_stepfail
      return
    end
    integrator.last_stepfail = true
    integrator.accept_step = false
  elseif integrator.opts.adaptive
    @fastmath integrator.q11 = integrator.EEst^integrator.opts.beta1
    @fastmath integrator.q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    @fastmath integrator.q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma))
    @fastmath integrator.dtnew = integrator.dt/integrator.q
    integrator.isout = integrator.opts.isoutofdomain(integrator.u,integrator.p,ttmp)
    integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0) || (integrator.opts.force_dtmin && integrator.dt <= integrator.opts.dtmin)
    if integrator.accept_step # Accepted
      integrator.last_stepfail = false
      integrator.tprev = integrator.t
      if typeof(integrator.t)<:AbstractFloat && !isempty(integrator.opts.tstops)
        tstop = top(integrator.opts.tstops)
        @fastmath abs(ttmp - tstop) < 10eps(integrator.t) ? (integrator.t = tstop) : (integrator.t = ttmp)
      else
        integrator.t = ttmp
      end
      calc_dt_propose!(integrator)
      handle_callbacks!(integrator)
    end
  else # Non adaptive
    integrator.tprev = integrator.t
    if typeof(integrator.t)<:AbstractFloat && !isempty(integrator.opts.tstops)
      tstop = top(integrator.opts.tstops)
      # For some reason 10eps(integrator.t) is slow here
      # TODO: Allow higher precision but profile
      @fastmath abs(ttmp - tstop) < 10eps(max(integrator.t,tstop)) ? (integrator.t = tstop) : (integrator.t = ttmp)
    else
      integrator.t = ttmp
    end
    integrator.last_stepfail = false
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    handle_callbacks!(integrator)
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    @logmsg(-1,
    integrator.opts.progress_name,
    _id = :StochasticDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress=integrator.t/integrator.sol.prob.tspan[2])
  end
end

@inline function calc_dt_propose!(integrator)
  integrator.qold = max(integrator.EEst,integrator.opts.qoldinit)
  if integrator.tdir > 0
    integrator.dtpropose = min(integrator.opts.dtmax,integrator.dtnew)
  else
    integrator.dtpropose = max(integrator.opts.dtmax,integrator.dtnew)
  end
  if integrator.tdir > 0
    integrator.dtpropose = max(integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
  else
    integrator.dtpropose = min(integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
  end
end

@inline function solution_endpoint_match_cur_integrator!(integrator)
  if integrator.opts.save_end && (integrator.saveiter == 0 || integrator.sol.t[integrator.saveiter] != integrator.t)
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if integrator.opts.save_idxs === nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
  end
  if integrator.W.curt != integrator.t
    accept_step!(integrator.W,integrator.dt,false)
  end
  save_noise!(integrator.W)
end

@inline function postamble!(integrator::SDEIntegrator)
  solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  if integrator.opts.progress
    @logmsg(-1,
    integrator.opts.progress_name,
    _id = :StochasticDiffEq,
    message=integrator.opts.progress_message(integrator.dt,integrator.u,integrator.p,integrator.t),
    progress="done")
  end
  return nothing
end

@inline function handle_callbacks!(integrator)
  discrete_callbacks = integrator.opts.callback.discrete_callbacks
  continuous_callbacks = integrator.opts.callback.continuous_callbacks
  atleast_one_callback = false

  continuous_modified = false
  discrete_modified = false
  saved_in_cb = false
  if !(typeof(continuous_callbacks)<:Tuple{})
    time,upcrossing,event_occurred,idx,counter =
              DiffEqBase.find_first_continuous_callback(integrator,continuous_callbacks...)
    if event_occurred
      integrator.event_last_time = idx
      continuous_modified,saved_in_cb = DiffEqBase.apply_callback!(integrator,continuous_callbacks[idx],time,upcrossing)
    else
      integrator.event_last_time = 0
    end
  end
  if !(typeof(discrete_callbacks)<:Tuple{})
    discrete_modified,saved_in_cb = DiffEqBase.apply_discrete_callback!(integrator,discrete_callbacks...)
  end
  if !saved_in_cb
    savevalues!(integrator)
  end

  integrator.u_modified = continuous_modified || discrete_modified
  if integrator.u_modified
    handle_callback_modifiers!(integrator)
  end
end

@inline function handle_callback_modifiers!(integrator::SDEIntegrator)
  #integrator.reeval_fsal = true
end

@inline function apply_step!(integrator)
  if isinplace(integrator.sol.prob)
    recursivecopy!(integrator.uprev,integrator.u)
  else
    integrator.uprev = integrator.u
  end
  integrator.dt = integrator.dtpropose
  modify_dt_for_tstops!(integrator)
  accept_step!(integrator.W,integrator.dt)
  integrator.dt = integrator.W.dt
  integrator.sqdt = @fastmath sqrt(abs(integrator.dt)) # It can change dt, like in RSwM1
end

@inline function handle_tstop!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    t = integrator.t
    ts_top = top(tstops)
    if t == ts_top
      pop!(tstops)
      integrator.just_hit_tstop = true
    elseif integrator.tdir*t > integrator.tdir*ts_top
      if !integrator.dtchangeable
        change_t_via_interpolation!(integrator, pop!(tstops), Val{true})
        integrator.just_hit_tstop = true
      else
        error("Something went wrong. Integrator stepped past tstops but the algorithm was dtchangeable. Please report this error.")
      end
    end
  end
end

@inline function update_noise!(integrator,scaling_factor=integrator.sqdt)
  if isinplace(integrator.noise)
    integrator.noise(integrator.ΔW,integrator)
    rmul!(integrator.ΔW,scaling_factor)
    if alg_needs_extra_process(integrator.alg)
      integrator.noise(integrator.ΔZ,integrator)
      rmul!(integrator.ΔZ,scaling_factor)
    end
  else
    if typeof(integrator.u) <: AbstractArray
      integrator.ΔW .= scaling_factor.*integrator.noise(size(integrator.u),integrator)
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZ .= scaling_factor.*integrator.noise(size(integrator.u),integrator)
      end
    else
      integrator.ΔW = scaling_factor*integrator.noise(integrator)
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZ = scaling_factor*integrator.noise(integrator)
      end
    end
  end
end

@inline function generate_tildes(integrator,add1,add2,scaling)
  if isinplace(integrator.noise)
    integrator.noise(integrator.ΔWtilde,integrator)
    if add1 != 0
      #@.. integrator.ΔWtilde = add1 + scaling*integrator.ΔWtilde
      @tight_loop_macros for i in eachinex(integrator.u)
        @inbounds integrator.ΔWtilde[i] = add1[i] + scaling*integrator.ΔWtilde[i]
      end
    else
      #@.. integrator.ΔWtilde = scaling*integrator.ΔWtilde
      @tight_loop_macros for i in eachinex(integrator.u)
        @inbounds integrator.ΔWtilde[i] = scaling*integrator.ΔWtilde[i]
      end
    end
    if alg_needs_extra_process(integrator.alg)
      integrator.noise(integrator.ΔZtilde,integrator)
      if add2 != 0
        #@.. integrator.ΔZtilde = add2 + scaling*integrator.ΔZtilde
        @tight_loop_macros for i in eachinex(integrator.u)
          @inbounds integrator.ΔZtilde[i] = add2[i] + scaling*integrator.ΔZtilde[i]
        end
      else
        #@.. integrator.ΔZtilde = scaling*integrator.ΔZtilde
        @tight_loop_macros for i in eachinex(integrator.u)
          @inbounds integrator.ΔZtilde[i] = scaling*integrator.ΔZtilde[i]
        end
      end
    end
  else
    if typeof(integrator.u) <: AbstractArray
      if add1 != 0
        integrator.ΔWtilde = add1 .+ scaling.*integrator.noise(size(integrator.u),integrator)
      else
        integrator.ΔWtilde = scaling.*integrator.noise(size(integrator.u),integrator)
      end
      if alg_needs_extra_process(integrator.alg)
        if add2 != 0
          integrator.ΔZtilde = add2 .+ scaling.*integrator.noise(size(integrator.u),integrator)
        else
          integrator.ΔZtilde = scaling.*integrator.noise(size(integrator.u),integrator)
        end
      end
    else
      integrator.ΔWtilde = add1 + scaling*integrator.noise(integrator)
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZtilde = add2 + scaling*integrator.noise(integrator)
      end
    end
  end
end

@inline initialize!(integrator,cache::StochasticDiffEqCache,f=integrator.f) = nothing
