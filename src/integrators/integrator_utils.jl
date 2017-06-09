@inline function loopheader!(integrator::SDEIntegrator)
  # Apply right after iterators / callbacks

  # Accept or reject the step
  if integrator.iter > 0
    if (integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive
      apply_step!(integrator)
    elseif integrator.opts.adaptive && !integrator.accept_step
      if integrator.isout
        integrator.dtnew = integrator.dt*integrator.opts.qmin
      else
        integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
      end
      fix_dtnew_at_bounds!(integrator)
      reject_step!(integrator.W,integrator.dtnew)
      integrator.dt = integrator.dtnew
      integrator.sqdt = sqrt(abs(integrator.dt))
    end
  end

  integrator.iter += 1
  choose_algorithm!(integrator,integrator.cache)
end

@inline function fix_dtnew_at_bounds!(integrator)
  integrator.dtnew = integrator.tdir*min(abs(integrator.opts.dtmax),abs(integrator.dtnew))
  integrator.dtnew = integrator.tdir*max(abs(integrator.dtnew),abs(integrator.opts.dtmin))
end

@inline function modify_dt_for_tstops!(integrator)
  tstops = integrator.opts.tstops
  if !isempty(tstops)
    if integrator.opts.adaptive
      if integrator.tdir > 0
        integrator.dt = min(abs(integrator.dt),abs(top(tstops)-integrator.t)) # step! to the end
      else
        integrator.dt = -min(abs(integrator.dt),abs(top(tstops)-integrator.t))
      end
    elseif integrator.dtcache == zero(integrator.t) && integrator.dtchangeable # Use integrator.opts.tstops
      integrator.dt = integrator.tdir*abs(top(tstops)-integrator.t)
    elseif integrator.dtchangeable # always try to step! with dtcache, but lower if a tstops
      integrator.dt = integrator.tdir*min(abs(integrator.dtcache),abs(top(tstops)-integrator.t)) # step! to the end
    end
  end
end

@def sde_exit_condtions begin
  if integrator.iter > integrator.opts.maxiters
    if integrator.opts.verbose
      warn("Max Iters Reached. Aborting")
    end
    postamble!(integrator)
    integrator.sol.retcode = :MaxIters
    return integrator.sol
  end
  if !integrator.opts.force_dtmin && integrator.opts.adaptive && abs(integrator.dt) <= abs(integrator.opts.dtmin)
    if integrator.opts.verbose
      warn("dt <= dtmin. Aborting. If you would like to force continuation with dt=dtmin, set force_dtmin=true")
    end
    postamble!(integrator)
    integrator.sol.retcode = :DtLessThanMin
    return integrator.sol
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.t,integrator.u)
    if integrator.opts.verbose
      warn("Instability detected. Aborting")
    end
    postamble!(integrator)
    integrator.sol.retcode = :Unstable
    return integrator.sol
  end
end


@inline function savevalues!(integrator::SDEIntegrator)
  while !isempty(integrator.opts.saveat) && integrator.tdir*top(integrator.opts.saveat) <= integrator.tdir*integrator.t # Perform saveat
    integrator.saveiter += 1
    curt = pop!(integrator.opts.saveat)
    if integrator.opts.saveat!=integrator.t # If <t, interpolate
      Θ = (curt - integrator.tprev)/integrator.dt
      val = sde_interpolant(Θ,integrator,integrator.opts.save_idxs,Val{0}) # out of place, but force copy later
      save_val = val
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,save_val,Val{false})
      if typeof(integrator.alg) <: StochasticDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    else # ==t, just save
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      if integrator.opts.save_idxs == nothing
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      else
        copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
      end
      if typeof(alg) <: StochasticDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
    end
  end
  if integrator.opts.save_everystep && integrator.iter%integrator.opts.timeseries_steps==0
    integrator.saveiter += 1
    if integrator.opts.save_idxs == nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    #if typeof(integrator.alg) <: StochasticDiffEqCompositeAlgorithm
    #  copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    #end
  end
end

@inline function loopfooter!(integrator::SDEIntegrator)
  if integrator.opts.adaptive
    @fastmath integrator.q11 = integrator.EEst^integrator.opts.beta1
    @fastmath integrator.q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    @fastmath integrator.q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma))
    @fastmath integrator.dtnew = integrator.dt/integrator.q
    ttmp = integrator.t + integrator.dt
    integrator.isout = integrator.opts.isoutofdomain(ttmp,integrator.u)
    integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0) || (integrator.opts.force_dtmin && integrator.dt <= integrator.opts.dtmin)
    if integrator.accept_step # Accepted
      integrator.tprev = integrator.t
      if typeof(integrator.t)<:AbstractFloat && !isempty(integrator.opts.tstops)
        tstop = top(integrator.opts.tstops)
        abs(ttmp - tstop) < 10eps(integrator.t) ? (integrator.t = tstop) : (integrator.t = ttmp)
      else
        integrator.t = ttmp
      end
      calc_dt_propose!(integrator)
      handle_callbacks!(integrator)
    end
  else # Non adaptive
    integrator.tprev = integrator.t
    integrator.t = integrator.t + integrator.dt
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    handle_callbacks!(integrator)
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    Juno.msg(integrator.prog,integrator.opts.progress_message(integrator.dt,integrator.t,integrator.u))
    Juno.progress(integrator.prog,integrator.t/integrator.T)
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
  if integrator.sol.t[integrator.saveiter] != integrator.t
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if integrator.opts.save_idxs == nothing
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    else
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u[integrator.opts.save_idxs],Val{false})
    end
  end
  if integrator.W.t[end] != integrator.t
    accept_step!(integrator.W,integrator.dt,false)
  end
end

@inline function postamble!(integrator)
  solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  !(typeof(integrator.prog)<:Void) && Juno.done(integrator.prog)
  return nothing
end

@inline function handle_callbacks!(integrator)
  discrete_callbacks = integrator.opts.callback.discrete_callbacks
  continuous_callbacks = integrator.opts.callback.continuous_callbacks
  atleast_one_callback = false

  continuous_modified = false
  discrete_modified = false
  if !(typeof(continuous_callbacks)<:Tuple{})
    time,upcrossing,idx,counter = find_first_continuous_callback(integrator,continuous_callbacks...)
    if time != zero(typeof(integrator.t)) && upcrossing != 0 # if not, then no events
      atleast_one_callback = true
      continuous_modified = apply_callback!(integrator,continuous_callbacks[idx],time,upcrossing)
    end
  end
  if !(typeof(discrete_callbacks)<:Tuple{})
    atleast_one_callback = true
    discrete_modified = apply_discrete_callback!(integrator,discrete_callbacks...)
  end
  if !atleast_one_callback
    update_running_noise!(integrator)
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
  integrator.sqdt = sqrt(abs(integrator.dt)) # It can change dt, like in RSwM1
end

@inline function update_running_noise!(integrator)
end

@inline function perform_rswm_rejection!(integrator)
  integrator.q = integrator.dtnew/integrator.dt
  if adaptive_alg(integrator.alg.rswm)==:RSwM1 || adaptive_alg(integrator.alg.rswm)==:RSwM2
    generate_tildes(integrator,integrator.q*integrator.ΔW,integrator.q*integrator.ΔZ,sqrt(abs((1-integrator.q)*integrator.dtnew)))
    cutLength = integrator.dt-integrator.dtnew
    if cutLength > integrator.alg.rswm.discard_length
      push!(integrator.S₁,(cutLength,integrator.ΔW-integrator.ΔWtilde,integrator.ΔZ-integrator.ΔZtilde))
    end
    if length(integrator.S₁) > integrator.sol.maxstacksize
        integrator.sol.maxstacksize = length(integrator.S₁)
    end
    if isinplace(integrator.sol.prob)
      copy!(integrator.ΔW,integrator.ΔWtilde)
      if alg_needs_extra_process(integrator.alg)
        copy!(integrator.ΔZ,integrator.ΔZtilde)
      end
    else
      integrator.ΔW = integrator.ΔWtilde
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZ = integrator.ΔZtilde
      end
    end
    integrator.dt = integrator.dtnew
    integrator.sqdt = sqrt(abs(integrator.dt))
  else # RSwM3
    if !(isinplace(integrator.sol.prob))
      dttmp = zero(integrator.t); integrator.ΔWtmp = zero(integrator.ΔWtmp)
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZtmp = zero(integrator.ΔZtmp)
      end
    else
      dttmp = zero(integrator.t); fill!(integrator.ΔWtmp,zero(eltype(integrator.ΔWtmp)))
      if alg_needs_extra_process(integrator.alg)
        fill!(integrator.ΔZtmp,zero(eltype(integrator.ΔZtmp)))
      end
    end
    if length(integrator.S₂) > integrator.sol.maxstacksize2
      integrator.sol.maxstacksize2= length(integrator.S₂)
    end
    while !isempty(integrator.S₂)
      L₁,L₂,L₃ = pop!(integrator.S₂)
      if dttmp + L₁ < (1-integrator.q)*integrator.dt #while the backwards movement is less than chop off
        dttmp += L₁
        if isinplace(integrator.sol.prob)
          #=
          integrator.ΔWtmp .+= L₂
          if alg_needs_extra_process(integrator.alg)
            integrator.ΔZtmp .+= L₃
          end
          =#
          @tight_loop_macros for i in eachindex(integrator.ΔW)
            @inbounds integrator.ΔWtmp[i] += L₂[i]
          end
          if alg_needs_extra_process(integrator.alg)
            @tight_loop_macros for i in eachindex(integrator.ΔW)
              @inbounds integrator.ΔZtmp[i] += L₃[i]
            end
          end
        else
          integrator.ΔWtmp += L₂
          if alg_needs_extra_process(integrator.alg)
            integrator.ΔZtmp += L₃
          end
        end
        push!(integrator.S₁,(L₁,L₂,L₃))
      else
        push!(integrator.S₂,(L₁,L₂,L₃))
        break
      end
    end # end while
    dtK = integrator.dt - dttmp
    qK = integrator.q*integrator.dt/dtK
    if isinplace(integrator.sol.prob)
      #@. integrator.ΔWtmp = integrator.ΔW - integrator.ΔWtmp
      @tight_loop_macros for i in eachindex(integrator.u)
        @inbounds integrator.ΔWtmp[i] = integrator.ΔW[i] - integrator.ΔWtmp[i]
      end
      if alg_needs_extra_process(integrator.alg)
        #@. integrator.ΔZtmp = integrator.ΔZ - integrator.ΔZtmp
        @tight_loop_macros for i in eachindex(integrator.u)
          @inbounds integrator.ΔZtmp[i] = integrator.ΔZ[i] - integrator.ΔZtmp[i]
        end
      end
    else
      integrator.ΔWtmp = integrator.ΔW - integrator.ΔWtmp
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZtmp = integrator.ΔZ - integrator.ΔZtmp
      end
    end
    generate_tildes(integrator,qK*integrator.ΔWtmp,qK*integrator.ΔZtmp,sqrt(abs((1-qK)*qK*dtK)))
    cutLength = (1-qK)*dtK
    if cutLength > integrator.alg.rswm.discard_length
      push!(integrator.S₁,(cutLength,integrator.ΔWtmp-integrator.ΔWtilde,integrator.ΔZtmp-integrator.ΔZtilde))
    end
    if length(integrator.S₁) > integrator.sol.maxstacksize
        integrator.sol.maxstacksize = length(integrator.S₁)
    end
    integrator.dt = integrator.dtnew
    integrator.sqdt = sqrt(abs(integrator.dt))
    if isinplace(integrator.sol.prob)
      copy!(integrator.ΔW,integrator.ΔWtilde)
      if alg_needs_extra_process(integrator.alg)
        copy!(integrator.ΔZ,integrator.ΔZtilde)
      end
    else
      integrator.ΔW = integrator.ΔWtilde
      if alg_needs_extra_process(integrator.alg)
        integrator.ΔZ = integrator.ΔZtilde
      end
    end
  end
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
    scale!(integrator.ΔW,scaling_factor)
    if alg_needs_extra_process(integrator.alg)
      integrator.noise(integrator.ΔZ,integrator)
      scale!(integrator.ΔZ,scaling_factor)
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
      #@. integrator.ΔWtilde = add1 + scaling*integrator.ΔWtilde
      @tight_loop_macros for i in eachinex(integrator.u)
        @inbounds integrator.ΔWtilde[i] = add1[i] + scaling*integrator.ΔWtilde[i]
      end
    else
      #@. integrator.ΔWtilde = scaling*integrator.ΔWtilde
      @tight_loop_macros for i in eachinex(integrator.u)
        @inbounds integrator.ΔWtilde[i] = scaling*integrator.ΔWtilde[i]
      end
    end
    if alg_needs_extra_process(integrator.alg)
      integrator.noise(integrator.ΔZtilde,integrator)
      if add2 != 0
        #@. integrator.ΔZtilde = add2 + scaling*integrator.ΔZtilde
        @tight_loop_macros for i in eachinex(integrator.u)
          @inbounds integrator.ΔZtilde[i] = add2[i] + scaling*integrator.ΔZtilde[i]
        end
      else
        #@. integrator.ΔZtilde = scaling*integrator.ΔZtilde
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
