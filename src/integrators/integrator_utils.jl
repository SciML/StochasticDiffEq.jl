@inline function loopheader!(integrator::SDEIntegrator)
  # Apply right after iterators / callbacks

  # Accept or reject the step
  if integrator.iter > 0
    if (integrator.opts.adaptive && integrator.accept_step) || !integrator.opts.adaptive
      apply_step!(integrator)
    elseif integrator.opts.adaptive && !integrator.accept_step
      perform_rswm_rejection!(integrator)
    end
  end

  integrator.iter += 1
  choose_algorithm!(integrator,integrator.cache)
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
  integrator.sqdt = sqrt(abs(integrator.dt))
end

@inline choose_algorithm!(integrator,cache::StochasticDiffEqCache) = nothing

@def sde_exit_condtions begin
  if integrator.iter > integrator.opts.maxiters
    if integrator.opts.verbose
      warn("Max Iters Reached. Aborting")
    end
    postamble!(integrator)
    return nothing
  end
  if integrator.dt == 0
    if integrator.opts.verbose
      warn("dt == 0. Aborting")
    end
    postamble!(integrator)
    return nothing
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.t,integrator.u)
    if integrator.opts.verbose
      warn("Instability detected. Aborting")
    end
    postamble!(integrator)
    return nothing
  end
end

@inline function savevalues!(integrator::SDEIntegrator)
  while !isempty(integrator.opts.saveat) && integrator.tdir*top(integrator.opts.saveat) <= integrator.tdir*integrator.t # Perform saveat
    integrator.saveiter += 1
    curt = pop!(integrator.opts.saveat)
    if integrator.opts.saveat!=integrator.t # If <t, interpolate
      Θ = (curt - integrator.tprev)/integrator.dt
      val = sde_interpolant(Θ,integrator)
      copyat_or_push!(integrator.sol.t,integrator.saveiter,curt)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,val)
      if typeof(integrator.alg) <: StochasticEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
      if integrator.opts.save_noise
        copyat_or_push!(integrator.sol.W,integrator.saveiter,integrator.W)
      end
    else # ==t, just save
      copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
      copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
      if typeof(alg) <: StochasticDiffEqCompositeAlgorithm
        copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
      end
      if integrator.opts.save_noise
        copyat_or_push!(integrator.sol.W,integrator.saveiter,integrator.W)
      end
    end
  end
  if integrator.opts.save_timeseries && integrator.iter%integrator.opts.timeseries_steps==0
    integrator.saveiter += 1
    copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    copyat_or_push!(integrator.sol.t,integrator.saveiter,integrator.t)
    if typeof(integrator.alg) <: StochasticDiffEqCompositeAlgorithm
      copyat_or_push!(integrator.sol.alg_choice,integrator.saveiter,integrator.cache.current)
    end
    if integrator.opts.save_noise
      copyat_or_push!(integrator.sol.W,integrator.saveiter,integrator.W)
    end
  end
end

function loopfooter!(integrator::SDEIntegrator)
  if integrator.opts.adaptive
    integrator.q11 = integrator.EEst^integrator.opts.beta1
    integrator.q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    integrator.q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma))
    integrator.dtnew = integrator.dt/integrator.q
    ttmp = integrator.t + integrator.dt
    integrator.isout = integrator.opts.isoutofdomain(ttmp,integrator.u)
    integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0)
    if integrator.accept_step # Accepted
      integrator.t = ttmp
      calc_dt_propose!(integrator)
      update_running_noise!(integrator)
      savevalues!(integrator)
    end
  else # Non adaptive
    integrator.t = integrator.t + integrator.dt
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt
    update_running_noise!(integrator)
    savevalues!(integrator)
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
    copyat_or_push!(integrator.sol.u,integrator.saveiter,integrator.u)
    if integrator.opts.save_noise
      copyat_or_push!(integrator.sol.W,integrator.saveiter,integrator.W)
    end
  end
end

function postamble!(integrator)
  solution_endpoint_match_cur_integrator!(integrator)
  resize!(integrator.sol.t,integrator.saveiter)
  resize!(integrator.sol.u,integrator.saveiter)
  if integrator.opts.save_noise
    resize!(integrator.sol.W,integrator.saveiter)
  end
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
  if typeof(integrator.u) <: AbstractArray
    recursivecopy!(integrator.uprev,integrator.u)
  else
    integrator.uprev = integrator.u
  end
  # Setup next step
  if integrator.opts.adaptive
    if adaptive_alg(integrator.alg.rswm)==:RSwM3
      ResettableStacks.reset!(integrator.S₂) #Empty integrator.S₂
    end
    if adaptive_alg(integrator.alg.rswm)==:RSwM1
      if !isempty(integrator.S₁)
        integrator.dt,integrator.ΔW,integrator.ΔZ = pop!(integrator.S₁)
        integrator.sqdt = sqrt(abs(integrator.dt))
      else # Stack is empty
        integrator.dt = integrator.dtpropose
        modify_dt_for_tstops!(integrator)
        update_noise!(integrator)
      end
    elseif adaptive_alg(integrator.alg.rswm)==:RSwM2 || adaptive_alg(integrator.alg.rswm)==:RSwM3
      integrator.dt = integrator.dtpropose
      modify_dt_for_tstops!(integrator)
      if !(typeof(integrator.u) <: AbstractArray)
        dttmp = 0.0; integrator.ΔW = 0.0; integrator.ΔZ = 0.0
      else
        dttmp = 0.0; integrator.ΔW = zeros(size(integrator.u)...); integrator.ΔZ = zeros(size(integrator.u)...)
      end
      while !isempty(integrator.S₁)
        L₁,L₂,L₃ = pop!(integrator.S₁)
        qtmp = (integrator.dt-dttmp)/L₁
        if qtmp>1
          dttmp+=L₁
          integrator.ΔW+=L₂
          integrator.ΔZ+=L₃
          if adaptive_alg(integrator.alg.rswm)==:RSwM3
            push!(integrator.S₂,(L₁,L₂,L₃))
          end
        else #Popped too far
          ΔWtilde,ΔZtilde = generate_tildes(integrator,qtmp*L₂,qtmp*L₃,sqrt(abs((1-qtmp)*qtmp*L₁)))
          integrator.ΔW += ΔWtilde
          integrator.ΔZ += ΔZtilde
          if (1-qtmp)*L₁ > integrator.alg.rswm.discard_length
            push!(integrator.S₁,((1-qtmp)*L₁,L₂-ΔWtilde,L₃-ΔZtilde))
            if adaptive_alg(integrator.alg.rswm)==:RSwM3 && qtmp*L₁ > integrator.alg.rswm.discard_length
              push!(integrator.S₂,(qtmp*L₁,ΔWtilde,ΔZtilde))
            end
          end
          break
        end
      end #end while empty
      dtleft = integrator.dt - dttmp
      if dtleft != 0 #Stack emptied
        ΔWtilde,ΔZtilde = generate_tildes(integrator,0,0,sqrt(abs(dtleft)))
        integrator.ΔW += ΔWtilde
        integrator.ΔZ += ΔZtilde
        if adaptive_alg(integrator.alg.rswm)==:RSwM3
          push!(integrator.S₂,(dtleft,ΔWtilde,ΔZtilde))
        end
      end
    end # End RSwM2 and RSwM3
  else # Not adaptive
    modify_dt_for_tstops!(integrator)
    update_noise!(integrator)
  end
end

@inline function update_running_noise!(integrator)
  if integrator.opts.save_noise
    if typeof(integrator.u) <: AbstractArray
      for i in eachindex(integrator.u)
        integrator.W[i] = integrator.W[i] + integrator.ΔW[i]
      end
    else
      integrator.W = integrator.W + integrator.ΔW
    end
    if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
      if typeof(integrator.u) <: AbstractArray
        for i in eachindex(integrator.u)
          integrator.Z[i] = integrator.Z[i] + integrator.ΔZ[i]
        end
      else
        integrator.Z = integrator.Z + integrator.ΔZ
      end
    end
  end
end

@inline function perform_rswm_rejection!(integrator)
  if integrator.isout
    integrator.dtnew = integrator.dt*integrator.opts.qmin
  else
    integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
  end
  integrator.q = integrator.dtnew/integrator.dt
  if adaptive_alg(integrator.alg.rswm)==:RSwM1 || adaptive_alg(integrator.alg.rswm)==:RSwM2
    ΔWtmp,ΔZtmp = generate_tildes(integrator,integrator.q*integrator.ΔW,integrator.q*integrator.ΔZ,sqrt(abs((1-integrator.q)*integrator.dtnew)))
    cutLength = integrator.dt-integrator.dtnew
    if cutLength > integrator.alg.rswm.discard_length
      push!(integrator.S₁,(cutLength,integrator.ΔW-ΔWtmp,integrator.ΔZ-ΔZtmp))
    end
    if length(integrator.S₁) > integrator.sol.maxstacksize
        integrator.sol.maxstacksize = length(integrator.S₁)
    end
    integrator.ΔW = ΔWtmp
    integrator.ΔZ = ΔZtmp
    integrator.dt = integrator.dtnew
    integrator.sqdt = sqrt(integrator.dt)
  else # RSwM3
    if !(typeof(integrator.u) <: AbstractArray)
      dttmp = 0.0; ΔWtmp = 0.0; ΔZtmp = 0.0
    else
      dttmp = 0.0; ΔWtmp = zeros(size(integrator.u)...); ΔZtmp = zeros(size(integrator.u)...)
    end
    if length(integrator.S₂) > integrator.sol.maxstacksize2
      integrator.sol.maxstacksize2= length(integrator.S₂)
    end
    while !isempty(integrator.S₂)
      L₁,L₂,L₃ = pop!(integrator.S₂)
      if dttmp + L₁ < (1-integrator.q)*integrator.dt #while the backwards movement is less than chop off
        dttmp += L₁
        ΔWtmp += L₂
        ΔZtmp += L₃
        push!(integrator.S₁,(L₁,L₂,L₃))
      else
        push!(integrator.S₂,(L₁,L₂,L₃))
        break
      end
    end # end while
    dtK = integrator.dt - dttmp
    K₂ = integrator.ΔW - ΔWtmp
    K₃ = integrator.ΔZ - ΔZtmp
    qK = integrator.q*integrator.dt/dtK
    ΔWtilde,ΔZtilde = generate_tildes(integrator,qK*K₂,qK*K₃,sqrt(abs((1-qK)*qK*dtK)))
    cutLength = (1-qK)*dtK
    if cutLength > integrator.alg.rswm.discard_length
      push!(integrator.S₁,(cutLength,K₂-ΔWtilde,K₃-ΔZtilde))
    end
    if length(integrator.S₁) > integrator.sol.maxstacksize
        integrator.sol.maxstacksize = length(integrator.S₁)
    end
    integrator.dt = integrator.dtnew
    integrator.sqdt = sqrt(integrator.dt)
    integrator.ΔW = ΔWtilde
    integrator.ΔZ = ΔZtilde
    integrator.sqdt = sqrt(abs(integrator.dt))
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
    elseif t > ts_top
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
    integrator.noise(integrator.ΔW)
    integrator.ΔW .*= scaling_factor
    if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
      integrator.noise(integrator.ΔZ)
      integrator.ΔZ .*= scaling_factor
    end
  else
    if (typeof(integrator.u) <: AbstractArray)
      integrator.ΔW = scaling_factor*integrator.noise(size(integrator.u))
      if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
        integrator.ΔZ = scaling_factor*integrator.noise(size(integrator.u))
      end
    else
      integrator.ΔW = scaling_factor*integrator.noise()
      if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
        integrator.ΔZ = scaling_factor*integrator.noise()
      end
    end
  end
end

@inline function generate_tildes(integrator,add1,add2,scaling)
  local ΔWtilde::typeof(integrator.ΔW)
  local ΔZtilde::typeof(integrator.ΔZ)
  if isinplace(integrator.noise)
    ΔWtilde = similar(integrator.ΔW)
    integrator.noise(ΔWtilde)
    if add1 != 0
      ΔWtilde .= add1 .+ scaling.*ΔWtilde
    else
      ΔWtilde .= scaling.*ΔWtilde
    end
    if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
      ΔZtilde = similar(integrator.ΔZ)
      integrator.noise(ΔZtilde)
      if add2 != 0
        ΔZtilde .= add2 .+ scaling.*ΔZtilde
      else
        ΔZtilde .= scaling.*ΔZtilde
      end
    end
  else
    if (typeof(integrator.u) <: AbstractArray)
      ΔWtilde = add1 .+ scaling.*integrator.noise(size(integrator.u))
      if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
        ΔZtilde = add2 .+ scaling.*integrator.noise(size(integrator.u))
      end
    else
      ΔWtilde = add1 + scaling*integrator.noise()
      if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
        ΔZtilde = add2 + scaling*integrator.noise()
      end
    end
  end
  ΔWtilde,ΔZtilde
end
