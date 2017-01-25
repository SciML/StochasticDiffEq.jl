@def sde_preamble begin
  local T::tType
  @unpack T,rands,cache = integrator
end

@def sde_loopheader begin
  integrator.iter += 1
  if integrator.iter > integrator.opts.maxiters
    warn("Max Iters Reached. Aborting")
    @sde_postamble
  end
  if integrator.dt == 0
    warn("dt == 0. Aborting")
    @sde_postamble
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.t,integrator.u)
    warn("Instability detected. Aborting")
    @sde_postamble
  end
end

function savevalues!(integrator::SDEIntegrator)
  if integrator.opts.save_timeseries && integrator.iter%integrator.opts.timeseries_steps==0
    push!(integrator.sol.u,copy(integrator.u))
    push!(integrator.sol.t,integrator.t)
    if integrator.opts.save_noise
      push!(integrator.sol.W,copy(integrator.W))
    end
  end
end

@def sde_loopfooter begin
  if integrator.opts.adaptive
    integrator.q11 = integrator.EEst^integrator.opts.beta1
    integrator.q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    integrator.q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma))
    dtnew = integrator.dt/integrator.q
    ttmp = integrator.t + integrator.dt
    integrator.isout = integrator.opts.isoutofdomain(ttmp,integrator.u)
    integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0)
    if integrator.accept_step # Accepted
      integrator.t = ttmp
      integrator.qold = max(integrator.EEst,integrator.opts.qoldinit)
      if integrator.tdir > 0
        integrator.dtpropose = min(integrator.opts.dtmax,dtnew)
      else
        integrator.integrator.dtpropose = max(integrator.opts.dtmax,dtnew)
      end
      if integrator.tdir > 0
        integrator.dtpropose = max(integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
      else
        integrator.integrator.dtpropose = min(integrator.integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
      end

      update_running_noise!(integrator)
      if uType <: AbstractArray
        recursivecopy!(integrator.uprev,integrator.u)
      else
        integrator.uprev = integrator.u
      end
      savevalues!(integrator)
      # Setup next step
      if adaptive_alg(integrator.alg.rswm)==:RSwM3
        ResettableStacks.reset!(integrator.S₂) #Empty integrator.S₂
      end
      if adaptive_alg(integrator.alg.rswm)==:RSwM1
        if !isempty(integrator.S₁)
          integrator.dt,integrator.ΔW,integrator.ΔZ = pop!(integrator.S₁)
          integrator.sqdt = sqrt(integrator.dt)
        else # Stack is empty
          c = min(integrator.opts.dtmax,dtnew)
          integrator.dt = max(min(c,abs(T-integrator.t)),integrator.opts.dtmin)#abs to fix complex sqrt issue at end
          #integrator.dt = min(c,abs(T-integrator.t))
          integrator.sqdt = sqrt(integrator.dt)
          integrator.ΔW = integrator.sqdt*next(rands)
          integrator.ΔZ = integrator.sqdt*next(rands)
        end
      elseif adaptive_alg(integrator.alg.rswm)==:RSwM2 || adaptive_alg(integrator.alg.rswm)==:RSwM3
        c = min(integrator.opts.dtmax,dtnew)
        integrator.dt = max(min(c,abs(T-integrator.t)),integrator.opts.dtmin) #abs to fix complex sqrt issue at end
        integrator.sqdt = sqrt(integrator.dt)
        if !(uType <: AbstractArray)
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
            ΔWtilde = qtmp*L₂ + sqrt((1-qtmp)*qtmp*L₁)*next(rands)
            ΔZtilde = qtmp*L₃ + sqrt((1-qtmp)*qtmp*L₁)*next(rands)
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
          ΔWtilde = sqrt(dtleft)*next(rands)
          ΔZtilde = sqrt(dtleft)*next(rands)
          integrator.ΔW += ΔWtilde
          integrator.ΔZ += ΔZtilde
          if adaptive_alg(integrator.alg.rswm)==:RSwM3
            push!(integrator.S₂,(dtleft,ΔWtilde,ΔZtilde))
          end
        end
      end # End RSwM2 and RSwM3
    else #Rejection
      dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
      integrator.q = dtnew/integrator.dt
      if adaptive_alg(integrator.alg.rswm)==:RSwM1 || adaptive_alg(integrator.alg.rswm)==:RSwM2
        ΔWtmp = integrator.q*integrator.ΔW + sqrt((1-integrator.q)*dtnew)*next(rands)
        ΔZtmp = integrator.q*integrator.ΔZ + sqrt((1-integrator.q)*dtnew)*next(rands)
        cutLength = integrator.dt-dtnew
        if cutLength > integrator.alg.rswm.discard_length
          push!(integrator.S₁,(cutLength,integrator.ΔW-ΔWtmp,integrator.ΔZ-ΔZtmp))
        end
        if length(integrator.S₁) > integrator.sol.maxstacksize
            integrator.sol.maxstacksize = length(integrator.S₁)
        end
        integrator.ΔW = ΔWtmp
        integrator.ΔZ = ΔZtmp
        integrator.dt = dtnew
      else # RSwM3
        if !(uType <: AbstractArray)
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
        ΔWtilde = qK*K₂ + sqrt((1-qK)*qK*dtK)*next(rands)
        ΔZtilde = qK*K₃ + sqrt((1-qK)*qK*dtK)*next(rands)
        cutLength = (1-qK)*dtK
        if cutLength > integrator.alg.rswm.discard_length
          push!(integrator.S₁,(cutLength,K₂-ΔWtilde,K₃-ΔZtilde))
        end
        if length(integrator.S₁) > integrator.sol.maxstacksize
            integrator.sol.maxstacksize = length(integrator.S₁)
        end
        integrator.dt = dtnew
        integrator.ΔW = ΔWtilde
        integrator.ΔZ = ΔZtilde
      end
    end
  else # Non adaptive
    integrator.t = integrator.t + integrator.dt
    integrator.accept_step = true
    integrator.dtpropose = integrator.dt

    if typeof(integrator.u) <: AbstractArray
      recursivecopy!(integrator.uprev,integrator.u)
    else
      integrator.uprev = integrator.u
    end
    update_running_noise!(integrator)
    integrator.ΔW = integrator.sqdt*next(integrator.rands)
    if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
      integrator.ΔZ = integrator.sqdt*next(integrator.rands)
    end
    savevalues!(integrator)
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    Juno.msg(integrator.prog,integrator.opts.progress_message(integrator.dt,integrator.t,integrator.u))
    Juno.progress(integrator.prog,integrator.t/T)
  end
end

@def sde_postamble begin
  if integrator.sol.t[end] != integrator.t
    push!(integrator.sol.t,integrator.t)
    push!(integrator.sol.u,integrator.u)
    if integrator.opts.save_noise
      push!(integrator.sol.W,integrator.W)
    end
  end
  integrator.opts.progress && Juno.done(integrator.prog)
  return nothing
end

function update_running_noise!(integrator)
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
