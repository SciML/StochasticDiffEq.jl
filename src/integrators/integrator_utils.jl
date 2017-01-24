@def sde_preamble begin
  local T::tType
  local ΔW::randType
  local ΔZ::randType
  @unpack T,rands,W,Z,cache = integrator

  if uType <: AbstractArray
    EEsttmp = zeros(integrator.u)
  end
  ΔW = integrator.sqdt*next(rands) # Take one first
  ΔZ = integrator.sqdt*next(rands) # Take one first
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

@def sde_savevalues begin
  if integrator.opts.save_timeseries && integrator.iter%integrator.opts.timeseries_steps==0
    push!(integrator.sol.u,copy(integrator.u))
    push!(integrator.sol.t,integrator.t)
    push!(integrator.sol.W,copy(W))
  end
end

@def sde_loopfooter begin
  if integrator.opts.adaptive
    integrator.q11 = EEst^integrator.opts.beta1
    q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
    dtnew = integrator.dt/q
    ttmp = integrator.t + integrator.dt
    #integrator.isout = integrator.opts.isoutofdomain(ttmp,integrator.u)
    #integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0)
    if EEst <= 1 # Accepted
      acceptedIters += 1
      integrator.t = ttmp
      integrator.qold = max(EEst,integrator.opts.qoldinit)
      #if integrator.tdir > 0
        dtpropose = min(integrator.opts.dtmax,dtnew)
      #else
      #  integrator.dtpropose = max(integrator.opts.dtmax,dtnew)
      #end
      #if integrator.tdir > 0
        dtpropose = max(dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
      #else
      #  integrator.dtpropose = min(integrator.dtpropose,integrator.opts.dtmin) #abs to fix complex sqrt issue at end
      #end


      if uType <: AbstractArray
        for i in eachindex(integrator.u)
          W[i] = W[i] + ΔW[i]
          Z[i] = Z[i] + ΔZ[i]
        end
      else
        W = W + ΔW
        Z = Z + ΔZ
      end
      if uType <: AbstractArray
        recursivecopy!(integrator.uprev,integrator.u)
      else
        integrator.uprev = integrator.u
      end
      if adaptive_alg(integrator.alg.rswm)==:RSwM3
        ResettableStacks.reset!(S₂) #Empty S₂
      end
      @sde_savevalues
      # Setup next step
      if adaptive_alg(integrator.alg.rswm)==:RSwM1
        if !isempty(S₁)
          integrator.dt,ΔW,ΔZ = pop!(S₁)
          integrator.sqdt = sqrt(integrator.dt)
        else # Stack is empty
          c = min(integrator.opts.dtmax,dtnew)
          integrator.dt = max(min(c,abs(T-integrator.t)),integrator.opts.dtmin)#abs to fix complex sqrt issue at end
          #integrator.dt = min(c,abs(T-integrator.t))
          integrator.sqdt = sqrt(integrator.dt)
          ΔW = integrator.sqdt*next(rands)
          ΔZ = integrator.sqdt*next(rands)
        end
      elseif adaptive_alg(integrator.alg.rswm)==:RSwM2 || adaptive_alg(integrator.alg.rswm)==:RSwM3
        c = min(integrator.opts.dtmax,dtnew)
        integrator.dt = max(min(c,abs(T-integrator.t)),integrator.opts.dtmin) #abs to fix complex sqrt issue at end
        integrator.sqdt = sqrt(integrator.dt)
        if !(uType <: AbstractArray)
          dttmp = 0.0; ΔW = 0.0; ΔZ = 0.0
        else
          dttmp = 0.0; ΔW = zeros(size(integrator.u)...); ΔZ = zeros(size(integrator.u)...)
        end
        while !isempty(S₁)
          L₁,L₂,L₃ = pop!(S₁)
          qtmp = (integrator.dt-dttmp)/L₁
          if qtmp>1
            dttmp+=L₁
            ΔW+=L₂
            ΔZ+=L₃
            if adaptive_alg(integrator.alg.rswm)==:RSwM3
              push!(S₂,(L₁,L₂,L₃))
            end
          else #Popped too far
            ΔWtilde = qtmp*L₂ + sqrt((1-qtmp)*qtmp*L₁)*next(rands)
            ΔZtilde = qtmp*L₃ + sqrt((1-qtmp)*qtmp*L₁)*next(rands)
            ΔW += ΔWtilde
            ΔZ += ΔZtilde
            if (1-qtmp)*L₁ > integrator.alg.rswm.discard_length
              push!(S₁,((1-qtmp)*L₁,L₂-ΔWtilde,L₃-ΔZtilde))
              if adaptive_alg(integrator.alg.rswm)==:RSwM3 && qtmp*L₁ > integrator.alg.rswm.discard_length
                push!(S₂,(qtmp*L₁,ΔWtilde,ΔZtilde))
              end
            end
            break
          end
        end #end while empty
        dtleft = integrator.dt - dttmp
        if dtleft != 0 #Stack emptied
          ΔWtilde = sqrt(dtleft)*next(rands)
          ΔZtilde = sqrt(dtleft)*next(rands)
          ΔW += ΔWtilde
          ΔZ += ΔZtilde
          if adaptive_alg(integrator.alg.rswm)==:RSwM3
            push!(S₂,(dtleft,ΔWtilde,ΔZtilde))
          end
        end
      end # End RSwM2 and RSwM3
    else #Rejection
      dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
      q = dtnew/integrator.dt
      if adaptive_alg(integrator.alg.rswm)==:RSwM1 || adaptive_alg(integrator.alg.rswm)==:RSwM2
        ΔWtmp = q*ΔW + sqrt((1-q)*dtnew)*next(rands)
        ΔZtmp = q*ΔZ + sqrt((1-q)*dtnew)*next(rands)
        cutLength = integrator.dt-dtnew
        if cutLength > integrator.alg.rswm.discard_length
          push!(S₁,(cutLength,ΔW-ΔWtmp,ΔZ-ΔZtmp))
        end
        if length(S₁) > integrator.sol.maxstacksize
            integrator.sol.maxstacksize = length(S₁)
        end
        ΔW = ΔWtmp
        ΔZ = ΔZtmp
        integrator.dt = dtnew
      else # RSwM3
        if !(uType <: AbstractArray)
          dttmp = 0.0; ΔWtmp = 0.0; ΔZtmp = 0.0
        else
          dttmp = 0.0; ΔWtmp = zeros(size(integrator.u)...); ΔZtmp = zeros(size(integrator.u)...)
        end
        if length(S₂) > integrator.sol.maxstacksize2
          integrator.sol.maxstacksize2= length(S₂)
        end
        while !isempty(S₂)
          L₁,L₂,L₃ = pop!(S₂)
          if dttmp + L₁ < (1-q)*integrator.dt #while the backwards movement is less than chop off
            dttmp += L₁
            ΔWtmp += L₂
            ΔZtmp += L₃
            push!(S₁,(L₁,L₂,L₃))
          else
            push!(S₂,(L₁,L₂,L₃))
            break
          end
        end # end while
        dtK = integrator.dt - dttmp
        K₂ = ΔW - ΔWtmp
        K₃ = ΔZ - ΔZtmp
        qK = q*integrator.dt/dtK
        ΔWtilde = qK*K₂ + sqrt((1-qK)*qK*dtK)*next(rands)
        ΔZtilde = qK*K₃ + sqrt((1-qK)*qK*dtK)*next(rands)
        cutLength = (1-qK)*dtK
        if cutLength > integrator.alg.rswm.discard_length
          push!(S₁,(cutLength,K₂-ΔWtilde,K₃-ΔZtilde))
        end
        if length(S₁) > integrator.sol.maxstacksize
            integrator.sol.maxstacksize = length(S₁)
        end
        integrator.dt = dtnew
        ΔW = ΔWtilde
        ΔZ = ΔZtilde
      end
    end
  else # Non adaptive
    integrator.t = integrator.t + integrator.dt

    if typeof(integrator.u) <: AbstractArray
      recursivecopy!(integrator.uprev,integrator.u)
    else
      integrator.uprev = integrator.u
    end

    if uType <: AbstractArray
      for i in eachindex(integrator.u)
        W[i] = W[i] + ΔW[i]
      end
    else
      W = W + ΔW
    end
    ΔW = integrator.sqdt*next(rands)
    if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
      if uType <: AbstractArray
        for i in eachindex(integrator.u)
          Z[i] = Z[i] + ΔZ[i]
        end
      else
        Z = Z + ΔZ
      end
      ΔZ = integrator.sqdt*next(rands)
    end
    @sde_savevalues
  end
  if integrator.opts.progress && integrator.iter%integrator.opts.progress_steps==0
    Juno.msg(integrator.prog,integrator.opts.progress_message(integrator.dt,integrator.t,integrator.u))
    Juno.progress(integrator.prog,integrator.t/T)
  end
end

@def sde_adaptiveprelim begin
  if integrator.opts.adaptive
    S₁ = DataStructures.Stack{}(Tuple{typeof(integrator.t),typeof(W),typeof(Z)})
    acceptedIters = 0
    if adaptive_alg(integrator.alg.rswm)==:RSwM3
      S₂ = ResettableStacks.ResettableStack{}(Tuple{typeof(integrator.t),typeof(W),typeof(Z)})
    end
  end
end

@def sde_postamble begin
  if integrator.sol.t[end] != integrator.t
    push!(integrator.sol.t,integrator.t)
    push!(integrator.sol.u,integrator.u)
    push!(integrator.sol.W,W)
  end
  integrator.opts.progress && Juno.done(integrator.prog)
  return nothing
end
