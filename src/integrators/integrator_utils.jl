type SDEIntegrator{T1,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5,OType}
  f::F4
  g::F5
  u::uType
  t::tType
  dt::tType
  T::tType
  alg::T1
  timeseries::Vector{uType}
  Ws::Vector{randType}
  ts::Vector{tType}
  adaptivealg::Symbol
  δ::uEltypeNoUnits
  internalnorm::F
  discard_length::tType
  unstable_check::F3
  rands::ChunkedArray{uEltypeNoUnits,Nm1,N}
  sqdt::tType
  W::randType
  Z::randType
  opts::OType
  qold::tTypeNoUnits
  q11::tTypeNoUnits
end

@def sde_preamble begin
  local u::uType
  local t::tType
  local dt::tType
  local T::tType
  local ΔW::randType
  local ΔZ::randType
  @unpack u,t,dt,T,timeseries,Ws,ts,adaptivealg,δ,discard_length,unstable_check,rands,sqdt,W,Z = integrator

  integrator.opts.progress && (prog = Juno.ProgressBar(name=integrator.opts.progress_name))
  if uType <: AbstractArray
    EEsttmp = zeros(u)
  end
  if uType <: AbstractArray
    utmp = zeros(u)
  end
  iter = 0
  max_stack_size = 0
  max_stack_size2 = 0
  ΔW = sqdt*next(rands) # Take one first
  ΔZ = sqdt*next(rands) # Take one first
end

@def sde_sritableaupreamble begin
  local c₀::Vector{uEltypeNoUnits}
  local c₁::Vector{uEltypeNoUnits}
  local A₀::Matrix{uEltypeNoUnits}
  local A₁::Matrix{uEltypeNoUnits}
  local B₀::Matrix{uEltypeNoUnits}
  local B₁::Matrix{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local β₁::Vector{uEltypeNoUnits}
  local β₂::Vector{uEltypeNoUnits}
  local β₃::Vector{uEltypeNoUnits}
  local β₄::Vector{uEltypeNoUnits}
end

@def sde_sratableaupreamble begin
  local c₀::Vector{uEltypeNoUnits}
  local c₁::Vector{uEltypeNoUnits}
  local A₀::Matrix{uEltypeNoUnits}
  local B₀::Matrix{uEltypeNoUnits}
  local α::Vector{uEltypeNoUnits}
  local β₁::Vector{uEltypeNoUnits}
  local β₂::Vector{uEltypeNoUnits}
end

@def sde_loopheader begin
  iter += 1
  if iter > integrator.opts.maxiters
    warn("Max Iters Reached. Aborting")
    @sde_postamble
  end
  if dt == 0
    warn("dt == 0. Aborting")
    @sde_postamble
  end
  if unstable_check(dt,t,u)
    warn("Instability detected. Aborting")
    @sde_postamble
  end
end

@def sde_savevalues begin
  if integrator.opts.save_timeseries && iter%integrator.opts.timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
    push!(Ws,copy(W))
  end
end

@def sde_loopfooter begin
  if integrator.opts.adaptive
    integrator.q11 = EEst^integrator.opts.beta1
    q = integrator.q11/(integrator.qold^integrator.opts.beta2)
    q = max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),q/integrator.opts.gamma))
    dtnew = dt/q
    ttmp = t + dt
    #integrator.isout = integrator.opts.isoutofdomain(ttmp,integrator.u)
    #integrator.accept_step = (!integrator.isout && integrator.EEst <= 1.0)
    if EEst <= 1 # Accepted
      acceptedIters += 1
      t = ttmp
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
        for i in eachindex(u)
          W[i] = W[i] + ΔW[i]
          Z[i] = Z[i] + ΔZ[i]
        end
      else
        W = W + ΔW
        Z = Z + ΔZ
      end
      if uType <: AbstractArray
        recursivecopy!(u,utmp)
      else
        u = utmp
      end
      if adaptivealg==:RSwM3
        ResettableStacks.reset!(S₂) #Empty S₂
      end
      @sde_savevalues
      # Setup next step
      if adaptivealg==:RSwM1
        if !isempty(S₁)
          dt,ΔW,ΔZ = pop!(S₁)
          sqdt = sqrt(dt)
        else # Stack is empty
          c = min(integrator.opts.dtmax,dtnew)
          dt = max(min(c,abs(T-t)),integrator.opts.dtmin)#abs to fix complex sqrt issue at end
          #dt = min(c,abs(T-t))
          sqdt = sqrt(dt)
          ΔW = sqdt*next(rands)
          ΔZ = sqdt*next(rands)
        end
      elseif adaptivealg==:RSwM2 || adaptivealg==:RSwM3
        c = min(integrator.opts.dtmax,dtnew)
        dt = max(min(c,abs(T-t)),integrator.opts.dtmin) #abs to fix complex sqrt issue at end
        sqdt = sqrt(dt)
        if !(uType <: AbstractArray)
          dttmp = 0.0; ΔW = 0.0; ΔZ = 0.0
        else
          dttmp = 0.0; ΔW = zeros(size(u)...); ΔZ = zeros(size(u)...)
        end
        while !isempty(S₁)
          L₁,L₂,L₃ = pop!(S₁)
          qtmp = (dt-dttmp)/L₁
          if qtmp>1
            dttmp+=L₁
            ΔW+=L₂
            ΔZ+=L₃
            if adaptivealg==:RSwM3
              push!(S₂,(L₁,L₂,L₃))
            end
          else #Popped too far
            ΔWtilde = qtmp*L₂ + sqrt((1-qtmp)*qtmp*L₁)*next(rands)
            ΔZtilde = qtmp*L₃ + sqrt((1-qtmp)*qtmp*L₁)*next(rands)
            ΔW += ΔWtilde
            ΔZ += ΔZtilde
            if (1-qtmp)*L₁ > discard_length
              push!(S₁,((1-qtmp)*L₁,L₂-ΔWtilde,L₃-ΔZtilde))
              if adaptivealg==:RSwM3 && qtmp*L₁ > discard_length
                push!(S₂,(qtmp*L₁,ΔWtilde,ΔZtilde))
              end
            end
            break
          end
        end #end while empty
        dtleft = dt - dttmp
        if dtleft != 0 #Stack emptied
          ΔWtilde = sqrt(dtleft)*next(rands)
          ΔZtilde = sqrt(dtleft)*next(rands)
          ΔW += ΔWtilde
          ΔZ += ΔZtilde
          if adaptivealg==:RSwM3
            push!(S₂,(dtleft,ΔWtilde,ΔZtilde))
          end
        end
      end # End RSwM2 and RSwM3
    else #Rejection
      dtnew = dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
      q = dtnew/dt
      if adaptivealg==:RSwM1 || adaptivealg==:RSwM2
        ΔWtmp = q*ΔW + sqrt((1-q)*dtnew)*next(rands)
        ΔZtmp = q*ΔZ + sqrt((1-q)*dtnew)*next(rands)
        cutLength = dt-dtnew
        if cutLength > discard_length
          push!(S₁,(cutLength,ΔW-ΔWtmp,ΔZ-ΔZtmp))
        end
        if length(S₁) > max_stack_size
            max_stack_size = length(S₁)
        end
        ΔW = ΔWtmp
        ΔZ = ΔZtmp
        dt = dtnew
      else # RSwM3
        if !(uType <: AbstractArray)
          dttmp = 0.0; ΔWtmp = 0.0; ΔZtmp = 0.0
        else
          dttmp = 0.0; ΔWtmp = zeros(size(u)...); ΔZtmp = zeros(size(u)...)
        end
        if length(S₂) > max_stack_size2
          max_stack_size2= length(S₂)
        end
        while !isempty(S₂)
          L₁,L₂,L₃ = pop!(S₂)
          if dttmp + L₁ < (1-q)*dt #while the backwards movement is less than chop off
            dttmp += L₁
            ΔWtmp += L₂
            ΔZtmp += L₃
            push!(S₁,(L₁,L₂,L₃))
          else
            push!(S₂,(L₁,L₂,L₃))
            break
          end
        end # end while
        dtK = dt - dttmp
        K₂ = ΔW - ΔWtmp
        K₃ = ΔZ - ΔZtmp
        qK = q*dt/dtK
        ΔWtilde = qK*K₂ + sqrt((1-qK)*qK*dtK)*next(rands)
        ΔZtilde = qK*K₃ + sqrt((1-qK)*qK*dtK)*next(rands)
        cutLength = (1-qK)*dtK
        if cutLength > discard_length
          push!(S₁,(cutLength,K₂-ΔWtilde,K₃-ΔZtilde))
        end
        if length(S₁) > max_stack_size
            max_stack_size = length(S₁)
        end
        dt = dtnew
        ΔW = ΔWtilde
        ΔZ = ΔZtilde
      end
    end
  else # Non adaptive
    t = t + dt
    if uType <: AbstractArray
      for i in eachindex(u)
        W[i] = W[i] + ΔW[i]
      end
    else
      W = W + ΔW
    end
    ΔW = sqdt*next(rands)
    if !(typeof(integrator.alg) <: EM) || !(typeof(integrator.alg) <: RKMil)
      if uType <: AbstractArray
        for i in eachindex(u)
          Z[i] = Z[i] + ΔZ[i]
        end
      else
        Z = Z + ΔZ
      end
      ΔZ = sqdt*next(rands)
    end
    @sde_savevalues
  end
  if integrator.opts.progress && iter%integrator.opts.progress_steps==0
    Juno.msg(prog,integrator.opts.progress_message(dt,t,u))
    Juno.progress(prog,t/T)
  end
end

@def sde_adaptiveprelim begin
  max_stack_size = 0
  max_stack_size2 = 0
  if integrator.opts.adaptive
    S₁ = DataStructures.Stack{}(Tuple{typeof(t),typeof(W),typeof(Z)})
    acceptedIters = 0
    if adaptivealg==:RSwM3
      S₂ = ResettableStacks.ResettableStack{}(Tuple{typeof(t),typeof(W),typeof(Z)})
    end
  end
end

@def sde_postamble begin
  if ts[end] != t
    push!(ts,t)
    push!(timeseries,u)
    push!(Ws,W)
  end
  integrator.opts.progress && Juno.done(prog)
  return u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end
