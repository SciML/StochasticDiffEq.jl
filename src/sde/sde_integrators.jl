immutable SDEIntegrator{T1,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}
  f::F4
  g::F5
  u::uType
  t::tType
  dt::tType
  T::tType
  alg::T1
  maxiters::Int
  timeseries::Vector{uType}
  Ws::Vector{randType}
  ts::Vector{tType}
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  adaptivealg::Symbol
  δ::uEltypeNoUnits
  γ::uEltypeNoUnits
  abstol::uEltype
  reltol::uEltypeNoUnits
  qmax::uEltypeNoUnits
  dtmax::tType
  dtmin::tType
  internalnorm::F
  discard_length::tType
  progressbar::Bool
  progressbar_name::String
  progress_steps::Int
  progress_message::F2
  unstable_check::F3
  rands::ChunkedArray{uEltypeNoUnits,Nm1,N}
  sqdt::tType
  W::randType
  Z::randType
  tableau::tableauType
end

@def sde_preamble begin
  local u::uType
  local t::tType
  local dt::tType
  local T::tType
  local ΔW::randType
  local ΔZ::randType
  @unpack f,g,u,t,dt,T,alg,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,dtmax,dtmin,internalnorm,discard_length,progressbar,progressbar_name,progress_steps,progress_message,unstable_check,rands,sqdt,W,Z,tableau = integrator

  progressbar && (prog = ProgressBar(name=progressbar_name))
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
  if iter > maxiters
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
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,copy(u))
    push!(ts,t)
    push!(Ws,copy(W))
  end
end

@def sde_loopfooter begin
  if adaptive
    standard = abs(1/(γ*EEst))^(2)
    if isinf(standard)
        q = qmax
    else
       q = min(qmax,max(standard,eps()))
    end
    if q > 1
      acceptedIters += 1
      t = t + dt
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
          c = min(dtmax,q*dt)
          dt = max(min(c,abs(T-t)),dtmin)#abs to fix complex sqrt issue at end
          #dt = min(c,abs(T-t))
          sqdt = sqrt(dt)
          ΔW = sqdt*next(rands)
          ΔZ = sqdt*next(rands)
        end
      elseif adaptivealg==:RSwM2 || adaptivealg==:RSwM3
        c = min(dtmax,q*dt)
        dt = max(min(c,abs(T-t)),dtmin) #abs to fix complex sqrt issue at end
        sqdt = sqrt(dt)
        if !(uType <: AbstractArray)
          dttmp = 0.0; ΔW = 0.0; ΔZ = 0.0
        else
          dttmp = 0.0; ΔW = zeros(u); ΔZ = zeros(u)
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
      if adaptivealg==:RSwM1 || adaptivealg==:RSwM2
        ΔWtmp = q*ΔW + sqrt((1-q)*q*dt)*next(rands)
        ΔZtmp = q*ΔZ + sqrt((1-q)*q*dt)*next(rands)
        cutLength = dt-q*dt
        if cutLength > discard_length
          push!(S₁,(cutLength,ΔW-ΔWtmp,ΔZ-ΔZtmp))
        end
        if length(S₁) > max_stack_size
            max_stack_size = length(S₁)
        end
        ΔW = ΔWtmp
        ΔZ = ΔZtmp
        dt = q*dt
      else # RSwM3
        if !(uType <: AbstractArray)
          dttmp = 0.0; ΔWtmp = 0.0; ΔZtmp = 0.0
        else
          dttmp = 0.0; ΔWtmp = zeros(u); ΔZtmp = zeros(u)
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
        dt = q*dt
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
    if !(typeof(alg) <: EM) || !(typeof(alg) <: RKMil)
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
  if progressbar && iter%progress_steps==0
    msg(prog,progress_message(dt,t,u))
    progress(prog,t/T)
  end
end

@def sde_adaptiveprelim begin
  max_stack_size = 0
  max_stack_size2 = 0
  if adaptive
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
  progressbar && done(prog)
  u,t,W,timeseries,ts,Ws,max_stack_size,max_stack_size2
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{EM,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @inbounds while t<T
    @sde_loopheader
    u = u + dt.*f(t,u) + g(t,u).*ΔW
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{EM,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  utmp1 = zeros(u); utmp2 = zeros(u)
  @inbounds while t<T
    @sde_loopheader
    f(t,u,utmp1)
    g(t,u,utmp2)
    for i in eachindex(u)
      u[i] = u[i] + dt*utmp1[i] + utmp2[i]*ΔW[i]
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRI,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages = length(α)
  H0 = Vector{typeof(u)}(0)
  H1 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,zeros(u))
    push!(H1,zeros(u))
  end
  #TODO Reduce memory
  A0temp::uType = zeros(u); A1temp::uType = zeros(u)
  B0temp::uType = zeros(u); B1temp::uType = zeros(u)
  A0temp2::uType = zeros(u); A1temp2::uType = zeros(u)
  B0temp2::uType = zeros(u); B1temp2::uType = zeros(u)
  atemp::uType = zeros(u); btemp::uType = zeros(u)
  E₁::uType = zeros(u); E₂::uType = zeros(u); E₁temp::uType = zeros(u)
  ftemp::uType = zeros(u); gtemp::uType = zeros(u)
  chi1::randType = similar(ΔW); chi2::randType = similar(ΔW); chi3::randType = similar(ΔW)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = .5*(ΔW[i].^2 - dt)/sqdt #I_(1,1)/sqrt(h)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
      chi3[i] = 1/6 * (ΔW[i].^3 - 3*ΔW[i]*dt)/dt #I_(1,1,1)/h
    end
    for i=1:stages
      H0[i][:]=zero(uEltype)
      H1[i][:]=zero(uEltype)
    end
    for i = 1:stages
      A0temp[:]=zero(uEltype)
      B0temp[:]=zero(uEltype)
      A1temp[:]=zero(uEltype)
      B1temp[:]=zero(uEltype)
      for j = 1:i-1
        f(t + c₀[j]*dt,H0[j],ftemp)
        g(t + c₁[j]*dt,H1[j],gtemp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftemp[k]
          B0temp[k] += B₀[i,j]*gtemp[k]
          A1temp[k] += A₁[i,j]*ftemp[k]
          B1temp[k] += B₁[i,j]*gtemp[k]
        end
      end
      H0[i] = u + A0temp*dt + B0temp.*chi2
      H1[i] = u + A1temp*dt + B1temp*sqdt
    end
    atemp[:]=zero(uEltype)
    btemp[:]=zero(uEltype)
    E₂[:]=zero(uEltype)
    E₁temp[:]=zero(uEltype)
    for i = 1:stages
      f(t+c₀[i]*dt,H0[i],ftemp)
      g(t+c₁[i]*dt,H1[i],gtemp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftemp[j]
        btemp[j] += (β₁[i]*ΔW[j] + β₂[i]*chi1[j])*gtemp[j]
        E₂[j]    += (β₃[i]*chi2[j] + β₄[i]*chi3[j])*gtemp[j]
      end
      if i<3 #1 or 2
        for j in eachindex(u)
          E₁temp[j] += ftemp[j]
        end
      end
    end

    for i in eachindex(u)
      E₁[i] = dt*E₁temp[i]
    end

    if adaptive
      for i in eachindex(u)
        EEsttmp[i] = (δ*E₁[i]+E₂[i])/(abstol + u[i]*reltol)
      end
      EEst = internalnorm(EEsttmp)
      for i in eachindex(u)
        utmp[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    else
      for i in eachindex(u)
        u[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRIW1,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  chi1::randType = similar(ΔW)
  chi2::randType = similar(ΔW)
  chi3::randType = similar(ΔW)
  fH01o4::uType = zeros(u)
  g₁o2::uType = zeros(u)
  H0::uType = zeros(u)
  H11::uType = zeros(u)
  H12::uType = zeros(u)
  H13::uType = zeros(u)
  g₂o3::uType = zeros(u)
  Fg₂o3::uType = zeros(u)
  g₃o3::uType = zeros(u)
  Tg₃o3::uType = zeros(u)
  mg₁::uType = zeros(u)
  E₁::uType = zeros(u)
  E₂::uType = zeros(u)
  fH01::uType = zeros(u); fH02::uType = zeros(u)
  g₁::uType = zeros(u); g₂::uType = zeros(u); g₃::uType = zeros(u); g₄::uType = zeros(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    for i in eachindex(u)
      chi1[i] = (ΔW[i].^2 - dt)/2sqdt #I_(1,1)/sqrt(h)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      chi3[i] = (ΔW[i].^3 - 3ΔW[i]*dt)/6dt #I_(1,1,1)/h
    end
    f(t,u,fH01);fH01*=dt
    g(t,u,g₁)
    dto4 = dt/4
    for i in eachindex(u)
      fH01o4[i] = fH01[i]/4
      g₁o2[i] = g₁[i]/2
      H0[i] =  u[i] + 3*(fH01o4[i]  + chi2[i]*g₁o2[i])
      H11[i] = u[i] + fH01o4[i]   + sqdt*g₁o2[i]
      H12[i] = u[i] + fH01[i]     - sqdt*g₁[i]
    end
    g(t+dto4,H11,g₂)
    g(t+dt,H12,g₃)
    for i in eachindex(u)
      H13[i] = u[i] + fH01o4[i] + sqdt*(-5g₁[i] + 3g₂[i] + g₃[i]/2)
    end

    g(t+dto4,H13,g₄)
    f(t+3dto4,H0,fH02); fH02*=dt
    for i in eachindex(u)
      g₂o3[i] = g₂[i]/3
      Fg₂o3[i] = 4g₂o3[i]
      g₃o3[i] = g₃[i]/3
      Tg₃o3[i] = 2g₃o3[i]
      mg₁[i] = -g₁[i]
      E₁[i] = fH01[i]+fH02[i]
      E₂[i] = chi2[i]*(2g₁[i] - Fg₂o3[i] - Tg₃o3[i]) + chi3[i]*(2mg₁[i] + 5g₂o3[i] - Tg₃o3[i] + g₄[i])
    end

    if adaptive
      for i in eachindex(u)
        EEsttmp[i] = (δ*E₁[i]+E₂[i])/(abstol + u[i]*reltol)
      end
      EEst = internalnorm(EEsttmp)
      for i in eachindex(u)
        utmp[i] = u[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
      end
    else
      for i in eachindex(u)
        u[i] = u[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRIW1,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  local H0::uType
  @sde_adaptiveprelim
  local fH01::uType; local g₁::uType
  local fH01o4::uType; local g₁o2::uType
  local H11::uType; local H12::uType
  local g₂::uType; local g₃::uType; local g₄::uType
  local H13::uType; local fH02::uType
  local g₂o3::uType; local Fg₂o3::uType
  local g₃o3::uType; local Tg₃o3::uType
  local mg₁::uType; local E₁::uType; local E₂::uType
  @inbounds while t<T
    @sde_loopheader

    chi1 = (ΔW.^2 - dt)/2sqdt #I_(1,1)/sqrt(h)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = (ΔW.^3 - 3ΔW*dt)/6dt #I_(1,1,1)/h
    fH01 = dt*f(t,u)

    g₁ = g(t,u)
    fH01o4 = fH01/4
    dto4 = dt/4
    g₁o2 = g₁/2
    H0 =  u + 3*(fH01o4  + chi2.*g₁o2)
    H11 = u + fH01o4   + sqdt*g₁o2
    H12 = u + fH01     - sqdt*g₁
    g₂ = g(t+dto4,H11)
    g₃ = g(t+dt,H12)
    H13 = u + fH01o4 + sqdt*(-5g₁ + 3g₂ + g₃/2)


    g₄ = g(t+dto4,H13)
    fH02 = dt*f(t+3dto4,H0)

    g₂o3 = g₂/3
    Fg₂o3 = 4g₂o3
    g₃o3 = g₃/3
    Tg₃o3 = 2g₃o3
    mg₁ = -g₁
    E₁ = fH01+fH02
    E₂ = chi2.*(2g₁ - Fg₂o3 - Tg₃o3) + chi3.*(2mg₁ + 5g₂o3 - Tg₃o3 + g₄)

    if adaptive
      EEst = abs((δ*E₁+E₂)/(abstol + u*reltol))
      utmp = u + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
    else
      u = u + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRI,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = tableau
  stages::Int = length(α)
  H0 = Array{typeof(u)}(stages)
  H1 = Array{typeof(u)}(stages)
  local A0temp::uType; local A1temp::uType
  local B0temp::uType; local B1temp::uType
  local atemp::uType;  local btemp::uType
  local E₁::uType; local E₂::uType
  local E₁temp::uType; local ftemp::uType
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - dt)/sqdt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*dt)/dt #I_(1,1,1)/h

    H0[:]=zero(typeof(u))
    H1[:]=zero(typeof(u))
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      A1temp = zero(u)
      B1temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*f(t + c₀[j]*dt,H0[j])
        B0temp += B₀[i,j]*g(t + c₁[j]*dt,H1[j])
        A1temp += A₁[i,j]*f(t + c₀[j]*dt,H0[j])
        B1temp += B₁[i,j]*g(t + c₁[j]*dt,H1[j])
      end
      H0[i] = u + A0temp*dt + B0temp.*chi2
      H1[i] = u + A1temp*dt + B1temp*sqdt
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)
    for i = 1:stages
      ftemp = f(t+c₀[i]*dt,H0[i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW + β₂[i]*chi1).*g(t+c₁[i]*dt,H1[i])
      E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*g(t+c₁[i]*dt,H1[i])
      if i<3 #1 or 2
        E₁temp += ftemp
      end
    end
    E₁ = dt*E₁temp


    if adaptive
      EEst = abs((δ*E₁+E₂)./(abstol + u*reltol))
      utmp = u + dt*atemp + btemp + E₂
    else
      u = u + dt*atemp + btemp + E₂
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{RKMil,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  du1::uType = zeros(u); du2::uType = zeros(u)
  K::uType = zeros(u); utilde::uType = zeros(u); L::uType = zeros(u)
  @inbounds while t<T
    @sde_loopheader
    f(t,u,du1)
    g(t,u,L)
    for i in eachindex(u)
      K[i] = u[i] + dt*du1[i]
      utilde[i] = K[i] + L[i]*sqdt
    end
    g(t,utilde,du2)
    for i in eachindex(u)
      u[i] = K[i]+L[i]*ΔW[i]+(du2[i]-L[i])./(2sqdt).*(ΔW[i].^2 - dt)
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{RKMil,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  local L::uType; local K::uType; local utilde::uType
  @inbounds while t<T
    @sde_loopheader

    K = u + dt.*f(t,u)
    L = g(t,u)
    utilde = K + L*sqdt
    u = K+L*ΔW+(g(t,utilde)-g(t,u))/(2sqdt)*(ΔW^2 - dt)

    @sde_loopfooter
  end
  @sde_postamble
end


function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRA1,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  H0 = Array{uEltype}(size(u)...,2)
  local k₁::uType; local k₂::uType; local E₁::uType; local E₂::uType
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    gpdt = g(t+dt,u)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    k₁ = dt*f(t,u)
    k₂ = dt*f(t+3dt/4,u+3k₁/4 + 3chi2*g(t+dt,u)/2)
    E₁ = k₁ + k₂
    E₂ = chi2.*(g(t,u)-gpdt) #Only for additive!

    if adaptive
      EEst = abs((δ*E₁+E₂)./(abstol + u*reltol))
      utmp = u + k₁/3 + 2k₂/3 + E₂ + ΔW*gpdt
    else
      u = u + k₁/3 + 2k₂/3 + E₂ + ΔW*gpdt
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRA1,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble

  H0 = Array{uEltype}(size(u)...,2)
  chi2::randType = similar(ΔW)
  tmp1::uType = zeros(u)
  E₁::uType = zeros(u); gt::uType = zeros(u); gpdt::uType = zeros(u)
  E₂::uType = zeros(u); k₁::uType = zeros(u); k₂::uType = zeros(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    g(t,u,gt)
    g(t+dt,u,gpdt)
    f(t,u,k₁); k₁*=dt
    for i in eachindex(u)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      tmp1[i] = u[i]+3k₁[i]/4 + 3chi2[i]*gpdt[i]/2
    end

    f(t+3dt/4,tmp1,k₂); k₂*=dt

    for i in eachindex(u)
      E₁[i] = k₁[i] + k₂[i]
      E₂[i] = chi2[i]*(gt[i]-gpdt[i]) #Only for additive!
    end

    if adaptive
      for i in eachindex(u)
        EEsttmp[i] = (δ*E₁[i]+E₂[i])/(abstol + u[i]*reltol)
      end
      EEst = internalnorm(EEsttmp)
      for i in eachindex(u)
        utmp[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*gpdt[i]
      end
    else
      for i in eachindex(u)
        u[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*gpdt[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRA,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = tableau
  stages::Int = length(α)
  H0 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,zeros(u))
  end
  A0temp::uType = zeros(u); B0temp::uType = zeros(u)
  ftmp::uType = zeros(u); gtmp::uType = zeros(u); chi2::uType = zeros(u)
  atemp::uType = zeros(u); btemp::uType = zeros(u); E₂::uType = zeros(u); E₁temp::uType = zeros(u)
  E₁::uType = zeros(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    for i in eachindex(u)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
    end
    for i in 1:stages
      H0[i][:]=zero(uEltype)
    end
    for i = 1:stages
      A0temp[:] = zero(uEltype)
      B0temp[:] = zero(uEltype)
      for j = 1:i-1
        f(t + c₀[j]*dt,H0[j],ftmp)
        g(t + c₁[j]*dt,H0[j],gtmp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftmp[k]
          B0temp[k] += B₀[i,j]*gtmp[k]
        end
      end
      for j in eachindex(u)
        H0[i][j] = u[j] + A0temp[j]*dt + B0temp[j]*chi2[j]
      end
    end
    atemp[:] = zero(uEltype)
    btemp[:] = zero(uEltype)
    E₂[:]    = zero(uEltype)
    E₁temp[:]= zero(uEltype)

    for i = 1:stages
      f(t+c₀[i]*dt,H0[i],ftmp)
      g(t+c₁[i]*dt,H0[i],gtmp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftmp[j]
        btemp[j] += (β₁[i]*ΔW[j])*gtmp[j]
        E₂[j]    += (β₂[i]*chi2[j])*gtmp[j]
        E₁temp[j] += ftmp[j]
      end
    end
    for i in eachindex(u)
      E₁[i] = dt*E₁temp[i]
    end

    if adaptive
      for i in eachindex(u)
        EEsttmp[i] = (δ*E₁[i]+E₂[i])/(abstol + u[i]*reltol)
      end
      EEst = internalnorm(EEsttmp)
      for i in eachindex(u)
        utmp[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    else
      for i in eachindex(u)
        u[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRA,uType,uEltype,Nm1,N,tType,tableauType,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @sde_sratableaupreamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = tableau
  stages::Int = length(α)
  H0 = Array{uEltype}(stages)
  local atemp::uType; local btemp::uType
  local E₂::uType; local E₁::uType; local E₁temp::uType
  local ftemp::uType; local A0temp::uType; local B0temp::uType
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*f(t + c₀[j]*dt,H0[j])
        B0temp += B₀[i,j]*g(t + c₁[j]*dt,H0[j]) #H0[..,i] argument ignored
      end
      H0[i] = u + A0temp*dt + B0temp.*chi2
    end

    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)

    for i = 1:stages
      ftemp = f(t+c₀[i]*dt,H0[i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW ).*g(t+c₁[i]*dt,H0[i]) #H0[i] argument ignored
      E₂    += (β₂[i]*chi2).*g(t+c₁[i]*dt,H0[i]) #H0[i] argument ignored
    end

    if adaptive
      E₁temp += ftemp
      E₁ = dt*E₁temp
      EEst = abs((δ*E₁+E₂)./(abstol + u*reltol))
      utmp = u + dt*atemp + btemp + E₂
    else
      u = u + dt*atemp + btemp + E₂
    end
    @sde_loopfooter
  end
  @sde_postamble
end
