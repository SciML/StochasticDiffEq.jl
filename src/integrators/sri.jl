function sde_solve{algType<:SRI,uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = alg.tableau
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
        utmp[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
      for i in eachindex(u)
        EEsttmp[i] = (δ*E₁[i]+E₂[i])/(abstol + max(abs(u[i]),abs(utmp[i]))*reltol)
      end
      EEst = internalnorm(EEsttmp)
    else
      for i in eachindex(u)
        u[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRIW1,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
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
    f(t,u,fH01)
    for i in eachindex(u)
      fH01[i] = dt*fH01[i]
    end
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
    f(t+3dto4,H0,fH02)
    for i in eachindex(u)
      fH02[i] = fH02[i]*dt
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
        utmp[i] = u[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
      end
      for i in eachindex(u)
        EEsttmp[i] = (δ*E₁[i]+E₂[i])/(abstol + max(abs(u[i]),abs(utmp[i]))*reltol)
      end
      EEst = internalnorm(EEsttmp)
    else
      for i in eachindex(u)
        u[i] = u[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{SRIW1,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
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
      utmp = u + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
      EEst = abs((δ*E₁+E₂)/(abstol + max(abs(u),abs(utmp))*reltol))
    else
      u = u + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRI,uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @sde_sritableaupreamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = alg.tableau
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
      utmp = u + dt*atemp + btemp + E₂
      EEst = abs((δ*E₁+E₂)./(abstol + max(abs(u),abs(utmp))*reltol))
    else
      u = u + dt*atemp + btemp + E₂
    end

    @sde_loopfooter
  end
  @sde_postamble
end
