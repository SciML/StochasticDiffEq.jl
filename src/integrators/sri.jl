function sde_solve{algType<:SRI,uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = integrator.alg.tableau
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
      chi1[i] = .5*(ΔW[i].^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
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
        integrator.f(t + c₀[j]*dt,H0[j],ftemp)
        integrator.g(t + c₁[j]*dt,H1[j],gtemp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftemp[k]
          B0temp[k] += B₀[i,j]*gtemp[k]
          A1temp[k] += A₁[i,j]*ftemp[k]
          B1temp[k] += B₁[i,j]*gtemp[k]
        end
      end
      H0[i] = uprev + A0temp*dt + B0temp.*chi2
      H1[i] = uprev + A1temp*dt + B1temp*integrator.sqdt
    end
    atemp[:]=zero(uEltype)
    btemp[:]=zero(uEltype)
    E₂[:]=zero(uEltype)
    E₁temp[:]=zero(uEltype)
    for i = 1:stages
      integrator.f(t+c₀[i]*dt,H0[i],ftemp)
      integrator.g(t+c₁[i]*dt,H1[i],gtemp)
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

    if integrator.opts.adaptive
      for i in eachindex(u)
        u[i] = uprev[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
      for i in eachindex(u)
        EEsttmp[i] = (integrator.opts.delta*E₁[i]+E₂[i])/(integrator.opts.abstol + max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
      end
      EEst = integrator.opts.internalnorm(EEsttmp)
    else
      for i in eachindex(u)
        u[i] = uprev[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRIW1,uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  chi1::randType = similar(ΔW)
  chi2::randType = similar(ΔW)
  chi3::randType = similar(ΔW)
  fH01o4::uType = zeros(uprev)
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
      chi1[i] = (ΔW[i].^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      chi3[i] = (ΔW[i].^3 - 3ΔW[i]*dt)/6dt #I_(1,1,1)/h
    end
    integrator.f(t,uprev,fH01)
    for i in eachindex(u)
      fH01[i] = dt*fH01[i]
    end
    integrator.g(t,uprev,g₁)
    dto4 = dt/4
    for i in eachindex(u)
      fH01o4[i] = fH01[i]/4
      g₁o2[i] = g₁[i]/2
      H0[i] =  uprev[i] + 3*(fH01o4[i]  + chi2[i]*g₁o2[i])
      H11[i] = uprev[i] + fH01o4[i]   + integrator.sqdt*g₁o2[i]
      H12[i] = uprev[i] + fH01[i]     - integrator.sqdt*g₁[i]
    end
    integrator.g(t+dto4,H11,g₂)
    integrator.g(t+dt,H12,g₃)
    for i in eachindex(u)
      H13[i] = uprev[i] + fH01o4[i] + integrator.sqdt*(-5g₁[i] + 3g₂[i] + g₃[i]/2)
    end

    integrator.g(t+dto4,H13,g₄)
    integrator.f(t+3dto4,H0,fH02)
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

    if integrator.opts.adaptive
      for i in eachindex(u)
        u[i] = uprev[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
      end
      for i in eachindex(u)
        EEsttmp[i] = (integrator.opts.delta*E₁[i]+E₂[i])/(integrator.opts.abstol + max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
      end
      EEst = integrator.opts.internalnorm(EEsttmp)
    else
      for i in eachindex(u)
        u[i] = uprev[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRIW1,uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi1 = (ΔW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = (ΔW.^3 - 3ΔW*dt)/6dt #I_(1,1,1)/h
    fH01 = dt*integrator.f(t,uprev)

    g₁ = integrator.g(t,uprev)
    fH01o4 = fH01/4
    dto4 = dt/4
    g₁o2 = g₁/2
    H0 =  uprev + 3*(fH01o4  + chi2.*g₁o2)
    H11 = uprev + fH01o4   + integrator.sqdt*g₁o2
    H12 = uprev + fH01     - integrator.sqdt*g₁
    g₂ = integrator.g(t+dto4,H11)
    g₃ = integrator.g(t+dt,H12)
    H13 = uprev + fH01o4 + integrator.sqdt*(-5g₁ + 3g₂ + g₃/2)


    g₄ = integrator.g(t+dto4,H13)
    fH02 = dt*integrator.f(t+3dto4,H0)

    g₂o3 = g₂/3
    Fg₂o3 = 4g₂o3
    g₃o3 = g₃/3
    Tg₃o3 = 2g₃o3
    mg₁ = -g₁
    E₁ = fH01+fH02
    E₂ = chi2.*(2g₁ - Fg₂o3 - Tg₃o3) + chi3.*(2mg₁ + 5g₂o3 - Tg₃o3 + g₄)

    if integrator.opts.adaptive
      u = uprev + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
      EEst = abs((integrator.opts.delta*E₁+E₂)/(integrator.opts.abstol + max(abs(uprev),abs(u))*integrator.opts.reltol))
    else
      u = uprev + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRI,uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = integrator.alg.tableau
  stages::Int = length(α)
  H0 = Array{typeof(u)}(stages)
  H1 = Array{typeof(u)}(stages)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
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
        A0temp += A₀[i,j]*integrator.f(t + c₀[j]*dt,H0[j])
        B0temp += B₀[i,j]*integrator.g(t + c₁[j]*dt,H1[j])
        A1temp += A₁[i,j]*integrator.f(t + c₀[j]*dt,H0[j])
        B1temp += B₁[i,j]*integrator.g(t + c₁[j]*dt,H1[j])
      end
      H0[i] = uprev + A0temp*dt + B0temp.*chi2
      H1[i] = uprev + A1temp*dt + B1temp*integrator.sqdt
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)
    for i = 1:stages
      ftemp = integrator.f(t+c₀[i]*dt,H0[i])
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW + β₂[i]*chi1).*integrator.g(t+c₁[i]*dt,H1[i])
      E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*integrator.g(t+c₁[i]*dt,H1[i])
      if i<3 #1 or 2
        E₁temp += ftemp
      end
    end
    E₁ = dt*E₁temp


    if integrator.opts.adaptive
      u = uprev + dt*atemp + btemp + E₂
      EEst = abs((integrator.opts.delta*E₁+E₂)./(integrator.opts.abstol + max(abs(uprev),abs(u))*integrator.opts.reltol))
    else
      u = uprev + dt*atemp + btemp + E₂
    end

    @sde_loopfooter
  end
  @sde_postamble
end
