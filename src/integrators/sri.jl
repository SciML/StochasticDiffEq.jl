@inline function perform_step!(integrator,cache::SRICache,f=integrator.f)
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,error_terms = cache.tab
  @unpack H0,H1,A0temp,A1temp,B0temp,B1temp,A0temp2,A1temp2,B0temp2,B1temp2,atemp,btemp,E₁,E₂,E₁temp,ftemp,gtemp,chi1,chi2,chi3,tmp = cache
  @unpack t,dt,uprev,u,ΔW,ΔZ = integrator
  for i in eachindex(u)
    chi1[i] = .5*(ΔW[i].^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
    chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
    chi3[i] = 1/6 * (ΔW[i].^3 - 3*ΔW[i]*dt)/dt #I_(1,1,1)/h
  end
  for i=1:stages
    fill!(H0[i],zero(eltype(integrator.u)))
    fill!(H1[i],zero(eltype(integrator.u)))
  end
  for i = 1:stages
    fill!(A0temp,zero(eltype(integrator.u)))
    fill!(B0temp,zero(eltype(integrator.u)))
    fill!(A1temp,zero(eltype(integrator.u)))
    fill!(B1temp,zero(eltype(integrator.u)))
    for j = 1:i-1
      integrator.f(@muladd(t + c₀[j]*dt),H0[j],ftemp)
      integrator.g(@muladd(t + c₁[j]*dt),H1[j],gtemp)
      for k in eachindex(u)
        A0temp[k] = @muladd A0temp[k] + A₀[j,i]*ftemp[k]
        B0temp[k] = @muladd B0temp[k] + B₀[j,i]*gtemp[k]
        A1temp[k] = @muladd A1temp[k] + A₁[j,i]*ftemp[k]
        B1temp[k] = @muladd B1temp[k] + B₁[j,i]*gtemp[k]
      end
    end
    H0[i] = uprev + A0temp*dt + B0temp.*chi2
    H1[i] = uprev + A1temp*dt + B1temp*integrator.sqdt
  end
  fill!(atemp,zero(eltype(integrator.u)))
  fill!(btemp,zero(eltype(integrator.u)))
  fill!(E₂,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))
  for i = 1:stages
    integrator.f(@muladd(t+c₀[i]*dt),H0[i],ftemp)
    integrator.g(@muladd(t+c₁[i]*dt),H1[i],gtemp)
    for j in eachindex(u)
      atemp[j] = @muladd atemp[j] + α[i]*ftemp[j]
      btemp[j] = @muladd btemp[j] + (β₁[i]*ΔW[j] + β₂[i]*chi1[j])*gtemp[j]
      E₂[j]    = @muladd E₂[j]    + (β₃[i]*chi2[j] + β₄[i]*chi3[j])*gtemp[j]
    end
    if i <= error_terms
      for j in eachindex(u)
        E₁temp[j] += ftemp[j]
      end
    end
  end

  for i in eachindex(u)
    E₁[i] = dt*E₁temp[i]
  end

  for i in eachindex(u)
    u[i] = uprev[i] + @muladd(dt*atemp[i] + btemp[i]) + E₂[i]
  end

  if integrator.opts.adaptive
    for i in eachindex(u)
      tmp[i] = @muladd(integrator.opts.delta*E₁[i]+E₂[i])/@muladd(integrator.opts.abstol + max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SRIW1Cache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW,ΔZ = integrator
  @unpack chi1,chi2,chi3,fH01o4,g₁o2,H0,H11,H12,H13,g₂o3,Fg₂o3,g₃o3,Tg₃o3,mg₁,E₁,E₂,fH01,fH02,g₁,g₂,g₃,g₄,tmp = cache
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
    H0[i] =  @muladd uprev[i] + 3*(fH01o4[i]  + chi2[i]*g₁o2[i])
    H11[i] = @muladd uprev[i] + fH01o4[i]   + integrator.sqdt*g₁o2[i]
    H12[i] = @muladd uprev[i] + fH01[i]     - integrator.sqdt*g₁[i]
  end
  integrator.g(t+dto4,H11,g₂)
  integrator.g(t+dt,H12,g₃)
  for i in eachindex(u)
    H13[i] = @muladd uprev[i] + fH01o4[i] + integrator.sqdt*(-5g₁[i] + 3g₂[i] + g₃[i]/2)
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
    E₂[i] = @muladd chi2[i]*(2g₁[i] - Fg₂o3[i] - Tg₃o3[i]) + chi3[i]*(2mg₁[i] + 5g₂o3[i] - Tg₃o3[i] + g₄[i])
  end

  for i in eachindex(u)
    u[i] = @muladd uprev[i] +  (fH01[i] + 2fH02[i])/3 + ΔW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
  end

  if integrator.opts.adaptive
    for i in eachindex(u)
      tmp[i] = @muladd(integrator.opts.delta*E₁[i]+E₂[i])/@muladd(integrator.opts.abstol + max(abs(uprev[i]),abs(u[i]))*integrator.opts.reltol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SRIW1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW,ΔZ = integrator
  chi1 = (ΔW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
  chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
  chi3 = (ΔW.^3 - 3ΔW*dt)/6dt #I_(1,1,1)/h
  fH01 = dt*integrator.f(t,uprev)

  g₁ = integrator.g(t,uprev)
  fH01o4 = fH01/4
  dto4 = dt/4
  g₁o2 = g₁/2
  H0 =  @muladd uprev + 3*(fH01o4  + chi2.*g₁o2)
  H11 = @muladd uprev + fH01o4   + integrator.sqdt*g₁o2
  H12 = @muladd uprev + fH01     - integrator.sqdt*g₁
  g₂ = integrator.g(t+dto4,H11)
  g₃ = integrator.g(t+dt,H12)
  H13 = @muladd uprev + fH01o4 + integrator.sqdt*(-5g₁ + 3g₂ + g₃/2)


  g₄ = integrator.g(t+dto4,H13)
  fH02 = dt*integrator.f(t+3dto4,H0)

  g₂o3 = g₂/3
  Fg₂o3 = 4g₂o3
  g₃o3 = g₃/3
  Tg₃o3 = 2g₃o3
  mg₁ = -g₁
  E₁ = fH01+fH02
  E₂ = @muladd chi2.*(2g₁ - Fg₂o3 - Tg₃o3) + chi3.*(2mg₁ + 5g₂o3 - Tg₃o3 + g₄)

  u = uprev + (fH01 + 2fH02)/3 + ΔW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)/(@muladd(integrator.opts.abstol + max(abs(uprev),abs(u))*integrator.opts.reltol)))
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SRIConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW,ΔZ = integrator
  error_terms = integrator.alg.error_terms
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,H0,H1,error_terms = cache
  chi1 = .5*(ΔW.^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
  chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
  chi3 = 1/6 * (ΔW.^3 - 3*ΔW*dt)/dt #I_(1,1,1)/h

  fill!(H0,zero(typeof(u)))
  fill!(H1,zero(typeof(u)))
  for i = 1:stages
    A0temp = zero(u)
    B0temp = zero(u)
    A1temp = zero(u)
    B1temp = zero(u)
    for j = 1:i-1
      A0temp = @muladd A0temp + A₀[j,i]*integrator.f(t + c₀[j]*dt,H0[j])
      B0temp = @muladd B0temp + B₀[j,i]*integrator.g(t + c₁[j]*dt,H1[j])
      A1temp = @muladd A1temp + A₁[j,i]*integrator.f(t + c₀[j]*dt,H0[j])
      B1temp = @muladd B1temp + B₁[j,i]*integrator.g(t + c₁[j]*dt,H1[j])
    end
    H0[i] = @muladd uprev + A0temp*dt + B0temp.*chi2
    H1[i] = @muladd uprev + A1temp*dt + B1temp*integrator.sqdt
  end
  atemp = zero(u)
  btemp = zero(u)
  E₂    = zero(u)
  E₁temp= zero(u)
  for i = 1:stages
    ftemp = integrator.f(t+c₀[i]*dt,H0[i])
    atemp = @muladd atemp + α[i]*ftemp
    btemp = @muladd btemp + (β₁[i]*ΔW + β₂[i]*chi1).*integrator.g(t+c₁[i]*dt,H1[i])
    E₂    = @muladd E₂    + (β₃[i]*chi2 + β₄[i]*chi3).*integrator.g(t+c₁[i]*dt,H1[i])
    if i <= error_terms #1 or 2
      E₁temp += ftemp
    end
  end
  E₁ = dt*E₁temp

  u = @muladd(uprev + dt*atemp) + btemp + E₂

  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max(abs(uprev),abs(u))*integrator.opts.reltol))
  end
  @pack integrator = t,dt,u
end
