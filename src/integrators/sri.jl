#=
@inline function perform_step!(integrator,cache::SRICache,f=integrator.f)
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,error_terms = cache.tab
  @unpack H0,H1,A0temp,A1temp,B0temp,B1temp,A0temp2,A1temp2,B0temp2,B1temp2,atemp,btemp,E₁,E₂,E₁temp,ftemp,gtemp,chi1,chi2,chi3,tmp = cache
  @unpack t,dt,uprev,u,W = integrator
  @. chi1 = .5*(W.dW.^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
  @. chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  @. chi3 = 1/6 * (W.dW.^3 - 3*W.dW*dt)/dt #I_(1,1,1)/h
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
      @. A0temp = @muladd A0temp + A₀[j,i]*ftemp
      @. B0temp = @muladd B0temp + B₀[j,i]*gtemp
      @. A1temp = @muladd A1temp + A₁[j,i]*ftemp
      @. B1temp = @muladd B1temp + B₁[j,i]*gtemp
    end
    @. H0[i] = uprev + A0temp*dt + B0temp*chi2
    @. H1[i] = uprev + A1temp*dt + B1temp*integrator.sqdt
  end
  fill!(atemp,zero(eltype(integrator.u)))
  fill!(btemp,zero(eltype(integrator.u)))
  fill!(E₂,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))
  for i = 1:stages
    integrator.f(@muladd(t+c₀[i]*dt),H0[i],ftemp)
    integrator.g(@muladd(t+c₁[i]*dt),H1[i],gtemp)
    @. atemp = @muladd atemp + α[i]*ftemp
    @. btemp = @muladd btemp + (β₁[i]*W.dW + β₂[i]*chi1)*gtemp
    @. E₂    = @muladd E₂    + (β₃[i]*chi2 + β₄[i]*chi3)*gtemp
    if i <= error_terms
      @. E₁temp += ftemp
    end
  end

  @. E₁ = dt*E₁temp
  @. u = uprev + @muladd(dt*atemp + btemp) + E₂

  if integrator.opts.adaptive
    @. tmp = @muladd(integrator.opts.delta*E₁+E₂)/@muladd(integrator.opts.abstol + max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRICache,f=integrator.f)
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,error_terms = cache.tab
  @unpack H0,H1,A0temp,A1temp,B0temp,B1temp,A0temp2,A1temp2,B0temp2,B1temp2,atemp,btemp,E₁,E₂,E₁temp,ftemp,gtemp,chi1,chi2,chi3,tmp = cache
  @unpack t,dt,uprev,u,W = integrator
  @tight_loop_macros for i in eachindex(u)
    @inbounds chi1[i] = .5*(W.dW[i].^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
    @inbounds chi2[i] = .5*(W.dW[i] + W.dZ[i]/sqrt(3)) #I_(1,0)/h
    @inbounds chi3[i] = 1/6 * (W.dW[i].^3 - 3*W.dW[i]*dt)/dt #I_(1,1,1)/h
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
      @tight_loop_macros for k in eachindex(u)
        @inbounds A0temp[k] = @muladd A0temp[k] + A₀[j,i]*ftemp[k]
        @inbounds B0temp[k] = @muladd B0temp[k] + B₀[j,i]*gtemp[k]
        @inbounds A1temp[k] = @muladd A1temp[k] + A₁[j,i]*ftemp[k]
        @inbounds B1temp[k] = @muladd B1temp[k] + B₁[j,i]*gtemp[k]
      end
    end
    @tight_loop_macros for k in eachindex(u)
      @inbounds H0[i][k] = uprev[k] + A0temp[k]*dt + B0temp[k]*chi2[k]
      @inbounds H1[i][k] = uprev[k] + A1temp[k]*dt + B1temp[k]*integrator.sqdt
    end
  end
  fill!(atemp,zero(eltype(integrator.u)))
  fill!(btemp,zero(eltype(integrator.u)))
  fill!(E₂,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))
  for i = 1:stages
    integrator.f(@muladd(t+c₀[i]*dt),H0[i],ftemp)
    integrator.g(@muladd(t+c₁[i]*dt),H1[i],gtemp)
    @tight_loop_macros for j in eachindex(u)
      @inbounds atemp[j] = @muladd atemp[j] + α[i]*ftemp[j]
      @inbounds btemp[j] = @muladd btemp[j] + (β₁[i]*W.dW[j] + β₂[i]*chi1[j])*gtemp[j]
      @inbounds E₂[j]    = @muladd E₂[j]    + (β₃[i]*chi2[j] + β₄[i]*chi3[j])*gtemp[j]
    end
    if i <= error_terms
      @tight_loop_macros for j in eachindex(u)
        @inbounds E₁temp[j] += ftemp[j]
      end
    end
  end

  @tight_loop_macros for i in eachindex(u)
    @inbounds E₁[i] = dt*E₁temp[i]
  end

  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i] + @muladd(dt*atemp[i] + btemp[i]) + E₂[i]
  end

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

#=
@inline function perform_step!(integrator,cache::SRIW1Cache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack chi1,chi2,chi3,fH01o4,g₁o2,H0,H11,H12,H13,g₂o3,Fg₂o3,g₃o3,Tg₃o3,mg₁,E₁,E₂,fH01,fH02,g₁,g₂,g₃,g₄,tmp = cache
  @. chi1 = (W.dW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
  @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  @. chi3 = (W.dW.^3 - 3W.dW*dt)/6dt #I_(1,1,1)/h
  integrator.f(t,uprev,fH01)
  @. fH01 = dt*fH01
  integrator.g(t,uprev,g₁)
  dto4 = dt/4
  @. fH01o4 = fH01/4
  @. g₁o2 = g₁/2
  @. H0 =  @muladd uprev + 3*(fH01o4  + chi2*g₁o2)
  @. H11 = @muladd uprev + fH01o4   + integrator.sqdt*g₁o2
  @. H12 = @muladd uprev + fH01     - integrator.sqdt*g₁
  integrator.g(t+dto4,H11,g₂)
  integrator.g(t+dt,H12,g₃)
  @. H13 = @muladd uprev + fH01o4 + integrator.sqdt*(-5g₁ + 3g₂ + g₃/2)

  integrator.g(t+dto4,H13,g₄)
  integrator.f(t+3dto4,H0,fH02)

  @. fH02 = fH02*dt
  @. g₂o3 = g₂/3
  @. Fg₂o3 = 4g₂o3
  @. g₃o3 = g₃/3
  @. Tg₃o3 = 2g₃o3
  @. mg₁ = -g₁
  @. E₁ = fH01+fH02
  @. E₂ = @muladd chi2*(2g₁ - Fg₂o3 - Tg₃o3) + chi3*(2mg₁ + 5g₂o3 - Tg₃o3 + g₄)

  @. u = @muladd uprev +  (fH01 + 2fH02)/3 + W.dW*(mg₁ + Fg₂o3 + Tg₃o3) + chi1*(mg₁ + Fg₂o3 - g₃o3) + E₂

  if integrator.opts.adaptive
    @. tmp = @muladd(integrator.opts.delta*E₁+E₂)/@muladd(integrator.opts.abstol + max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRIW1Cache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack chi1,chi2,chi3,fH01o4,g₁o2,H0,H11,H12,H13,g₂o3,Fg₂o3,g₃o3,Tg₃o3,mg₁,E₁,E₂,fH01,fH02,g₁,g₂,g₃,g₄,tmp = cache

  if typeof(W.dW) <: Union{SArray,Number}
    chi1 = @. (W.dW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
    chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = @. (W.dW.^3 - 3W.dW*dt)/6dt #I_(1,1,1)/h
  else
    @. chi1 = (W.dW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
    @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
    @. chi3 = (W.dW.^3 - 3W.dW*dt)/6dt #I_(1,1,1)/h
  end

  integrator.f(t,uprev,fH01)
  @tight_loop_macros for i in eachindex(u)
    @inbounds fH01[i] = dt*fH01[i]
  end
  integrator.g(t,uprev,g₁)
  dto4 = dt/4
  @tight_loop_macros for i in eachindex(u)
    @inbounds fH01o4[i] = fH01[i]/4
    @inbounds g₁o2[i] = g₁[i]/2
    @inbounds H0[i] =  @muladd uprev[i] + 3*(fH01o4[i]  + chi2[i]*g₁o2[i])
    @inbounds H11[i] = @muladd uprev[i] + fH01o4[i]   + integrator.sqdt*g₁o2[i]
    @inbounds H12[i] = @muladd uprev[i] + fH01[i]     - integrator.sqdt*g₁[i]
  end
  integrator.g(t+dto4,H11,g₂)
  integrator.g(t+dt,H12,g₃)
  @tight_loop_macros for i in eachindex(u)
    @inbounds H13[i] = @muladd uprev[i] + fH01o4[i] + integrator.sqdt*(-5g₁[i] + 3g₂[i] + g₃[i]/2)
  end

  integrator.g(t+dto4,H13,g₄)
  integrator.f(t+3dto4,H0,fH02)
  @tight_loop_macros for i in eachindex(u)
    @inbounds fH02[i] = fH02[i]*dt
    @inbounds g₂o3[i] = g₂[i]/3
    @inbounds Fg₂o3[i] = 4g₂o3[i]
    @inbounds g₃o3[i] = g₃[i]/3
    @inbounds Tg₃o3[i] = 2g₃o3[i]
    @inbounds mg₁[i] = -g₁[i]
    @inbounds E₁[i] = fH01[i]+fH02[i]
    @inbounds E₂[i] = @muladd chi2[i]*(2g₁[i] - Fg₂o3[i] - Tg₃o3[i]) + chi3[i]*(2mg₁[i] + 5g₂o3[i] - Tg₃o3[i] + g₄[i])
  end

  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = @muladd uprev[i] +  (fH01[i] + 2fH02[i])/3 + W.dW[i]*(mg₁[i] + Fg₂o3[i] + Tg₃o3[i]) + chi1[i]*(mg₁[i] + Fg₂o3[i] - g₃o3[i]) + E₂[i]
  end

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

#=
@inline function perform_step!(integrator,cache::SRIW1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  chi1 = @. (W.dW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
  chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  chi3 = @. (W.dW.^3 - 3W.dW*dt)/6dt #I_(1,1,1)/h
  fH01 = dt*integrator.f(t,uprev)

  g₁ = integrator.g(t,uprev)
  fH01o4 = fH01/4
  dto4 = dt/4
  g₁o2 = g₁/2
  H0 =  @. @muladd uprev + 3*(fH01o4  + chi2.*g₁o2)
  H11 = @. @muladd uprev + fH01o4   + integrator.sqdt*g₁o2
  H12 = @. @muladd uprev + fH01     - integrator.sqdt*g₁
  g₂ = integrator.g(t+dto4,H11)
  g₃ = integrator.g(t+dt,H12)
  H13 = @. @muladd uprev + fH01o4 + integrator.sqdt*(-5g₁ + 3g₂ + g₃/2)


  g₄ = integrator.g(t+dto4,H13)
  fH02 = dt*integrator.f(t+3dto4,H0)

  g₂o3 = g₂/3
  Fg₂o3 = 4g₂o3
  g₃o3 = g₃/3
  Tg₃o3 = 2g₃o3
  mg₁ = -g₁
  E₁ = fH01+fH02
  E₂ = @. @muladd chi2.*(2g₁ - Fg₂o3 - Tg₃o3) + chi3.*(2mg₁ + 5g₂o3 - Tg₃o3 + g₄)

  u = @. uprev + (fH01 + 2fH02)/3 + W.dW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
  if integrator.opts.adaptive
    tmp = @. @muladd(integrator.opts.delta*E₁+E₂)/(@muladd(integrator.opts.abstol + max.(abs.(uprev),abs.(u))*integrator.opts.reltol))
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRIW1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  chi1 = (W.dW.^2 - dt)/2integrator.sqdt #I_(1,1)/sqrt(h)
  chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  chi3 = (W.dW.^3 - 3W.dW*dt)/6dt #I_(1,1,1)/h
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

  u = uprev + (fH01 + 2fH02)/3 + W.dW.*(mg₁ + Fg₂o3 + Tg₃o3) + chi1.*(mg₁ + Fg₂o3 - g₃o3) + E₂
  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta.*E₁.+E₂)./(@muladd(integrator.opts.abstol .+ max.(abs.(uprev),abs.(u)).*integrator.opts.reltol)))
  end
  @pack integrator = t,dt,u
end

#=
@inline function perform_step!(integrator,cache::SRIConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  error_terms = integrator.alg.error_terms
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,H0,H1,error_terms = cache
  chi1 = @. .5*(W.dW.^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
  chi2 = @. .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  chi3 = @. 1/6 * (W.dW.^3 - 3*W.dW*dt)/dt #I_(1,1,1)/h

  fill!(H0,zero(typeof(u)))
  fill!(H1,zero(typeof(u)))
  for i = 1:stages
    A0temp = zero(u)
    B0temp = zero(u)
    A1temp = zero(u)
    B1temp = zero(u)
    for j = 1:i-1
      A0temp = @muladd A0temp .+ A₀[j,i].*integrator.f(t + c₀[j]*dt,H0[j])
      B0temp = @muladd B0temp .+ B₀[j,i].*integrator.g(t + c₁[j]*dt,H1[j])
      A1temp = @muladd A1temp .+ A₁[j,i].*integrator.f(t + c₀[j]*dt,H0[j])
      B1temp = @muladd B1temp .+ B₁[j,i].*integrator.g(t + c₁[j]*dt,H1[j])
    end
    H0[i] = @muladd uprev .+ A0temp.*dt .+ B0temp.*chi2
    H1[i] = @muladd uprev .+ A1temp.*dt .+ B1temp.*integrator.sqdt
  end
  atemp = zero(u)
  btemp = zero(u)
  E₂    = zero(u)
  E₁temp= zero(u)
  for i = 1:stages
    ftemp = integrator.f(t+c₀[i]*dt,H0[i])
    atemp = @muladd atemp .+ α[i].*ftemp
    btemp = @muladd btemp .+ (β₁[i].*W.dW + β₂[i].*chi1).*integrator.g(t+c₁[i]*dt,H1[i])
    E₂    = @muladd E₂    .+ (β₃[i].*chi2 + β₄[i].*chi3).*integrator.g(t+c₁[i]*dt,H1[i])
    if i <= error_terms #1 or 2
      E₁temp += ftemp
    end
  end
  E₁ = dt*E₁temp

  u = muladd.(dt,atemp,uprev) .+ btemp .+ E₂

  if integrator.opts.adaptive
    tmp = @. @muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(abs.(uprev),abs.(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRIConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  error_terms = integrator.alg.error_terms
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄,stages,H0,H1,error_terms = cache
  chi1 = .5*(W.dW.^2 - dt)/integrator.sqdt #I_(1,1)/sqrt(h)
  chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  chi3 = 1/6 * (W.dW.^3 - 3*W.dW*dt)/dt #I_(1,1,1)/h

  fill!(H0,zero(typeof(u)))
  fill!(H1,zero(typeof(u)))
  for i = 1:stages
    A0temp = zero(u)
    B0temp = zero(u)
    A1temp = zero(u)
    B1temp = zero(u)
    for j = 1:i-1
      @inbounds A0temp = @muladd A0temp + A₀[j,i]*integrator.f(t + c₀[j]*dt,H0[j])
      @inbounds B0temp = @muladd B0temp + B₀[j,i]*integrator.g(t + c₁[j]*dt,H1[j])
      @inbounds A1temp = @muladd A1temp + A₁[j,i]*integrator.f(t + c₀[j]*dt,H0[j])
      @inbounds B1temp = @muladd B1temp + B₁[j,i]*integrator.g(t + c₁[j]*dt,H1[j])
    end
    @inbounds H0[i] = @muladd uprev + A0temp*dt + B0temp.*chi2
    @inbounds H1[i] = @muladd uprev + A1temp*dt + B1temp*integrator.sqdt
  end
  atemp = zero(u)
  btemp = zero(u)
  E₂    = zero(u)
  E₁temp= zero(u)
  for i = 1:stages
    @inbounds ftemp = integrator.f(t+c₀[i]*dt,H0[i])
    @inbounds atemp = @muladd atemp + α[i]*ftemp
    @inbounds btemp = @muladd btemp + (β₁[i]*W.dW + β₂[i]*chi1).*integrator.g(t+c₁[i]*dt,H1[i])
    @inbounds E₂    = @muladd E₂    + (β₃[i]*chi2 + β₄[i]*chi3).*integrator.g(t+c₁[i]*dt,H1[i])
    if i <= error_terms #1 or 2
      E₁temp += ftemp
    end
  end
  E₁ = dt*E₁temp

  u = @muladd(uprev + dt*atemp) + btemp + E₂

  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(abs.(uprev),abs.(u))*integrator.opts.reltol))
  end
  @pack integrator = t,dt,u
end
