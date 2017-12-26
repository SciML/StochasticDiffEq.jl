#=
@inline function perform_step!(integrator,cache::SRA1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  gpdt = integrator.g(t+dt,uprev)
  chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  k₁ = @. dt*integrator.f(t,uprev)
  k₂ = @. dt*integrator.f(t+3dt/4,uprev+3k₁/4 + 3chi2*integrator.g(t+dt,uprev)/2)
  E₁ = @. k₁ + k₂
  E₂ = @. chi2.*(integrator.g(t,uprev)-gpdt) #Only for additive!

  if integrator.opts.adaptive
    u = @. uprev + k₁/3 + 2k₂/3 + E₂ + W.dW*gpdt
    tmp = @. @muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  else
    u = @. uprev + k₁/3 + 2k₂/3 + E₂ + W.dW*gpdt
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRA1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  gpdt = integrator.g(t+dt,uprev)
  chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  k₁ = dt*integrator.f(t,uprev)
  k₂ = dt*integrator.f(t+3dt/4,uprev+3k₁/4 + 3chi2*gpdt/2)
  E₁ = k₁ + k₂
  E₂ = chi2.*(integrator.g(t,uprev)-gpdt) #Only for additive!

  if integrator.opts.adaptive
    u = uprev + k₁/3 + 2k₂/3 + E₂ + W.dW*gpdt
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  else
    u = uprev + k₁/3 + 2k₂/3 + E₂ + W.dW*gpdt
  end
  @pack integrator = t,dt,u
end

#=
@inline function perform_step!(integrator,cache::SRA1Cache,f=integrator.f)
  @unpack chi2,tmp1,E₁,E₂,gt,k₁,k₂,gpdt,tmp = cache
  @unpack t,dt,uprev,u,W = integrator
  integrator.g(t,uprev,gt)
  integrator.g(t+dt,uprev,gpdt)
  integrator.f(t,uprev,k₁); k₁*=dt
  @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  @. tmp1 = uprev+3k₁/4 + 3chi2*gpdt/2

  integrator.f(t+3dt/4,tmp1,k₂); k₂*=dt

  @. E₁ = k₁ + k₂
  @. E₂ = chi2*(gt-gpdt) #Only for additive!

  @. u = uprev + k₁/3 + 2k₂/3 + E₂ + W.dW*gpdt

  if integrator.opts.adaptive
    @. tmp = @muladd(integrator.opts.delta*E₁+E₂)/@muladd(integrator.opts.abstol + max(integrator.opts.internalnorm(uprev),integrator.opts.internalnorm(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRA1Cache,f=integrator.f)
  @unpack chi2,tmp1,E₁,E₂,gt,k₁,k₂,gpdt,tmp = cache
  @unpack t,dt,uprev,u,W = integrator
  integrator.g(t,uprev,gt)
  integrator.g(t+dt,uprev,gpdt)
  integrator.f(t,uprev,k₁); k₁*=dt
  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  if is_diagonal_noise(integrator.sol.prob)
    @. E₁ = chi2*gpdt
  else
    A_mul_B!(E₁,gpdt,chi2)
  end

  for i in eachindex(u)
    @inbounds tmp1[i] = uprev[i]+3k₁[i]/4 + 3E₁[i]/2
  end

  integrator.f(t+3dt/4,tmp1,k₂); k₂*=dt

  for i in eachindex(u)
    @inbounds E₁[i] = k₁[i] + k₂[i]
    @inbounds E₂[i] = chi2[i]*(gt[i]-gpdt[i]) #Only for additive!
  end

  if is_diagonal_noise(integrator.sol.prob)
    @. tmp1 = W.dW*gpdt
  else
    A_mul_B!(tmp1,gpdt,W.dW)
  end

  for i in eachindex(u)
    @inbounds u[i] = uprev[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + tmp1[i]
  end

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SRA2ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack a21,b21,c02,c11,c12,α1,α2,beta12,beta21,beta22 = cache

  chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h

  g1 = integrator.g(t+c11*dt,uprev)
  k1 = integrator.f(t,uprev)
  H01 = uprev + dt*a21*k1 + chi2*b21*g1

  g2 = integrator.g(t+c12*dt,H01)
  k2 = integrator.f(t+c02*dt,H01)

  E₁ = dt*(α1*k1 + α2*k2)
  E₂ = W.dW*(beta12*g2) + chi2*(beta21*g1 + beta22*g2)

  u = uprev + E₁ + E₂

  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SRA2Cache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack chi2,tab,g1,g2,k1,k2,E₁,E₂,tmp = cache
  @unpack a21,b21,c02,c11,c12,α1,α2,beta12,beta21,beta22 = cache.tab
  H01 = E₁

  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  integrator.g(t+c11*dt,uprev,g1)
  integrator.f(t,uprev,k1)

  if is_diagonal_noise(integrator.sol.prob)
    @. H01 = uprev + dt*a21*k1 + chi2*b21*g1
  else
    A_mul_B!(E₁,g1,chi2)
    @. H01 = uprev + dt*a21*k1 + b21*E₁
  end

  integrator.g(t+c12*dt,H01,g2)
  integrator.f(t+c02*dt,H01,k2)

  if is_diagonal_noise(integrator.sol.prob)
    @. E₂ = W.dW*(beta12*g2) + chi2*(beta21*g1 + beta22*g2)
  else
    @. g1 = beta21*g1 + beta22*g2
    A_mul_B!(E₁,g1,chi2)
    g2 .*= beta12
    A_mul_B!(E₂,g2,W.dW)
    @. E₂ += E₁
  end

  @. E₁ = dt*(α1*k1 + α2*k2)
  @. u = uprev + E₁ + E₂

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::ThreeStageSRAConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23 = cache

  chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  g1 = integrator.g(t+c11*dt,uprev)

  k1 = integrator.f(t,uprev)
  H01 = uprev + dt*a21*k1 + chi2*b21*g1

  g2 = integrator.g(t+c12*dt,H01)
  k2 = integrator.f(t+c02*dt,H01)

  H02 = uprev + dt*(a31*k1 + a32*k2) + chi2*(b31*g1 + b32*g2)

  g3 = integrator.g(t+c13*dt,H02)
  k3 = integrator.f(t+c03*dt,H02)

  E₁ = dt*(α1*k1 + α2*k2 + α3*k3)
  E₂ = W.dW*(beta11*g1 + beta12*g2 + beta13*g3) + chi2*(beta21*g1 + beta22*g2 + beta23*g3)

  u = uprev + E₁ + E₂

  if integrator.opts.adaptive
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::ThreeStageSRACache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack chi2,tab,g1,g2,g3,k1,k2,k3,E₁,E₂,tmp,gtmp = cache
  @unpack a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23 = cache.tab

  H01 = E₁; H02 = E₁

  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  integrator.g(t+c11*dt,uprev,g1)
  integrator.f(t,uprev,k1)

  if is_diagonal_noise(integrator.sol.prob)
    @. H01 = uprev + dt*a21*k1 + chi2*b21*g1
  else
    A_mul_B!(E₁,g1,chi2)
    @. H01 = uprev + dt*a21*k1 + b21*E₁
  end

  integrator.g(t+c12*dt,H01,g2)
  integrator.f(t+c02*dt,H01,k2)

  if is_diagonal_noise(integrator.sol.prob)
    for i in eachindex(u)
      H02[i] = uprev[i] + dt*(a31*k1[i] + a32*k2[i]) + chi2[i]*(b31*g1[i] + b32*g2[i])
    end
  else
    @. gtmp = b31*g1 + b32*g2
    A_mul_B!(E₁,gtmp,chi2)
    for i in eachindex(u)
      H02[i] = uprev[i] + dt*(a31*k1[i] + a32*k2[i]) + E₁[i]
    end
  end

  integrator.g(t+c13*dt,H02,g3)
  integrator.f(t+c03*dt,H02,k3)

  if is_diagonal_noise(integrator.sol.prob)
    for i in eachindex(u)
      E₂[i] = W.dW[i]*(beta11*g1[i] + beta12*g2[i] + beta13*g3[i]) + chi2[i]*(beta21*g1[i] + beta22*g2[i] + beta23*g3[i])
    end
  else
    @. gtmp = beta21*g1 + beta22*g2 + beta23*g3
    A_mul_B!(E₁,gtmp,chi2)
    @. gtmp = beta11*g1 + beta12*g2 + beta13*g3
    A_mul_B!(E₂,gtmp,W.dW)
    @. E₂ += E₁
  end

  @. E₁ = dt*(α1*k1 + α2*k2 + α3*k3)
  u = uprev + E₁ + E₂

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

#=
@inline function perform_step!(integrator,cache::SRACache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,tmp = cache
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages = cache.tab
  @. chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  for i in 1:stages
    fill!(H0[i],zero(eltype(integrator.u)))
  end
  for i = 1:stages
    fill!(A0temp,zero(eltype(integrator.u)))
    fill!(B0temp,zero(eltype(integrator.u)))
    for j = 1:i-1
      integrator.f(@muladd(t + c₀[j]*dt),H0[j],ftmp)
      integrator.g(@muladd(t + c₁[j]*dt),H0[j],gtmp)
      @. A0temp = @muladd A0temp + A₀[j,i]*ftmp
      @. B0temp = @muladd B0temp + B₀[j,i]*gtmp
    end
    @. H0[i] = @muladd uprev + A0temp*dt + B0temp*chi2
  end
  fill!(atemp ,zero(eltype(integrator.u)))
  fill!(btemp ,zero(eltype(integrator.u)))
  fill!(E₂    ,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))

  for i = 1:stages
    integrator.f(@muladd(t+c₀[i]*dt),H0[i],ftmp)
    integrator.g(@muladd(t+c₁[i]*dt),H0[i],gtmp)
    @. atemp  =  @muladd atemp  + α[i]*ftmp
    @. btemp  =  @muladd btemp  + (β₁[i]*W.dW)*gtmp
    @. E₂     =  @muladd E₂     + (β₂[i]*chi2)*gtmp
    @. E₁temp =  E₁temp +  ftmp
  end
  @. E₁ = dt*E₁temp
  @. u = @muladd uprev + dt*atemp + btemp + E₂

  if integrator.opts.adaptive
    @. tmp = (integrator.opts.delta*E₁+E₂)/(integrator.opts.abstol + max(integrator.opts.internalnorm(uprev),integrator.opts.internalnorm(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRACache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  @unpack H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,tmp = cache
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages = cache.tab

  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = @. (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  for i in 1:stages
    fill!(H0[i],zero(eltype(integrator.u)))
  end
  for i = 1:stages
    fill!(A0temp,zero(eltype(integrator.u)))
    fill!(B0temp,zero(eltype(integrator.u)))
    for j = 1:i-1
      integrator.f(@muladd(t + c₀[j]*dt),H0[j],ftmp)
      integrator.g(@muladd(t + c₁[j]*dt),H0[j],gtmp)
      @. A0temp = @muladd A0temp + A₀[j,i]*ftmp
      if is_diagonal_noise(integrator.sol.prob)
        @. B0temp = @muladd B0temp + B₀[j,i]*gtmp*chi2
      else
        A_mul_B!(E₁temp,gtmp,chi2)
        @. B0temp = @muladd B0temp + B₀[j,i]*E₁temp
      end
    end
    @tight_loop_macros for j in eachindex(u)
      @inbounds H0[i][j] = @muladd uprev[j] + A0temp[j]*dt + B0temp[j]
    end
  end
  fill!(atemp ,zero(eltype(integrator.u)))
  fill!(btemp ,zero(eltype(integrator.u)))
  fill!(E₂    ,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))

  for i = 1:stages
    integrator.f(@muladd(t+c₀[i]*dt),H0[i],ftmp)
    integrator.g(@muladd(t+c₁[i]*dt),H0[i],gtmp)
    if is_diagonal_noise(integrator.sol.prob)
      @. btemp = @muladd btemp + β₁[i]*W.dW*gtmp
    else
      A_mul_B!(E₁temp,gtmp,W.dW)
      @. btemp = @muladd btemp + β₁[i]*E₁temp
    end
    if is_diagonal_noise(integrator.sol.prob)
      @. E₂ = @muladd E₂ + β₂[i]*chi2*gtmp
    else
      A_mul_B!(E₁temp,gtmp,chi2)
      @. E₂ = @muladd E₂ + β₂[i]*E₁temp
    end
    @tight_loop_macros for j in eachindex(u)
      @inbounds atemp[j]  =  @muladd atemp[j]  + α[i]*ftmp[j]
      @inbounds E₁temp[j] =  E₁temp[j] +  ftmp[j]
    end
  end
  @tight_loop_macros for i in eachindex(u)
    @inbounds E₁[i] = dt*E₁temp[i]
  end

  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = @muladd uprev[i] + dt*atemp[i] + btemp[i] + E₂[i]
  end

  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end

#=
@inline function perform_step!(integrator,cache::SRAConstantCache,f=integrator.f)
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages,H0 = cache
  @unpack t,dt,uprev,u,W = integrator
  chi2 = @. .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  H0 .= zero(zero(u))
  for i = 1:stages
    A0temp = zero(u)
    B0temp = zero(u)
    for j = 1:i-1
      A0temp = @muladd A0temp .+ A₀[j,i].*integrator.f(@muladd(t .+ c₀[j].*dt),H0[j])
      B0temp = @muladd B0temp .+ B₀[j,i].*integrator.g(@muladd(t .+ c₁[j].*dt),H0[j]) #H0[..,i] argument ignored
    end
    H0[i] = @. uprev + A0temp*dt + B0temp.*chi2
  end

  atemp = zero(u)
  btemp = zero(u)
  E₂    = zero(u)
  E₁temp= zero(u)

  for i = 1:stages
    ftemp = integrator.f(t+c₀[i]*dt,H0[i])
    E₁temp =  @. E₁temp +  ftemp
    atemp  =  @. @muladd atemp  + α[i]*ftemp
    btemp  =  @muladd btemp  .+ (β₁[i].*W.dW ).* integrator.g(@muladd(t.+c₁[i]*dt),H0[i]) #H0[i] argument ignored
    E₂     =  @muladd E₂     .+ (β₂[i].*chi2) .* integrator.g(@muladd(t.+c₁[i]*dt),H0[i]) #H0[i] argument ignored
  end

  u = @. @muladd uprev + dt*atemp + btemp + E₂

  if integrator.opts.adaptive
    E₁ = @. dt*E₁temp
    tmp = @. @muladd(integrator.opts.delta*E₁+E₂)/@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm(uprev),integrator.opts.internalnorm(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  @pack integrator = t,dt,u
end
=#

@inline function perform_step!(integrator,cache::SRAConstantCache,f=integrator.f)
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages,H0 = cache
  @unpack t,dt,uprev,u,W = integrator
  chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  H0[:]=zeros(stages)
  for i = 1:stages
    A0temp = zero(u)
    B0temp = zero(u)
    for j = 1:i-1
      A0temp = @muladd A0temp + A₀[j,i]*integrator.f(@muladd(t + c₀[j]*dt),H0[j])
      B0temp = @muladd B0temp + B₀[j,i]*integrator.g(@muladd(t + c₁[j]*dt),H0[j]) #H0[..,i] argument ignored
    end
    H0[i] = uprev + A0temp*dt + B0temp.*chi2
  end

  atemp = zero(u)
  btemp = zero(u)
  E₂    = zero(u)
  E₁temp= zero(u)

  for i = 1:stages
    ftemp = integrator.f(t+c₀[i]*dt,H0[i])
    E₁temp =  E₁temp +  ftemp
    atemp  =  @muladd atemp  + α[i]*ftemp
    btemp  =  @muladd btemp  + (β₁[i]*W.dW ).* integrator.g(@muladd(t+c₁[i]*dt),H0[i]) #H0[i] argument ignored
    E₂     =  @muladd E₂     + (β₂[i]*chi2).*integrator.g(@muladd(t+c₁[i]*dt),H0[i]) #H0[i] argument ignored
  end

  u = @muladd uprev + dt*atemp + btemp + E₂

  if integrator.opts.adaptive
    E₁ = dt*E₁temp
    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  end
  @pack integrator = t,dt,u
end
