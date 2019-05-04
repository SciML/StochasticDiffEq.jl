@muladd function perform_step!(integrator,cache::SRA1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  gpdt = integrator.g(uprev,p,t+dt)
  chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  k₁ = dt*integrator.f(uprev,p,t)
  k₂ = dt*integrator.f(uprev+3k₁/4 + 3chi2.*gpdt/2,p,t+3dt/4)
  E₂ = chi2.*(integrator.g(uprev,p,t)-gpdt) #Only for additive!
  u = @.. uprev + k₁ / 3 + 2 * k₂ / 3 + E₂ + W.dW * gpdt

  if integrator.opts.adaptive
    E₁ = k₁ .+ k₂

    resids = calculate_residuals(E₁, E₂, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end
  integrator.u = u
end

#=
@muladd function perform_step!(integrator,cache::SRA1Cache,f=integrator.f)
  @unpack chi2,tmp1,E₁,E₂,gt,k₁,k₂,gpdt,tmp = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.g(t,uprev,gt)
  integrator.g(t+dt,uprev,gpdt)
  integrator.f(t,uprev,k₁); k₁*=dt
  @.. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  @.. tmp1 = uprev+3k₁/4 + 3chi2*gpdt/2

  integrator.f(t+3dt/4,tmp1,k₂); k₂*=dt

  @.. E₂ = chi2*(gt-gpdt) #Only for additive!

  @.. u = uprev + k₁/3 + 2k₂/3 + E₂ + W.dW*gpdt

  if integrator.opts.adaptive
    @.. E₁ = k₁ + k₂

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
  integrator.u = u
end
=#

@muladd function perform_step!(integrator,cache::SRA1Cache,f=integrator.f)
  @unpack chi2,tmp1,E₁,E₂,gt,k₁,k₂,gpdt,tmp = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.g(gt,uprev,p,t)
  integrator.g(gpdt,uprev,p,t+dt)
  integrator.f(k₁,uprev,p,t); k₁*=dt
  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @.. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  if is_diagonal_noise(integrator.sol.prob)
    @.. E₁ = chi2*gpdt
  else
    mul!(E₁,gpdt,chi2)
  end

  for i in eachindex(u)
    @inbounds tmp1[i] = uprev[i]+3k₁[i]/4 + 3E₁[i]/2
  end

  integrator.f(k₂,tmp1,p,t+3dt/4); k₂*=dt

  for i in eachindex(u)
    @inbounds E₂[i] = chi2[i]*(gt[i]-gpdt[i]) #Only for additive!
  end

  if is_diagonal_noise(integrator.sol.prob)
    @.. tmp1 = W.dW*gpdt
  else
    mul!(tmp1,gpdt,W.dW)
  end

  for i in eachindex(u)
    @inbounds u[i] = uprev[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + tmp1[i]
  end

  if integrator.opts.adaptive
    @.. E₁ = k₁ + k₂

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
end

@muladd function perform_step!(integrator,cache::SRA2ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack a21,b21,c02,c11,c12,α1,α2,beta12,beta21,beta22 = cache

  chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h

  g1 = integrator.g(uprev,p,t+c11*dt)
  k1 = integrator.f(uprev,p,t)
  H01 = uprev + dt*a21*k1 + chi2*b21*g1

  g2 = integrator.g(H01,p,t+c12*dt)
  k2 = integrator.f(H01,p,t+c02*dt)

  E₁ = dt*(α1*k1 + α2*k2)
  E₂ = chi2*(beta21*g1 + beta22*g2)

  u = uprev + E₁ + E₂ + W.dW*(beta12*g2)

  if integrator.opts.adaptive
    E₁ = @.. dt * (k1 + k2)

    resids = calculate_residuals(E₁, E₂, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SRA2Cache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack chi2,tab,g1,g2,k1,k2,E₁,E₂,tmp = cache
  @unpack a21,b21,c02,c11,c12,α1,α2,beta12,beta21,beta22 = cache.tab
  H01 = E₁

  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @.. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  integrator.g(g1,uprev,p,t+c11*dt)
  integrator.f(k1,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. H01 = uprev + dt*a21*k1 + chi2*b21*g1
  else
    mul!(E₁,g1,chi2)
    @.. H01 = uprev + dt*a21*k1 + b21*E₁
  end

  integrator.g(g2,H01,p,t+c12*dt)
  integrator.f(k2,H01,p,t+c02*dt)

  if is_diagonal_noise(integrator.sol.prob)
    @.. E₂ = chi2*(beta21*g1 + beta22*g2)
    #@.. u = uprev + dt*(α1*k1 + α2*k2) + E₂ + W.dW*(beta12*g2)
    for i in eachindex(u)
      @inbounds u[i] = uprev[i] + dt*(α1*k1[i] + α2*k2[i]) + E₂[i] + W.dW[i]*(beta12*g2[i])
    end
  else
    @.. g1 = beta21*g1 + beta22*g2
    mul!(E₂,g1,chi2)
    g2 .*= beta12
    mul!(E₁,g2,W.dW)
    for i in eachindex(u)
      @inbounds u[i] = uprev[i] + dt*(α1*k1[i] + α2*k2[i]) + E₂[i] + E₁[i]
    end
  end

  if integrator.opts.adaptive
    @.. E₁ = dt * (k1 + k2)

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
end

@muladd function perform_step!(integrator,cache::ThreeStageSRAConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23 = cache

  chi2 = .5*(W.dW .+ W.dZ/sqrt(3)) #I_(1,0)/h

  g1 = integrator.g(uprev,p,t+c11*dt)
  k1 = integrator.f(uprev,p,t)

  H01 = uprev + dt*a21*k1 + b21*chi2.*g1

  g2 = integrator.g(H01,p,t+c12*dt)
  k2 = integrator.f(H01,p,t+c02*dt)

  H02 = uprev + dt*(a31*k1 + a32*k2) + chi2.*(b31*g1 + b32*g2)

  g3 = integrator.g(H02,p,t+c13*dt)
  k3 = integrator.f(H02,p,t+c03*dt)

  E₁ = dt*(α1*k1 + α2*k2 + α3*k3)
  E₂ = chi2.*(beta21*g1 + beta22*g2 + beta23*g3)

  u = uprev + E₁ + E₂ + W.dW.*(beta11*g1 + beta12*g2 + beta13*g3)

  if typeof(integrator.alg) <: StochasticCompositeAlgorithm && typeof(integrator.alg.algs[1]) <: SOSRA2
    ϱu = integrator.opts.internalnorm(k3 - k2, t)
    ϱd = integrator.opts.internalnorm(H02 - H01, t)
    integrator.eigen_est = ϱu/ϱd
  end

  if integrator.opts.adaptive
    E₁ = @.. dt * (k1 + k2 + k3)

    resids = calculate_residuals(E₁, E₂, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::ThreeStageSRACache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack chi2,tab,g1,g2,g3,k1,k2,k3,E₁,E₂,tmp,gtmp = cache
  @unpack a21,a31,a32,b21,b31,b32,c02,c03,c11,c12,c13,α1,α2,α3,beta11,beta12,beta13,beta21,beta22,beta23 = cache.tab

  H01 = E₁; H02 = E₁

  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @.. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  integrator.g(g1,uprev,p,t+c11*dt)
  integrator.f(k1,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. H01 = uprev + dt*a21*k1 + chi2*b21*g1
  else
    mul!(E₁,g1,chi2)
    @.. H01 = uprev + dt*a21*k1 + b21*E₁
  end

  integrator.g(g2,H01,p,t+c12*dt)
  integrator.f(k2,H01,p,t+c02*dt)

  if is_diagonal_noise(integrator.sol.prob)
    for i in eachindex(u)
      H02[i] = uprev[i] + dt*(a31*k1[i] + a32*k2[i]) + chi2[i]*(b31*g1[i] + b32*g2[i])
    end
  else
    @.. gtmp = b31*g1 + b32*g2
    mul!(E₁,gtmp,chi2)
    for i in eachindex(u)
      H02[i] = uprev[i] + dt*(a31*k1[i] + a32*k2[i]) + E₁[i]
    end
  end

  integrator.g(g3,H02,p,t+c13*dt)
  integrator.f(k3,H02,p,t+c03*dt)

  if is_diagonal_noise(integrator.sol.prob)
    @.. E₂ = chi2*(beta21*g1 + beta22*g2 + beta23*g3)
    #@.. u = uprev + dt*(α1*k1 + α2*k2 + α3*k3) + E₂ + W.dW*(beta11*g1 + beta12*g2 + beta13*g3)
    for i in eachindex(u)
      @inbounds u[i] = uprev[i] + dt*(α1*k1[i] + α2*k2[i] + α3*k3[i]) + E₂[i] + W.dW[i]*(beta11*g1[i] + beta12*g2[i] + beta13*g3[i])
    end
  else
    @.. gtmp = beta21*g1 + beta22*g2 + beta23*g3
    mul!(E₂,gtmp,chi2)
    @.. gtmp = beta11*g1 + beta12*g2 + beta13*g3
    mul!(E₁,gtmp,W.dW)
    for i in eachindex(u)
      @inbounds u[i] = uprev[i] + dt*(α1*k1[i] + α2*k2[i] + α3*k3[i]) + E₂[i] + E₁[i]
    end
  end

  if typeof(integrator.alg) <: StochasticCompositeAlgorithm && typeof(integrator.alg.algs[1]) <: SOSRA2
    @.. tmp = k3 - k2
    ϱu = integrator.opts.internalnorm(tmp, t)
    @.. tmp = H02 - H01
    ϱd = integrator.opts.internalnorm(tmp, t)
    integrator.eigen_est = ϱu/ϱd
  end

  if integrator.opts.adaptive
    @.. E₁ = dt * (k1 + k2 + k3)

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
end

#=
@muladd function perform_step!(integrator,cache::SRACache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,tmp = cache
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages = cache.tab
  @.. chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  for i in 1:stages
    fill!(H0[i],zero(eltype(integrator.u)))
  end
  for i = 1:stages
    fill!(A0temp,zero(eltype(integrator.u)))
    fill!(B0temp,zero(eltype(integrator.u)))
    for j = 1:i-1
      integrator.f((t + c₀[j]*dt),H0[j],ftmp)
      integrator.g((t + c₁[j]*dt),H0[j],gtmp)
      @.. A0temp = A0temp + A₀[j,i] * ftmp
      @.. B0temp = B0temp + B₀[j,i] * gtmp
    end
    @.. H0[i] = uprev + dt * A0temp + chi2 * B0temp
  end
  fill!(atemp ,zero(eltype(integrator.u)))
  fill!(btemp ,zero(eltype(integrator.u)))
  fill!(E₂    ,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))

  for i = 1:stages
    integrator.f((t+c₀[i]*dt),H0[i],ftmp)
    integrator.g((t+c₁[i]*dt),H0[i],gtmp)
    @.. atemp  =  atemp  + α[i] * ftmp
    @.. btemp  =  btemp  + (β₁[i] * W.dW) * gtmp
    @.. E₂     =  E₂     + (β₂[i] * chi2) * gtmp
    @.. E₁temp =  E₁temp +  ftmp
  end
  @.. u = uprev + dt * atemp + btemp + E₂

  if integrator.opts.adaptive
    @.. E₁ = dt*E₁temp

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
  integrator.u = u
end
=#

@muladd function perform_step!(integrator,cache::SRACache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack H0,A0temp,B0temp,ftmp,gtmp,chi2,atemp,btemp,E₁,E₁temp,E₂,tmp = cache
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages = cache.tab

  if typeof(W.dW) <: Union{SArray,Number}
    chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @.. chi2 = (W.dW + W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  for i in 1:stages
    fill!(H0[i],zero(eltype(integrator.u)))
  end
  for i = 1:stages
    fill!(A0temp,zero(eltype(integrator.u)))
    fill!(B0temp,zero(eltype(integrator.u)))
    for j = 1:i-1
      integrator.f(ftmp,H0[j],p,t + c₀[j]*dt)
      integrator.g(gtmp,H0[j],p,t + c₁[j]*dt)
      @.. A0temp = A0temp + A₀[j,i] * ftmp
      if is_diagonal_noise(integrator.sol.prob)
        @.. B0temp = B0temp + B₀[j,i] * chi2 * gtmp
      else
        mul!(E₁temp,gtmp,chi2)
        @.. B0temp = B0temp + B₀[j,i] * E₁temp
      end
    end

    @.. H0[i] = uprev + dt * A0temp + B0temp
  end
  fill!(atemp ,zero(eltype(integrator.u)))
  fill!(btemp ,zero(eltype(integrator.u)))
  fill!(E₂    ,zero(eltype(integrator.u)))
  fill!(E₁temp,zero(eltype(integrator.u)))

  for i = 1:stages
    integrator.f(ftmp,H0[i],p,t+c₀[i]*dt)
    integrator.g(gtmp,H0[i],p,t+c₁[i]*dt)
    if is_diagonal_noise(integrator.sol.prob)
      @.. btemp = btemp + β₁[i]*W.dW*gtmp
    else
      mul!(E₁temp,gtmp,W.dW)
      @.. btemp = btemp + β₁[i]*E₁temp
    end
    if is_diagonal_noise(integrator.sol.prob)
      @.. E₂ = E₂ + β₂[i]*chi2*gtmp
    else
      mul!(E₁temp,gtmp,chi2)
      @.. E₂ = E₂ + β₂[i]*E₁temp
    end

    @.. atemp  =  atemp  + α[i]*ftmp
    @.. E₁temp =  E₁temp +  ftmp
  end

  @.. u = uprev + dt*atemp + btemp + E₂

  if integrator.opts.adaptive
    @.. E₁ = dt * E₁temp

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
end

@muladd function perform_step!(integrator,cache::SRAConstantCache,f=integrator.f)
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂,stages,H0 = cache
  @unpack t,dt,uprev,u,W,p = integrator
  chi2 = .5*(W.dW + W.dZ/sqrt(3)) #I_(1,0)/h
  H0[:]=fill(zero(u),stages)

  for i = 1:stages
    A0temp = zero(u)
    B0temp = zero(u)
    for j = 1:i-1
      A0temp = A0temp .+ A₀[j,i] .* integrator.f(H0[j],p,t + c₀[j]*dt)
      B0temp = B0temp .+ B₀[j,i] .* integrator.g(H0[j],p,t + c₁[j]*dt) #H0[..,i] argument ignored
    end

    H0[i] = uprev + A0temp*dt + B0temp.*chi2
  end

  atemp = zero(u)
  btemp = zero(u)
  E₂    = zero(u)
  E₁temp= zero(u)

  for i = 1:stages
    ftmp = integrator.f(H0[i],p,t+c₀[i]*dt)
    E₁temp = E₁temp .+ ftmp
    atemp  = @.. atemp + α[i] * ftmp

    gtmp = integrator.g(H0[i], p, t + c₁[i] * dt) #H0[i] argument ignored
    btemp  = @.. btemp + (β₁[i] * W.dW) * gtmp
    E₂     = @.. E₂ + (β₂[i] * chi2) * gtmp
  end

  u = @.. uprev + dt * atemp + btemp + E₂

  if integrator.opts.adaptive
    E₁ = @.. dt * E₁temp

    resids = calculate_residuals(E₁, E₂, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end
  integrator.u = u
end
