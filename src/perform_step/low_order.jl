@muladd function perform_step!(integrator,cache::EMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator

  K = uprev .+ dt .* integrator.f(uprev,p,t)

  if is_split_step(integrator.alg)
    u_choice = K
  else
    u_choice = uprev
  end

  if !is_diagonal_noise(integrator.sol.prob) || typeof(W.dW) <: Number
    noise = integrator.g(u_choice,p,t)*W.dW
  else
    noise = integrator.g(u_choice,p,t).*W.dW
  end

  u = K + noise
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::EMCache,f=integrator.f)
  @unpack rtmp1,rtmp2 = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(rtmp1,uprev,p,t)

  @.. u = uprev + dt * rtmp1

  if is_split_step(integrator.alg)
    u_choice = u
  else
    u_choice = uprev
  end

  integrator.g(rtmp2,u_choice,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. rtmp2 *= W.dW
    @.. u += rtmp2
  else
    mul!(rtmp1,rtmp2,W.dW)
    @.. u += rtmp1
  end
end

@muladd function perform_step!(integrator,cache::EulerHeunConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  ftmp = integrator.f(uprev,p,t)
  gtmp = integrator.g(uprev,p,t)
  if is_diagonal_noise(integrator.sol.prob)
    noise = gtmp.*W.dW
  else
    noise = gtmp*W.dW
  end
  tmp = @.. uprev + dt * ftmp + noise
  gtmp2 = (gtmp .+ integrator.g(tmp,p,t+dt)) ./ 2
  if is_diagonal_noise(integrator.sol.prob)
    noise2 = gtmp2.*W.dW
  else
    noise2 = gtmp2*W.dW
  end
  u = uprev .+ (dt / 2) .* (ftmp .+ integrator.f(tmp,p,t+dt)) .+ noise2
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::EulerHeunCache,f=integrator.f)
  @unpack ftmp1,ftmp2,gtmp1,gtmp2,tmp,nrtmp = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(ftmp1,uprev,p,t)
  integrator.g(gtmp1,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. nrtmp=gtmp1*W.dW
  else
    mul!(nrtmp,gtmp1,W.dW)
  end

  @.. tmp = uprev + dt * ftmp1 + nrtmp

  integrator.f(ftmp2,tmp,p,t+dt)
  integrator.g(gtmp2,tmp,p,t+dt)

  if is_diagonal_noise(integrator.sol.prob)
    @.. nrtmp=(1/2)*W.dW*(gtmp1+gtmp2)
  else
    @.. gtmp1 = (1/2)*(gtmp1+gtmp2)
    mul!(nrtmp,gtmp1,W.dW)
  end

  dto2 = dt / 2
  @.. u = uprev + dto2 * (ftmp1 + ftmp2) + nrtmp
end

@muladd function perform_step!(integrator,cache::RandomEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  u = uprev .+ dt .* integrator.f(uprev,p,t,W.curW)
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::RandomEMCache,f=integrator.f)
  @unpack rtmp = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(rtmp,uprev,p,t,W.curW)
  @.. u = uprev + dt * rtmp
end

@muladd function perform_step!(integrator,cache::RKMilConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  du1 = integrator.f(uprev,p,t)
  K = @.. uprev + dt * du1
  L = integrator.g(uprev,p,t)
  mil_correction = zero(u)
  if alg_interpretation(integrator.alg) == :Ito
    utilde =  K + L*integrator.sqdt
    ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
    mil_correction = ggprime.*(W.dW.^2 .- dt)./2
  elseif alg_interpretation(integrator.alg) == :Stratonovich
    utilde = uprev + L*integrator.sqdt
    ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
    mil_correction = ggprime.*(W.dW.^2)./2
  end
  u = K+L.*W.dW+mil_correction
  if integrator.opts.adaptive
    du2 = integrator.f(K,p,t+dt)
    Ed = dt*(du2 - du1)/2
    En = W.dW.^3 .* ((du2-L)/(integrator.sqdt)).^2 / 6

    resids = calculate_residuals(Ed, En, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end
  integrator.u = u
end

#=
@muladd function perform_step!(integrator,cache::RKMilCache,f=integrator.f)
  @unpack du1,du2,K,tmp,L = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(du1,uprev,p,t)
  integrator.g(L,uprev,p,t)
  @.. K = uprev + dt * du1
  @.. tmp = K + integrator.sqdt * L
  integrator.g(du2,tmp,p,t)
  if alg_interpretation(integrator.alg) == :Ito
    @.. tmp = (du2-L)/(2integrator.sqdt)*(W.dW^2 - dt)
  elseif alg_interpretation(integrator.alg) == :Stratonovich
    @.. tmp = (du2-L)/(2integrator.sqdt)*(W.dW^2)
  end
  @.. u = K+L*W.dW + tmp
  if integrator.opts.adaptive
    @.. tmp = (tmp)/(integrator.opts.abstol + max(integrator.opts.internalnorm(uprev,t),integrator.opts.internalnorm(u,t))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp,t)
  end
  integrator.u = u
end
=#

@muladd function perform_step!(integrator,cache::RKMilCache,f=integrator.f)
  @unpack du1,du2,K,tmp,L = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(du1,uprev,p,t)
  integrator.g(L,uprev,p,t)
  @.. K = uprev + dt * du1
  @.. du2 = zero(eltype(u)) # This makes it safe to re-use the array
  if alg_interpretation(integrator.alg) == :Ito
    @.. tmp = K + integrator.sqdt * L
    integrator.g(du2,tmp,p,t)
    @.. tmp = (du2-L)/(2integrator.sqdt)*(W.dW.^2 - dt)
  elseif alg_interpretation(integrator.alg) == :Stratonovich
    @.. tmp = uprev + integrator.sqdt * L
    integrator.g(du2,tmp,p,t)
    @.. tmp = (du2-L)/(2integrator.sqdt)*(W.dW.^2)
  end
  @.. u = K+L*W.dW + tmp
  if integrator.opts.adaptive
    @.. tmp = integrator.opts.internalnorm(W.dW^3,t)*
             integrator.opts.internalnorm((du2-L)/(integrator.sqdt),t)^2 / 6
    integrator.f(du2,K,p,t+dt)
    @.. tmp += integrator.opts.internalnorm(integrator.opts.delta * dt * (du2 - du1)/2,t)

    calculate_residuals!(tmp, tmp, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
end

@muladd function perform_step!(integrator,cache::RKMilCommuteCache,f=integrator.f)
  @unpack du1,du2,K,gtmp,L = cache
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack I,mil_correction,Kj,Dgj,tmp = cache
  dW = W.dW; sqdt = integrator.sqdt
  f = integrator.f; g = integrator.g

  ggprime_norm = 0.0

  @.. mil_correction = zero(u)
  for i=1:length(dW),j=1:length(dW)
      I[j,i] = 0.5*dW[i]*dW[j]
      if alg_interpretation(integrator.alg) == :Ito
        j == i && (I[i,i] -= 0.5*dt) # Ito correction
      end
  end

  integrator.f(du1,uprev,p,t)
  integrator.g(L,uprev,p,t)

  for j = 1:length(uprev)
    @.. Kj = uprev + dt*du1 + sqdt*@view(L[:,j]) # This works too
    #Kj .= uprev .+ sqdt*L[:,j]
    g(gtmp,Kj,p,t)
    @.. Dgj = (gtmp - L)/sqdt
    if integrator.opts.adaptive
        ggprime_norm += integrator.opts.internalnorm(Dgj,t)
    end
    mul!(tmp,Dgj,@view(I[:,j]))
    mil_correction .+= tmp
  end
  mul!(tmp,L,dW)
  @.. u .= uprev + dt*du1 + tmp + mil_correction

  if integrator.opts.adaptive
      En = integrator.opts.internalnorm(W.dW,t)^3*ggprime_norm^2 / 6
      integrator.f(du2,K,p,t+dt)
      @.. tmp = integrator.opts.internalnorm(integrator.opts.delta * dt * (du2 - du1) / 2,t) + En

      calculate_residuals!(tmp, tmp, uprev, u, integrator.opts.abstol,
                           integrator.opts.reltol, integrator.opts.internalnorm, t)
      integrator.EEst = integrator.opts.internalnorm(tmp,t)

  end
end
