@muladd function perform_step!(integrator,cache::LambaEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  du1 = integrator.f(uprev,p,t)
  K = @muladd uprev + dt*du1

  if is_split_step(integrator.alg)
    L = integrator.g(uprev,p,t+dt)
  else
    L = integrator.g(K,p,t+dt)
  end

  mil_correction = zero(u)

  u = K+L*W.dW

  if integrator.opts.adaptive
    du2 = integrator.f(K,p,t+dt)
    Ed = dt*(du2 - du1)/2

    utilde =  K + L*integrator.sqdt
    ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
    En = ggprime.*(W.dW.^2 .- dt)./2

    integrator.EEst = integrator.opts.internalnorm((Ed + En)/((integrator.opts.abstol + max.(abs(uprev),abs(u))*integrator.opts.reltol)))
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::LambaEMCache,f=integrator.f)
  @unpack du1,du2,K,tmp,L,gtmp,dW_cache = cache
  @unpack t,dt,uprev,u,W,p = integrator

  integrator.f(du1,uprev,p,t)
  @. K = @muladd uprev + dt*du1

  if is_split_step(integrator.alg)
    integrator.g(L,K,p,t+dt)
  else
    integrator.g(L,uprev,p,t+dt)
  end

  if is_diagonal_noise(integrator.sol.prob)
    @. tmp=L*W.dW
  else
    A_mul_B!(tmp,L,W.dW)
  end

  @. u = K+tmp

  if integrator.opts.adaptive

    if !is_diagonal_noise(integrator.sol.prob)
      g_sized = norm(L,2)
    else
      g_sized = L
    end

    @. tmp = @muladd K + L*integrator.sqdt

    if !is_diagonal_noise(integrator.sol.prob)
      integrator.g(gtmp,tmp,p,t)
      g_sized2 = norm(gtmp,2)
      @. dW_cache = dW.^2 - dt
      diff_tmp = integrator.opts.internalnorm(dW_cache)
      En = (g_sized2-g_sized)/(2integrator.sqdt)*diff_tmp
      @. tmp = En
    else
      integrator.g(gtmp,tmp,p,t)
      @. tmp = (gtmp-L)/(2integrator.sqdt)*(W.dW.^2 - dt)
    end

    # Ed
    integrator.f(du2,K,p,t+dt)
    @. tmp += integrator.opts.internalnorm(dt*(du2 - du1)/2)


    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds tmp[i] = (tmp[i])/(atol + max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
end

@muladd function perform_step!(integrator,cache::LambaEulerHeunConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  du1 = integrator.f(uprev,p,t)
  K = @muladd uprev + dt*du1
  L = integrator.g(uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    noise = L.*W.dW
  else
    noise = L*W.dW
  end
  tmp = @muladd K+L*W.dW
  gtmp2 = (1/2).*(L.+integrator.g(tmp,p,t+dt))
  if is_diagonal_noise(integrator.sol.prob)
    noise2 = gtmp2.*W.dW
  else
    noise2 = gtmp2*W.dW
  end

  u = @muladd uprev + (1/2)*dt*(du1+integrator.f(tmp,p,t+dt)) + noise2

  if integrator.opts.adaptive
    du2 = integrator.f(K,p,t+dt)
    Ed = dt*(du2 - du1)/2

    utilde = uprev + L*integrator.sqdt
    ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
    En = ggprime.*(W.dW.^2)./2

    integrator.EEst = integrator.opts.internalnorm((Ed + En)/((integrator.opts.abstol + max.(abs(uprev),abs(u))*integrator.opts.reltol)))
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::LambaEulerHeunCache,f=integrator.f)
  @unpack du1,du2,K,tmp,L,gtmp,dW_cache = cache
  @unpack t,dt,uprev,u,W,p = integrator
  integrator.f(du1,uprev,p,t)
  integrator.g(L,uprev,p,t)
  @. K = @muladd uprev + dt*du1

  if is_diagonal_noise(integrator.sol.prob)
    @. tmp=L*W.dW
  else
    A_mul_B!(tmp,L,W.dW)
  end

  @. tmp = K+tmp

  integrator.f(du2,tmp,p,t+dt)
  integrator.g(gtmp,tmp,p,t+dt)

  if is_diagonal_noise(integrator.sol.prob)
    #@. nrtmp=(1/2)*W.dW*(gtmp1+gtmp2)
    @tight_loop_macros for i in eachindex(u)
      @inbounds dWo2 = (1/2)*W.dW[i]
      @inbounds tmp[i]=dWo2*(L[i]+gtmp[i])
    end
  else
    #@. gtmp1 = (1/2)*(gtmp1+gtmp2)
    @tight_loop_macros for i in eachindex(gtmp)
      @inbounds gtmp[i] = (1/2)*(L[i]+gtmp[i])
    end
    A_mul_B!(tmp,gtmp,W.dW)
  end

  dto2 = dt*(1/2)
  #@. u = @muladd uprev + dto2*(ftmp1+ftmp2) + nrtmp
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = @muladd uprev[i] + dto2*(du1[i]+du2[i]) + tmp[i]
  end

  if integrator.opts.adaptive

    if !is_diagonal_noise(integrator.sol.prob)
      g_sized = norm(L,2)
    else
      g_sized = L
    end

    @. tmp = @muladd uprev + L*integrator.sqdt

    if !is_diagonal_noise(integrator.sol.prob)
      integrator.g(gtmp,tmp,p,t)
      g_sized2 = norm(gtmp,2)
      @. dW_cache = dW.^2
      diff_tmp = integrator.opts.internalnorm(dW_cache)
      En = (g_sized2-g_sized)/(2integrator.sqdt)*diff_tmp
      @. tmp = En
    else
      integrator.g(gtmp,tmp,p,t)
      @. tmp = (gtmp-L)/(2integrator.sqdt)*(W.dW.^2)
    end

    # Ed
    integrator.f(du2,K,p,t+dt)
    @. tmp += integrator.opts.internalnorm(dt*(du2 - du1)/2)


    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds tmp[i] = (tmp[i])/(atol + max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
end
