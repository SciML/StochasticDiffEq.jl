@muladd function perform_step!(integrator,cache::WangLi3SMil_AConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  #stage 1
  k  = integrator.f(uprev,p,t)
  u  = uprev + dt*k

  #stage 2
  k  = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = u - (dt/2)*k₁

  #stage 3
  k  = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = u + k*W.dW + (W.dW^2/2)*k₁
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_ACache,f=integrator.f)
  @unpack tmp,k,k₁ = cache
  @unpack t,dt,uprev,u,W,p = integrator

  #stage 1
  integrator.f(k,uprev,p,t)
  @.. u  = uprev + dt*k

  #stage 2
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = u - (dt/2)*k₁

  #stage 3
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = u + k*W.dW + (W.dW^2/2)*k₁
  integrator.u = u
end
