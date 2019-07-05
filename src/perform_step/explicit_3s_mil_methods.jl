#methods from https://doi.org/10.1016/j.cam.2009.11.010
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
  u  = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k₁
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

@muladd function perform_step!(integrator,cache::WangLi3SMil_BConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  #stage 1
  k   = integrator.g(uprev,p,t)
  tmp = uprev + k*integrator.sqdt
  k₁  = integrator.g(tmp,p,t)
  k₁  = (k₁ - k)/integrator.sqdt
  u  = uprev .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k₁

  #stage 2
  k  = integrator.f(u,p,t)
  u  = u + dt*k

  #stage 3
  k   = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁  = integrator.g(tmp,p,t)
  k₁  = (k₁ - k)/integrator.sqdt
  u   = u - (dt/2)*k₁
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_BCache,f=integrator.f)
  @unpack tmp,k,k₁ = cache
  @unpack t,dt,uprev,u,W,p = integrator

    #stage 1
  integrator.g(k,uprev,p,t)
  @.. tmp = uprev + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁  = (k₁ - k)/integrator.sqdt
  @.. u   = uprev + W.dW*k + (W.dW^2/2)*k₁

  #stage 2
  integrator.f(k,u,p,t)
  @.. u  = u + dt*k

  #stage 3
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁  = (k₁ - k)/integrator.sqdt
  @.. u   = u - (dt/2)*k₁

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_CConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  #stage 1
  k  = integrator.f(uprev,p,t)
  u  = uprev + dt*k

  #stage 2
  k  = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k₁

  #stage 3
  k  = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = u - (dt/2)*k₁
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_CCache,f=integrator.f)
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
  @.. u  = u + k*W.dW + (W.dW^2/2)*k₁

  #stage 3
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = u - (dt/2)*k₁
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_DConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  #stage 1
  k  = integrator.g(uprev,p,t)
  tmp = uprev + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = uprev - (dt/2)*k₁

  #stage 2
  k  = integrator.f(u,p,t)
  u  = u + dt*k

  #stage 3
  k  = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k₁
 integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_DCache,f=integrator.f)
  @unpack tmp,k,k₁ = cache
  @unpack t,dt,uprev,u,W,p = integrator

  #stage 1
  integrator.g(k,uprev,p,t)
  @.. tmp = uprev + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = uprev - (dt/2)*k₁

  #stage 2
  integrator.f(k,u,p,t)
  @.. u  = u + dt*k

  #stage 3
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = u + k*W.dW + (W.dW^2/2)*k₁
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_EConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  #stage 1
  k  = integrator.g(uprev,p,t)
  tmp = uprev + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = uprev - (dt/2)*k₁

  #stage 2
  k  = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁ = integrator.g(tmp,p,t)
  k₁ = (k₁ - k)/(integrator.sqdt)
  u  = u .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k₁

  #stage 3
  k  = integrator.f(u,p,t)
  u  = u + dt*k
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_ECache,f=integrator.f)
  @unpack tmp,k,k₁ = cache
  @unpack t,dt,uprev,u,W,p = integrator

  #stage 1
  integrator.g(k,uprev,p,t)
  @.. tmp = uprev + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = uprev - (dt/2)*k₁

  #stage 2
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁ = (k₁ - k)/(integrator.sqdt)
  @.. u  = u + k*W.dW + (W.dW^2/2)*k₁

  #stage 3
  integrator.f(k,u,p,t)
  @.. u  = u + dt*k
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_FConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  #stage 1
  k   = integrator.g(uprev,p,t)
  tmp = uprev + k*integrator.sqdt
  k₁  = integrator.g(tmp,p,t)
  k₁  = (k₁ - k)/integrator.sqdt
  u  = uprev .+ k .* W.dW .+ 0.5 .* (W.dW .^ 2) .* k₁

  #stage 2
  k   = integrator.g(u,p,t)
  tmp = u + k*integrator.sqdt
  k₁  = integrator.g(tmp,p,t)
  k₁  = (k₁ - k)/integrator.sqdt
  u   = u - (dt/2)*k₁

  #stage 3
  k  = integrator.f(u,p,t)
  u  = u + dt*k
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::WangLi3SMil_FCache,f=integrator.f)
  @unpack tmp,k,k₁ = cache
  @unpack t,dt,uprev,u,W,p = integrator

    #stage 1
  integrator.g(k,uprev,p,t)
  @.. tmp = uprev + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁  = (k₁ - k)/integrator.sqdt
  @.. u   = uprev + W.dW*k + (W.dW^2/2)*k₁

  #stage 2
  integrator.g(k,u,p,t)
  @.. tmp = u + k*integrator.sqdt
  integrator.g(k₁,tmp,p,t)
  @.. k₁  = (k₁ - k)/integrator.sqdt
  @.. u   = u - (dt/2)*k₁

  #stage 3
  integrator.f(k,u,p,t)
  @.. u  = u + dt*k
  integrator.u = u
end
