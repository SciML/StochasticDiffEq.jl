@muladd function perform_step!(integrator,cache::TauLeapingConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p,P,c = integrator
  tmp = c(uprev, p, t, P.dW, nothing)
  integrator.u = uprev .+ tmp
end

@muladd function perform_step!(integrator,cache::TauLeapingCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p,P,c = integrator
  @unpack tmp = cache
  c(tmp, uprev, p, t, P.dW, nothing)
  @.. u = uprev + tmp
end
