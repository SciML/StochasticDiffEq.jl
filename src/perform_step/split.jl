@muladd function perform_step!(integrator,cache::SplitEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  u = dt*(integrator.f.f1(uprev,p,t) +
           integrator.f.f2(uprev,p,t)) +
           integrator.g(uprev,p,t).*W.dW + uprev
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SplitEMCache,f=integrator.f)
  @unpack rtmp1,rtmp2 = cache
  @unpack t,dt,uprev,u,W,p = integrator

  integrator.g(rtmp2,uprev,p,t)
  if is_diagonal_noise(integrator.sol.prob)
    rmul!(rtmp2,W.dW) # rtmp2 === rtmp3
    @.. u = uprev + rtmp2
  else
    mul!(rtmp1,rtmp2,W.dW)
    @.. u = uprev + rtmp1
  end

  integrator.f.f1(t,uprev,rtmp1)
  @.. u += dt*rtmp1
  integrator.f.f2(t,uprev,rtmp1)
  @.. u += dt*rtmp1
end
