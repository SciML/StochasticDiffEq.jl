@inline function perform_step!(integrator,cache::SplitEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  u = dt.*(integrator.f[1](t,uprev) .+
           integrator.f[2](t,uprev)) .+
           integrator.g(t,uprev).*W.dW .+ uprev
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SplitEMCache,f=integrator.f)
  @unpack rtmp1,rtmp2,rtmp3 = cache
  @unpack t,dt,uprev,u,W = integrator

  integrator.g(t,uprev,rtmp2)
  if is_diagonal_noise(integrator.sol.prob)
    rtmp2 .*=W.dW # rtmp2 === rtmp3
  else
    A_mul_B!(rtmp3,rtmp2,W.dW)
  end

  integrator.f[1](t,uprev,rtmp1)
  @. u = @muladd uprev + dt*rtmp1 + rtmp3
  integrator.f[2](t,uprev,rtmp1)
  @. u = @muladd u + dt*rtmp1

  @pack integrator = t,dt,u
end
