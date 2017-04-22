@inline function perform_step!(integrator,cache::SplitEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  u = @muladd uprev + dt.*(integrator.f[1](t,uprev) +
                           integrator.f[2](t,uprev))+
                           integrator.g(t,uprev).*W.dW
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::SplitEMCache,f=integrator.f)
  @unpack rtmp1,rtmp2,rtmp3 = cache
  @unpack t,dt,uprev,u,W = integrator

  integrator.g(t,uprev,rtmp2)
  if is_diagonal_noise(integrator.sol.prob)
    for i in eachindex(u)
      rtmp2[i]*=W.dW[i] # rtmp2 === rtmp3
    end
  else
    A_mul_B!(rtmp3,rtmp2,W.dW)
  end

  integrator.f[1](t,uprev,rtmp1)
  for i in eachindex(u)
    u[i] = @muladd uprev[i] + dt*rtmp1[i] + rtmp3[i]
  end
  integrator.f[2](t,uprev,rtmp1)
  for i in eachindex(u)
    u[i] = @muladd u[i] + dt*rtmp1[i]
  end

  @pack integrator = t,dt,u
end
