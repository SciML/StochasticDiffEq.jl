@inline function perform_step!(integrator,cache::EMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW = integrator
  u = @muladd uprev + dt.*integrator.f(t,uprev) + integrator.g(t,uprev).*ΔW
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::EMCache,f=integrator.f)
  @unpack rtmp1,rtmp2,rtmp3 = cache
  @unpack t,dt,uprev,u,ΔW = integrator
  integrator.f(t,uprev,rtmp1)
  integrator.g(t,uprev,rtmp2)

  if is_diagonal_noise(integrator.sol.prob)
    for i in eachindex(u)
      rtmp2[i]*=ΔW[i] # rtmp2 === rtmp3
    end
  else
    A_mul_B!(rtmp3,rtmp2,ΔW)
  end

  for i in eachindex(u)
    u[i] = @muladd uprev[i] + dt*rtmp1[i] + rtmp3[i]
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::EulerHeunConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW = integrator
  ftmp = integrator.f(t,uprev)
  gtmp = integrator.g(t,uprev)
  tmp = @muladd uprev + ftmp.*dt + gtmp.*ΔW
  u = @muladd uprev + (1/2)*dt.*(ftmp+integrator.f(t,tmp)) + (1/2)*(gtmp+integrator.g(t,tmp)).*ΔW
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::EulerHeunCache,f=integrator.f)
  @unpack rtmp1,rtmp2,rtmp3,rtmp4,tmp = cache
  @unpack t,dt,uprev,u,ΔW = integrator
  integrator.f(t,uprev,rtmp1)
  integrator.g(t,uprev,rtmp2)
  for i in eachindex(u)
    tmp[i] = @muladd uprev[i] + rtmp1[i]*dt + rtmp2[i]*ΔW[i]
  end
  integrator.f(t,tmp,rtmp3)
  integrator.g(t,tmp,rtmp4)
  dto2 = dt*(1/2)
  for i in eachindex(u)
    ΔWo2 = (1/2)*ΔW[i]
    u[i] = @muladd uprev[i] + dto2*(rtmp1[i]+rtmp3[i]) + ΔWo2*(rtmp2[i]+rtmp4[i])
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::RandomEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W = integrator
  u = @muladd uprev + dt.*integrator.f(t,uprev,W)
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::RandomEMCache,f=integrator.f)
  @unpack rtmp = cache
  @unpack t,dt,uprev,u,W = integrator
  integrator.f(t,uprev,W,rtmp)
  for i in eachindex(u)
    u[i] = @muladd uprev[i] + dt*rtmp[i]
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::RKMilConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW = integrator
  K = @muladd uprev + dt.*integrator.f(t,uprev)
  L = integrator.g(t,uprev)
  utilde = K + L*integrator.sqdt
  if alg_interpretation(integrator.alg) == :Ito
    u = @muladd K+L*ΔW+(integrator.g(t,utilde)-L)/(2integrator.sqdt)*(ΔW^2 - dt)
  elseif alg_interpretation(integrator.alg) == :Stratonovich
    u = K+L*ΔW+((integrator.g(t,utilde)-L)/(2integrator.sqdt))*(ΔW^2)
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::RKMilCache,f=integrator.f)
  @unpack du1,du2,K,tmp,L = cache
  @unpack t,dt,uprev,u,ΔW = integrator
  integrator.f(t,uprev,du1)
  integrator.g(t,uprev,L)
  for i in eachindex(u)
    K[i] = @muladd uprev[i] + dt*du1[i]
    tmp[i] = @muladd K[i] + L[i]*integrator.sqdt
  end
  integrator.g(t,tmp,du2)
  if alg_interpretation(integrator.alg) == :Ito
    for i in eachindex(u)
      u[i] = @muladd K[i]+L[i]*ΔW[i]+(du2[i]-L[i])/(2integrator.sqdt)*(ΔW[i].^2 - dt)
    end
  elseif alg_interpretation(integrator.alg) == :Stratonovich
    for i in eachindex(u)
      u[i] = @muladd K[i]+L[i]*ΔW[i]+(du2[i]-L[i])/(2integrator.sqdt)*(ΔW[i].^2)
    end
  end
  @pack integrator = t,dt,u
end
