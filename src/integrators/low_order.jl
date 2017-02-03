@inline function perform_step!(integrator,cache::EMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW = integrator
  u = @muladd uprev + dt.*integrator.f(t,uprev) + integrator.g(t,uprev).*ΔW
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::EMCache,f=integrator.f)
  @unpack tmp,utmp2 = cache
  @unpack t,dt,uprev,u,ΔW = integrator
  integrator.f(t,uprev,tmp)
  integrator.g(t,uprev,utmp2)
  for i in eachindex(u)
    u[i] = @muladd uprev[i] + dt*tmp[i] + utmp2[i]*ΔW[i]
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
  for i in eachindex(u)
    u[i] = @muladd K[i]+L[i]*ΔW[i]+(du2[i]-L[i])./(2integrator.sqdt).*(ΔW[i].^2 - dt)
  end
  @pack integrator = t,dt,u
end

@inline function perform_step!(integrator,cache::RKMilConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,ΔW = integrator
  K = @muladd uprev + dt.*integrator.f(t,uprev)
  L = integrator.g(t,uprev)
  utilde = K + L*integrator.sqdt
  u = @muladd K+L*ΔW+(integrator.g(t,utilde)-integrator.g(t,uprev))/(2integrator.sqdt)*(ΔW^2 - dt)
  @pack integrator = t,dt,u
end
