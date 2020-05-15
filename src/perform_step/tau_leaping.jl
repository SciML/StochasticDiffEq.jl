@muladd function perform_step!(integrator,cache::TauLeapingConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  K = uprev .+ integrator.P.dP
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::TauLeapingCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @.. u = uprev + dt * integrator.P.dP
end
