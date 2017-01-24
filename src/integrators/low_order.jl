function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{EM,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @inbounds while integrator.t < integrator.T
    @sde_loopheader
    @unpack t,dt,uprev,u = integrator
    u = uprev + dt.*integrator.f(t,uprev) + integrator.g(t,uprev).*ΔW
    @pack integrator = t,dt,u
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{EM,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @inbounds while integrator.t < integrator.T
    @sde_loopheader
    @unpack utmp1,utmp2 = cache
    @unpack t,dt,uprev,u = integrator
    integrator.f(t,uprev,utmp1)
    integrator.g(t,uprev,utmp2)
    for i in eachindex(u)
      u[i] = uprev[i] + dt*utmp1[i] + utmp2[i]*ΔW[i]
    end
    @pack integrator = t,dt,u
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{RKMil,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @inbounds while integrator.t < integrator.T
    @sde_loopheader
    @unpack du1,du2,K,utilde,L = cache
    @unpack t,dt,uprev,u = integrator
    integrator.f(t,uprev,du1)
    integrator.g(t,uprev,L)
    for i in eachindex(u)
      K[i] = uprev[i] + dt*du1[i]
      utilde[i] = K[i] + L[i]*integrator.sqdt
    end
    integrator.g(t,utilde,du2)
    for i in eachindex(u)
      u[i] = K[i]+L[i]*ΔW[i]+(du2[i]-L[i])./(2integrator.sqdt).*(ΔW[i].^2 - dt)
    end
    @pack integrator = t,dt,u
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{RKMil,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @inbounds while integrator.t < integrator.T
    @sde_loopheader
    @unpack t,dt,uprev,u = integrator
    K = uprev + dt.*integrator.f(t,uprev)
    L = integrator.g(t,uprev)
    utilde = K + L*integrator.sqdt
    u = K+L*ΔW+(integrator.g(t,utilde)-integrator.g(t,uprev))/(2integrator.sqdt)*(ΔW^2 - dt)
    @pack integrator = t,dt,u
    @sde_loopfooter
  end
  @sde_postamble
end
