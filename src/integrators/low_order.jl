function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{EM,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  @inbounds while t<T
    @sde_loopheader
    u = u + dt.*f(t,u) + g(t,u).*ΔW
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{EM,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  utmp1 = zeros(u); utmp2 = zeros(u)
  @inbounds while t<T
    @sde_loopheader
    f(t,u,utmp1)
    g(t,u,utmp2)
    for i in eachindex(u)
      u[i] = u[i] + dt*utmp1[i] + utmp2[i]*ΔW[i]
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{RKMil,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  du1::uType = zeros(u); du2::uType = zeros(u)
  K::uType = zeros(u); utilde::uType = zeros(u); L::uType = zeros(u)
  @inbounds while t<T
    @sde_loopheader
    f(t,u,du1)
    g(t,u,L)
    for i in eachindex(u)
      K[i] = u[i] + dt*du1[i]
      utilde[i] = K[i] + L[i]*sqdt
    end
    g(t,utilde,du2)
    for i in eachindex(u)
      u[i] = K[i]+L[i]*ΔW[i]+(du2[i]-L[i])./(2sqdt).*(ΔW[i].^2 - dt)
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5}(integrator::SDEIntegrator{RKMil,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,F,F2,F3,F4,F5})
  @sde_preamble
  local L::uType; local K::uType; local utilde::uType
  @inbounds while t<T
    @sde_loopheader

    K = u + dt.*f(t,u)
    L = g(t,u)
    utilde = K + L*sqdt
    u = K+L*ΔW+(g(t,utilde)-g(t,u))/(2sqdt)*(ΔW^2 - dt)

    @sde_loopfooter
  end
  @sde_postamble
end
