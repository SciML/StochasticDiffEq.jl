function sde_solve{algType<:SRA1,uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  H0 = Array{uEltype}(size(u)...,2)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    gpdt = integrator.g(t+dt,u)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    k₁ = dt*integrator.f(t,u)
    k₂ = dt*integrator.f(t+3dt/4,u+3k₁/4 + 3chi2*integrator.g(t+dt,u)/2)
    E₁ = k₁ + k₂
    E₂ = chi2.*(integrator.g(t,u)-gpdt) #Only for additive!

    if integrator.opts.adaptive
      utmp = u + k₁/3 + 2k₂/3 + E₂ + ΔW*gpdt
      EEst = abs((integrator.opts.delta*E₁+E₂)./(integrator.opts.abstol + max(abs(u),abs(utmp))*integrator.opts.reltol))
    else
      u = u + k₁/3 + 2k₂/3 + E₂ + ΔW*gpdt
    end

    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRA1,uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble

  H0 = Array{uEltype}(size(u)...,2)
  chi2::randType = similar(ΔW)
  tmp1::uType = zeros(u)
  E₁::uType = zeros(u); gt::uType = zeros(u); gpdt::uType = zeros(u)
  E₂::uType = zeros(u); k₁::uType = zeros(u); k₂::uType = zeros(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    integrator.g(t,u,gt)
    integrator.g(t+dt,u,gpdt)
    integrator.f(t,u,k₁); k₁*=dt
    for i in eachindex(u)
      chi2[i] = (ΔW[i] + ΔZ[i]/sqrt(3))/2 #I_(1,0)/h
      tmp1[i] = u[i]+3k₁[i]/4 + 3chi2[i]*gpdt[i]/2
    end

    integrator.f(t+3dt/4,tmp1,k₂); k₂*=dt

    for i in eachindex(u)
      E₁[i] = k₁[i] + k₂[i]
      E₂[i] = chi2[i]*(gt[i]-gpdt[i]) #Only for additive!
    end

    if integrator.opts.adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*gpdt[i]
      end
      for i in eachindex(u)
        EEsttmp[i] = (integrator.opts.delta*E₁[i]+E₂[i])/(integrator.opts.abstol + max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol)
      end
      EEst = integrator.opts.internalnorm(EEsttmp)
    else
      for i in eachindex(u)
        u[i] = u[i] + k₁[i]/3 + 2k₂[i]/3 + E₂[i] + ΔW[i]*gpdt[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRA,uType<:AbstractArray,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = integrator.alg.tableau
  stages::Int = length(α)
  H0 = Vector{typeof(u)}(0)
  for i = 1:stages
    push!(H0,zeros(u))
  end
  A0temp::uType = zeros(u); B0temp::uType = zeros(u)
  ftmp::uType = zeros(u); gtmp::uType = zeros(u); chi2::uType = zeros(u)
  atemp::uType = zeros(u); btemp::uType = zeros(u); E₂::uType = zeros(u); E₁temp::uType = zeros(u)
  E₁::uType = zeros(u)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader
    for i in eachindex(u)
      chi2[i] = .5*(ΔW[i] + ΔZ[i]/sqrt(3)) #I_(1,0)/h
    end
    for i in 1:stages
      H0[i][:]=zero(uEltype)
    end
    for i = 1:stages
      A0temp[:] = zero(uEltype)
      B0temp[:] = zero(uEltype)
      for j = 1:i-1
        integrator.f(t + c₀[j]*dt,H0[j],ftmp)
        integrator.g(t + c₁[j]*dt,H0[j],gtmp)
        for k in eachindex(u)
          A0temp[k] += A₀[i,j]*ftmp[k]
          B0temp[k] += B₀[i,j]*gtmp[k]
        end
      end
      for j in eachindex(u)
        H0[i][j] = u[j] + A0temp[j]*dt + B0temp[j]*chi2[j]
      end
    end
    atemp[:] = zero(uEltype)
    btemp[:] = zero(uEltype)
    E₂[:]    = zero(uEltype)
    E₁temp[:]= zero(uEltype)

    for i = 1:stages
      integrator.f(t+c₀[i]*dt,H0[i],ftmp)
      integrator.g(t+c₁[i]*dt,H0[i],gtmp)
      for j in eachindex(u)
        atemp[j] += α[i]*ftmp[j]
        btemp[j] += (β₁[i]*ΔW[j])*gtmp[j]
        E₂[j]    += (β₂[i]*chi2[j])*gtmp[j]
        E₁temp[j] += ftmp[j]
      end
    end
    for i in eachindex(u)
      E₁[i] = dt*E₁temp[i]
    end

    if integrator.opts.adaptive
      for i in eachindex(u)
        utmp[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
      for i in eachindex(u)
        EEsttmp[i] = (integrator.opts.delta*E₁[i]+E₂[i])/(integrator.opts.abstol + max(abs(u[i]),abs(utmp[i]))*integrator.opts.reltol)
      end
      EEst = integrator.opts.internalnorm(EEsttmp)
    else
      for i in eachindex(u)
        u[i] = u[i] + dt*atemp[i] + btemp[i] + E₂[i]
      end
    end
    @sde_loopfooter
  end
  @sde_postamble
end

function sde_solve{algType<:SRA,uType<:Number,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType}(integrator::SDEIntegrator{algType,uType,uEltype,Nm1,N,tType,tTypeNoUnits,uEltypeNoUnits,randType,rateType,solType,F4,F5,OType})
  @sde_preamble
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = integrator.alg.tableau
  stages::Int = length(α)
  H0 = Array{uEltype}(stages)
  @sde_adaptiveprelim
  @inbounds while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(stages)
    for i = 1:stages
      A0temp = zero(u)
      B0temp = zero(u)
      for j = 1:i-1
        A0temp += A₀[i,j]*integrator.f(t + c₀[j]*dt,H0[j])
        B0temp += B₀[i,j]*integrator.g(t + c₁[j]*dt,H0[j]) #H0[..,i] argument ignored
      end
      H0[i] = u + A0temp*dt + B0temp.*chi2
    end

    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)

    for i = 1:stages
      ftemp = integrator.f(t+c₀[i]*dt,H0[i])
      E₁temp += ftemp
      atemp += α[i]*ftemp
      btemp += (β₁[i]*ΔW ).*integrator.g(t+c₁[i]*dt,H0[i]) #H0[i] argument ignored
      E₂    += (β₂[i]*chi2).*integrator.g(t+c₁[i]*dt,H0[i]) #H0[i] argument ignored
    end

    if integrator.opts.adaptive
      E₁ = dt*E₁temp
      utmp = u + dt*atemp + btemp + E₂
      EEst = abs((integrator.opts.delta*E₁+E₂)./(integrator.opts.abstol + max(abs(u),abs(utmp))*integrator.opts.reltol))
    else
      u = u + dt*atemp + btemp + E₂
    end
    @sde_loopfooter
  end
  @sde_postamble
end
