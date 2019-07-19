@muladd function perform_step!(integrator, cache::SKenCarpConstantCache, f=integrator.f)
  @unpack t,dt,uprev,u,g,p = integrator
  @unpack uf, nlsolver = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab
  @unpack ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4,ebtilde1,ebtilde2,ebtilde3,ebtilde4 = cache.tab
  @unpack nb021,nb043 = cache.tab
  alg = unwrap_alg(integrator, true)

  sqrt3 = sqrt(3one(eltype(integrator.W.dW)))
  chi2 = (integrator.W.dW + integrator.W.dZ/sqrt3)/2 #I_(1,0)/h

  if typeof(integrator.f) <: SplitSDEFunction
    f = integrator.f.f1
    f2 = integrator.f.f2
  else
    f = integrator.f
  end

  # precalculations
  γdt = γ*dt

  # calculate W
  repeat_step = false
  if isnewton(nlsolver)
    uf.t = t
  end
  J = update_W!(integrator, cache, γdt, repeat_step)

  z₁ = dt*f( uprev,p,t)
  nlsolver.c = γ

  g1 = g(uprev,p,t)
  tmp = uprev + γ*z₁ + nb021*chi2.*g1
  if typeof(integrator.f) <: SplitSDEFunction
    # This assumes the implicit part is cheaper than the explicit part
    k1 = dt*f2(uprev,p,t)
    tmp += ea21*k1
  end
  nlsolver.tmp = tmp

  if alg.extrapolant == :min_correct
    z₂ = zero(z₁)
  elseif alg.extrapolant == :trivial
    z₂ = z₁
  end
  nlsolver.z = z₂

  z₂ = nlsolve!(integrator, cache)
  nlsolvefail(nlsolver) && return nothing


  ################################## Solve Step 3

  nlsolver.c = c3
  if typeof(integrator.f) <: SplitSDEFunction
    u = tmp + γ*z₂
    k2 = dt*f2(u,p,t + 2γ*dt)
    tmp = uprev + a31*z₁ + a32*z₂ + ea31*k1 + ea32*k2
  else
    tmp = uprev + a31*z₁ + a32*z₂
  end
  nlsolver.tmp = tmp

  if alg.extrapolant == :min_correct
    z₃ = zero(z₂)
  elseif alg.extrapolant == :trivial
    z₃ = z₂
  end
  nlsolver.z = z₃

  z₃ = nlsolve!(integrator, cache)
  nlsolvefail(nlsolver) && return nothing

  ################################## Solve Step 4

  nlsolver.c = one(nlsolver.c)

  # Note: Can use g1 since c13 = 0 => g3 == g1

  if typeof(integrator.f) <: SplitSDEFunction
    u = tmp + γ*z₃
    k3 = dt*f2(u,p,t + c3*dt)
    tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + ea41*k1 + ea42*k2 + ea43*k3 + nb043*chi2.*g1
  else
    tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + nb043*chi2.*g1
  end
  nlsolver.tmp = tmp

  if alg.extrapolant == :min_correct
    z₄ = zero(z₂)
  elseif alg.extrapolant == :trivial
    z₄ = z₂
  end
  nlsolver.z = z₄

  z₄ = nlsolve!(integrator, cache)
  nlsolvefail(nlsolver) && return nothing

  u = tmp + γ*z₄
  g4 = g(uprev,p,t+dt)

  E₂ = chi2.*(g1-g4)

  if typeof(integrator.f) <: SplitSDEFunction
    k4 = dt*f2(u,p,t+dt)
    u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + eb1*k1 + eb2*k2 + eb3*k3 + eb4*k4 + integrator.W.dW.*g4 + E₂
  else
    u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + integrator.W.dW.*g4 + E₂
  end

  ################################### Finalize

  if integrator.opts.adaptive

    #=
    if typeof(integrator.f) <: SplitSDEFunction
      tmp = btilde1*z₁  + btilde2*z₂  + btilde3*z₃ + btilde4*z₄ + ebtilde1*k1 + ebtilde2*k2 + ebtilde3*k3 + ebtilde4*k4 + chi2*(g1-g4)
    else
      tmp = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + chi2*(g1-g4)
    end
    if alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end

    =#

    E₁ = z₁ + z₂ + z₃ + z₄

    resids = calculate_residuals(E₁, E₂, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(resids,t)
  end

  integrator.u = u
end

@muladd function perform_step!(integrator, cache::SKenCarpCache, f=integrator.f)
  repeat_step=false
  @unpack t,dt,uprev,u,g,p = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,k1,k2,k3,k4,k,b,J,W,jac_config,tmp,atmp = cache
  @unpack g1,g4,chi2,nlsolver = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab
  @unpack ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4 = cache.tab
  @unpack ebtilde1,ebtilde2,ebtilde3,ebtilde4 = cache.tab
  @unpack nb021,nb043 = cache.tab
  alg = unwrap_alg(integrator, true)

  # Some aliases

  E₁ = g4
  E₂ = dz

  if typeof(integrator.f) <: SplitSDEFunction
    f = integrator.f.f1
    f2 = integrator.f.f2
  else
    f = integrator.f
  end

  sqrt3 = sqrt(3one(eltype(integrator.W.dW)))
  if typeof(integrator.W.dW) <: Union{SArray,Number}
    chi2 = (integrator.W.dW + integrator.W.dZ/sqrt3)/2 #I_(1,0)/h
  else
    @.. chi2 = (integrator.W.dW + integrator.W.dZ/sqrt3)/2 #I_(1,0)/h
  end

  # precalculations

  γdt = γ*dt

  update_W!(integrator, cache, γdt, repeat_step)

  if !repeat_step && !integrator.last_stepfail
    f(z₁, integrator.uprev, p, integrator.t)
    z₁ .*= dt
  end

  ##### Step 2

  # TODO: Add a cache so this isn't overwritten near the end, so it can not repeat on fail
  g(g1,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. z₄ = chi2*g1 # use z₄ as storage for the g1*chi2
  else
    mul!(z₄,g1,chi2) # use z₄ as storage for the g1*chi2
  end

  @.. tmp = uprev + γ*z₁ + nb021*z₄

  if alg.extrapolant == :min_correct
    @.. z₂ = zero(eltype(dz))
  elseif alg.extrapolant == :trivial
    @.. z₂ = z₁
  end
  nlsolver.z = z₂
  nlsolver.c = γ


  if typeof(integrator.f) <: SplitSDEFunction
    # This assumes the implicit part is cheaper than the explicit part
    if !repeat_step && !integrator.last_stepfail
      f2(k1,integrator.uprev,integrator.p,integrator.t)
      k1 .*= dt
    end
    @.. tmp += ea21*k1
  end

  z₂ = nlsolve!(integrator, cache)
  nlsolvefail(nlsolver) && return nothing

  ################################## Solve Step 3

  nlsolver.c = c3

  if typeof(integrator.f) <: SplitSDEFunction
    @.. u = tmp + γ*z₂
    f2(k2,u,p,t + 2γ*dt); k2 .*= dt
    @.. tmp = uprev + a31*z₁ + a32*z₂ + ea31*k1 + ea32*k2
  else
    @.. tmp = uprev + a31*z₁ + a32*z₂
  end

  if alg.extrapolant == :min_correct
    @.. z₃ = zero(eltype(dz))
  elseif alg.extrapolant == :trivial
    @.. z₃ = z₂
  end
  nlsolver.z = z₃

  z₃ = nlsolve!(integrator, cache)
  nlsolvefail(nlsolver) && return nothing

  ################################## Solve Step 4

  nlsolver.c = one(nlsolver.c)

  if typeof(integrator.f) <: SplitSDEFunction
    @.. u = tmp + γ*z₃
    f2(k3,u,p,t + c3*dt); k3 .*= dt
    # z₄ is storage for the g1*chi2 from earlier
    @.. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + ea41*k1 + ea42*k2 + ea43*k3 + nb043*z₄
  else
    @unpack α41,α42 = cache.tab
    # z₄ is storage for the g1*chi2
    @.. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + nb043*z₄
  end

  if alg.extrapolant == :min_correct
    @.. z₄ = zero(eltype(dz))
  elseif alg.extrapolant == :trivial
    @.. z₄ = z₂
  end
  nlsolver.z = z₄

  z₄ = nlsolve!(integrator, cache)
  nlsolvefail(nlsolver) && return nothing

  g(g4,u,p,t+dt)

  if typeof(integrator.f) <: SplitSDEFunction
    @.. u = tmp + γ*z₄
    f2(k4,u,p,t+dt); k4 .*= dt
    if is_diagonal_noise(integrator.sol.prob)
      @.. E₂ = chi2*(g1-g4)
      @.. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + eb1*k1 + eb2*k2 + eb3*k3 + eb4*k4 + integrator.W.dW*g4 + E₂
    else
      g1 .-= g4
      mul!(E₂,g1,chi2)
      mul!(tmp,g4,integrator.W.dW)
      @.. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + eb1*k1 + eb2*k2 + eb3*k3 + eb4*k4 + tmp + E₂
    end
  else
    if is_diagonal_noise(integrator.sol.prob)
      @.. E₂ = chi2*(g1-g4)
      @.. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + integrator.W.dW*g4 + E₂
    else
      g1 .-= g4
      mul!(E₂,g1,chi2)
      mul!(tmp,g4,integrator.W.dW)
      @.. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + tmp + E₂
    end
  end

  ################################### Finalize

  if integrator.opts.adaptive

    #=
    if typeof(integrator.f) <: SplitSDEFunction
      if is_diagonal_noise(integrator.sol.prob)
        #@.. dz = btilde1*z₁  + btilde2*z₂  + btilde3*z₃ + btilde4*z₄ + ebtilde1*k1 + ebtilde2*k2 + ebtilde3*k3 + ebtilde4*k4
        for i in eachindex(dz)
          @inbounds dz[i] = btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + ebtilde1*k1[i] + ebtilde2*k2[i] + ebtilde3*k3[i] + ebtilde4*k4[i] + chi2[i]*(g1[i]-g4[i])
        end
      else
        # dz already holds mul!(dz,g1,chi2)!
        #@.. dz += btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + ebtilde1*k1 + ebtilde2*k2 + ebtilde3*k3 + ebtilde4*k4
        for i in eachindex(dz)
          @inbounds dz[i] += btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + ebtilde1*k1[i] + ebtilde2*k2[i] + ebtilde3*k3[i] + ebtilde4*k4[i]
        end
      end
    else
      if is_diagonal_noise(integrator.sol.prob)
        #@.. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + chi2*(g1-g4)
        for i in eachindex(dz)
          @inbounds dz[i] = btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + chi2[i]*(g1[i]-g4[i])
        end
      else
        # dz already holds mul!(dz,g1,chi2)!
        @.. dz += btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
      end
    end
    if alg.smooth_est # From Shampine
      cache.linsolve(vec(tmp),W,vec(dz),false)
    else
      tmp .= dz
    end
    =#

    @.. E₁ = z₁ + z₂ + z₃ + z₄

    calculate_residuals!(tmp, E₁, E₂, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)
  end
end
