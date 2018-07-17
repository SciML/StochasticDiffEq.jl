@muladd function perform_step!(integrator, cache::SKenCarpConstantCache, f=integrator.f)

  repeat_step=false

  @unpack t,dt,uprev,u,g,p = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab
  @unpack ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4,ebtilde1,ebtilde2,ebtilde3,ebtilde4 = cache.tab
  @unpack nb021,nb043 = cache.tab
  alg = unwrap_alg(integrator, true)

  chi2 = (integrator.W.dW + integrator.W.dZ/sqrt(3))/2 #I_(1,0)/h

  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
    f2 = integrator.f.f2
  else
    f = integrator.f
  end

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  # calculate W
  uf.t = t
  repeat_step = false
  J, W = calc_W!(integrator, cache, γdt, repeat_step)

  z₁ = dt*f( uprev,p,t)

  ##### Step 2

  iter = 1
  tstep = t + 2*γdt

  g1 = g(uprev,p,t)

  tmp = uprev + γ*z₁ + chi2*nb021*g1

  if typeof(integrator.f) <: SplitFunction
    # This assumes the implicit part is cheaper than the explicit part
    k1 = dt*f2(uprev,p,t)
    tmp += ea21*k1
  end

  if alg.extrapolant == :min_correct
    z₂ = zero(z₁)
  elseif alg.extrapolant == :trivial
    z₂ = z₁
  end

  u = tmp + γ*z₂
  b = dt*f(u,p,tstep) - z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₂
    b = dt*f(u,p,tstep) - z₂
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ = z₂ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  iter = 1
  tstep = t + c3*dt

  if typeof(integrator.f) <: SplitFunction
    u = tmp + γ*z₂
    k2 = dt*f2(u,p,t + 2γ*dt)
    tmp = uprev + a31*z₁ + a32*z₂ + ea31*k1 + ea32*k2
  else
    tmp = uprev + a31*z₁ + a32*z₂
  end

  if alg.extrapolant == :min_correct
    z₃ = zero(dz)
  elseif alg.extrapolant == :trivial
    z₃ = z₂
  end

  u = tmp + γ*z₃
  b = dt*f(u,p,tstep) - z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₃
    b = dt*f(u,p,tstep) - z₃
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ = z₃ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  iter = 1
  tstep = t + dt

  # Note: Can use g1 since c13 = 0 => g3 == g1

  if typeof(integrator.f) <: SplitFunction
    u = tmp + γ*z₃
    k3 = dt*f2(u,p,t + c3*dt)
    tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + ea41*k1 + ea42*k2 + ea43*k3 + chi2*nb043*g1
  else
    tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + chi2*nb043*g1
  end

  if alg.extrapolant == :min_correct
    z₄ = zero(dz)
  elseif alg.extrapolant == :trivial
    z₄ = z₂
  end

  u = tmp + γ*z₄
  b = dt*f(u,p,tstep) - z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₄
    b = dt*f(u,p,tstep) - z₄
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ = z₄ + dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + γ*z₄
  g4 = g(uprev,p,t+dt)

  E₂ = chi2*(g1-g4)

  if typeof(integrator.f) <: SplitFunction
    k4 = dt*f2(u,p,t+dt)
    u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + eb1*k1 + eb2*k2 + eb3*k3 + eb4*k4 + integrator.W.dW*g4 + E₂
  else
    u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + integrator.W.dW*g4 + E₂
  end

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive

    #=
    if typeof(integrator.f) <: SplitFunction
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

    integrator.EEst = integrator.opts.internalnorm((integrator.opts.delta*E₁+E₂)./(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  end

  integrator.u = u
end

@muladd function perform_step!(integrator, cache::SKenCarpCache, f=integrator.f)
  repeat_step=false
  @unpack t,dt,uprev,u,g,p = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,k1,k2,k3,k4,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack g1,g4,chi2 = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab
  @unpack ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4 = cache.tab
  @unpack ebtilde1,ebtilde2,ebtilde3,ebtilde4 = cache.tab
  @unpack nb021,nb043 = cache.tab
  alg = unwrap_alg(integrator, true)

  # Some aliases

  E₁ = g4
  E₂ = dz

  if typeof(integrator.f) <: SplitFunction
    f = integrator.f.f1
    f2 = integrator.f.f2
  else
    f = integrator.f
  end

  if typeof(integrator.W.dW) <: Union{SArray,Number}
    chi2 = (integrator.W.dW + integrator.W.dZ/sqrt(3))/2 #I_(1,0)/h
  else
    @. chi2 = (integrator.W.dW + integrator.W.dZ/sqrt(3))/2 #I_(1,0)/h
  end

  # precalculations
  κtol = κ*tol

  γdt = γ*dt

  new_W = calc_W!(integrator, cache, γdt, repeat_step)

  if !repeat_step && !integrator.last_stepfail
    f(z₁, integrator.uprev, p, integrator.t)
    z₁ .*= dt
  end

  ##### Step 2



  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt

  # TODO: Add a cache so this isn't overwritten near the end, so it can not repeat on fail
  g(g1,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @. z₄ = chi2*g1 # use z₄ as storage for the g1*chi2
  else
    mul!(z₄,g1,chi2) # use z₄ as storage for the g1*chi2
  end

  @. tmp = uprev + γ*z₁ + nb021*z₄

  if alg.extrapolant == :min_correct
    @. z₂ = zero(dz)
  elseif alg.extrapolant == :trivial
    @. z₂ = z₁
  end



  if typeof(integrator.f) <: SplitFunction
    # This assumes the implicit part is cheaper than the explicit part
    if !repeat_step && !integrator.last_stepfail
      f2(k1,integrator.uprev,integrator.p,integrator.t)
      k1 .*= dt
    end
    @. tmp += ea21*k1
  end

  @. u = tmp + γ*z₂
  f(k,u,p,tstep)
  @. b = dt*k - z₂
  if has_invW(f)
    mul!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(k,u,p,tstep)
    @. b = dt*k - z₂
    if has_invW(f)
      mul!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt

  if typeof(integrator.f) <: SplitFunction
    @. u = tmp + γ*z₂
    f2(k2,u,p,t + 2γ*dt); k2 .*= dt
    #@. tmp = uprev + a31*z₁ + a32*z₂ + ea31*k1 + ea32*k2
    for i in eachindex(tmp)
      @inbounds tmp[i] = uprev[i] + a31*z₁[i] + a32*z₂[i] + ea31*k1[i] + ea32*k2[i]
    end
  else
    @. tmp = uprev + a31*z₁ + a32*z₂
  end

  if alg.extrapolant == :min_correct
    @. z₃ = zero(dz)
  elseif alg.extrapolant == :trivial
    @. z₃ = z₂
  end

  @. u = tmp + γ*z₃
  f(k,u,p,tstep)
  @. b = dt*k - z₃
  if has_invW(f)
    mul!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₃
    f(k,u,p,tstep)
    @. b = dt*k - z₃
    if has_invW(f)
      mul!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  iter = 1
  tstep = t + dt

  if typeof(integrator.f) <: SplitFunction
    @. u = tmp + γ*z₃
    f2(k3,u,p,t + c3*dt); k3 .*= dt
    # z₄ is storage for the g1*chi2 from earlier
    #@. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + ea41*k1 + ea42*k2 + ea43*k3
    for i in eachindex(tmp)
      @inbounds tmp[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + ea41*k1[i] + ea42*k2[i] + ea43*k3[i] + nb043*z₄[i]
    end
  else
    @unpack α41,α42 = cache.tab
    # z₄ is storage for the g1*chi2
    @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + nb043*z₄
  end

  if alg.extrapolant == :min_correct
    @. z₄ = zero(dz)
  elseif alg.extrapolant == :trivial
    @. z₄ = z₂
  end

  @. u = tmp + γ*z₄
  f(k,u,p,tstep)
  @. b = dt*k - z₄
  if has_invW(f)
    mul!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < alg.min_newton_iter) && iter < alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(k,u,p,tstep)
    @. b = dt*k - z₄
    if has_invW(f)
      mul!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ .+= dz
  end

  if (iter >= alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end


  g(g4,u,p,t+dt)

  if typeof(integrator.f) <: SplitFunction
    @. u = tmp + γ*z₄
    f2(k4,u,p,t+dt); k4 .*= dt
    if is_diagonal_noise(integrator.sol.prob)
      @. E₂ = chi2*(g1-g4)
      #@. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + eb1*k1 + eb2*k2 + eb3*k3 + eb4*k4 + integrator.W.dW*g4 + chi2*(g1-g4)
      for i in eachindex(u)
        @inbounds u[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + γ*z₄[i] + eb1*k1[i] + eb2*k2[i] + eb3*k3[i] + eb4*k4[i] + integrator.W.dW[i]*g4[i] + E₂[i]
      end
    else
      g1 .-= g4
      mul!(E₂,g1,chi2)
      mul!(tmp,g4,integrator.W.dW)
      for i in eachindex(u)
        @inbounds u[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + γ*z₄[i] + eb1*k1[i] + eb2*k2[i] + eb3*k3[i] + eb4*k4[i] + tmp[i] + E₂[i]
      end
    end
  else
    if is_diagonal_noise(integrator.sol.prob)
      @. E₂ = chi2*(g1-g4)
      #@. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + integrator.W.dW*g4 + chi2*(g1-g4)
      for i in eachindex(u)
        @inbounds u[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + γ*z₄[i] + integrator.W.dW[i]*g4[i] + E₂[i]
      end
    else
      g1 .-= g4
      mul!(E₂,g1,chi2)
      mul!(tmp,g4,integrator.W.dW)
      for i in eachindex(u)
        @inbounds u[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + γ*z₄[i] + tmp[i] + E₂[i]
      end
    end
  end

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive

    #=
    if typeof(integrator.f) <: SplitFunction
      if is_diagonal_noise(integrator.sol.prob)
        #@. dz = btilde1*z₁  + btilde2*z₂  + btilde3*z₃ + btilde4*z₄ + ebtilde1*k1 + ebtilde2*k2 + ebtilde3*k3 + ebtilde4*k4
        for i in eachindex(dz)
          @inbounds dz[i] = btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + ebtilde1*k1[i] + ebtilde2*k2[i] + ebtilde3*k3[i] + ebtilde4*k4[i] + chi2[i]*(g1[i]-g4[i])
        end
      else
        # dz already holds mul!(dz,g1,chi2)!
        #@. dz += btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + ebtilde1*k1 + ebtilde2*k2 + ebtilde3*k3 + ebtilde4*k4
        for i in eachindex(dz)
          @inbounds dz[i] += btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + ebtilde1*k1[i] + ebtilde2*k2[i] + ebtilde3*k3[i] + ebtilde4*k4[i]
        end
      end
    else
      if is_diagonal_noise(integrator.sol.prob)
        #@. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄ + chi2*(g1-g4)
        for i in eachindex(dz)
          @inbounds dz[i] = btilde1*z₁[i] + btilde2*z₂[i] + btilde3*z₃[i] + btilde4*z₄[i] + chi2[i]*(g1[i]-g4[i])
        end
      else
        # dz already holds mul!(dz,g1,chi2)!
        @. dz += btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
      end
    end

    if alg.smooth_est # From Shampine
      if has_invW(f)
        mul!(vec(tmp),W,vec(dz))
      else
        cache.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end

    =#

    @. E₁ = z₁ + z₂ + z₃ + z₄

    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
						  Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = (δ*E₁[i]+E₂[i])/(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

end
