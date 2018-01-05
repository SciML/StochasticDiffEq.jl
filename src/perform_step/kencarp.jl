@muladd function perform_step!(integrator, cache::RackKenCarpConstantCache)

  repeat_step=false

  @unpack t,dt,uprev,u,g = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab
  @unpack ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4,ebtilde1,ebtilde2,ebtilde3,ebtilde4 = cache.tab
  @unpack nb021,nb043 = cache.tab

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
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
  end

  z₁ = dt*f(t, uprev)

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = z₁

  iter = 1
  tstep = t + 2*γdt

  g1 = g(t,uprev)

  tmp = uprev + γ*z₁ + chi2*nb021*g1

  if typeof(integrator.f) <: SplitFunction
    # This assumes the implicit part is cheaper than the explicit part
    k1 = dt*f2(t,uprev)
    tmp += ea21*k1
  end

  u = tmp + γ*z₂
  b = dt*f(tstep,u) - z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ + dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₂
    b = dt*f(tstep,u) - z₂
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ = z₂ + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  iter = 1
  tstep = t + c3*dt

  if typeof(integrator.f) <: SplitFunction
    z₃ = z₂
    u = tmp + γ*z₂
    k2 = dt*f2(t + 2γ*dt, u)
    tmp = uprev + a31*z₁ + a32*z₂ + ea31*k1 + ea32*k2
  else
    # Guess is from Hermite derivative on z₁ and z₂
    #z₃ = α31*z₁ + α32*z₂
    z₃ = z₂
    tmp = uprev + a31*z₁ + a32*z₂
  end

  u = tmp + γ*z₃
  b = dt*f(tstep,u) - z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₃
    b = dt*f(tstep,u) - z₃
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ = z₃ + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  iter = 1
  tstep = t + dt

  # Note: Can use g1 since c13 = 0 => g3 == g1

  if typeof(integrator.f) <: SplitFunction
    z₄ = z₂
    u = tmp + γ*z₃
    k3 = dt*f2(t + c3*dt, u)
    tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + ea41*k1 + ea42*k2 + ea43*k3 + chi2*nb043*g1
  else
    @unpack α41,α42 = cache.tab
    #z₄ = α41*z₁ + α42*z₂
    z₄ = z₂
    tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + chi2*nb043*g1
  end

  u = tmp + γ*z₄
  b = dt*f(tstep,u) - z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ + dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = tmp + γ*z₄
    b = dt*f(tstep,u) - z₄
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ = z₄ + dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = tmp + γ*z₄
  g4 = g(t+dt,uprev)

  E₂ = chi2*(g1-g4)

  if typeof(integrator.f) <: SplitFunction
    k4 = dt*f2(t+dt, u)
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
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end

    =#

    E₁ = z₁ + z₂ + z₃ + z₄

    integrator.EEst = integrator.opts.internalnorm(@muladd(integrator.opts.delta*E₁+E₂)./@muladd(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  end

  integrator.u = u
end

@muladd function perform_step!(integrator, cache::RackKenCarpCache)
  repeat_step=false
  @unpack t,dt,uprev,u,g = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,k1,k2,k3,k4,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack g1,g4,chi2 = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab
  @unpack ea21,ea31,ea32,ea41,ea42,ea43,eb1,eb2,eb3,eb4 = cache.tab
  @unpack ebtilde1,ebtilde2,ebtilde3,ebtilde4 = cache.tab
  @unpack nb021,nb043 = cache.tab

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

  if has_invW(f)
    # skip calculation of inv(W) if step is repeated
    !repeat_step && f(Val{:invW},t,uprev,γdt,W) # W == inverse W
  else
    # skip calculation of J if step is repeated
    if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound)
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
      end
    end
    # skip calculation of W if step is repeated
    if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t)))
      new_W = true
      mass_matrix = integrator.sol.prob.mass_matrix
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
      end
    else
      new_W = false
    end
  end

  if !repeat_step && !integrator.last_stepfail
    f(integrator.t, integrator.uprev, z₁)
    z₁ .*= dt
  end

  ##### Step 2

  # TODO: Add extrapolation for guess
  @. z₂ = z₁
  #@. z₂ = zero(u)

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt

  # TODO: Add a cache so this isn't overwritten near the end, so it can not repeat on fail
  g(t,uprev,g1)

  if is_diagonal_noise(integrator.sol.prob)
    @. z₄ = chi2*g1 # use z₄ as storage for the g1*chi2
  else
    A_mul_B!(z₄,g1,chi2) # use z₄ as storage for the g1*chi2
  end

  @. tmp = uprev + γ*z₁ + nb021*z₄

  if typeof(integrator.f) <: SplitFunction
    # This assumes the implicit part is cheaper than the explicit part
    if !repeat_step && !integrator.last_stepfail
      f2(integrator.t, integrator.uprev, k1)
      k1 .*= dt
    end
    @. tmp += ea21*k1
  end

  @. u = tmp + γ*z₂
  f(tstep,u,k)
  @. b = dt*k - z₂
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₂ .+= dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₂
    f(tstep,u,k)
    @. b = dt*k - z₂
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₂ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt

  if typeof(integrator.f) <: SplitFunction
    z₃ .= z₂
    #@. z₃ = zero(u)
    @. u = tmp + γ*z₂
    f2(t + 2γ*dt, u, k2); k2 .*= dt
    #@. tmp = uprev + a31*z₁ + a32*z₂ + ea31*k1 + ea32*k2
    for i in eachindex(tmp)
      @inbounds tmp[i] = uprev[i] + a31*z₁[i] + a32*z₂[i] + ea31*k1[i] + ea32*k2[i]
    end
  else
    # Guess is from Hermite derivative on z₁ and z₂
    #@. z₃ = α31*z₁ + α32*z₂
    z₃ .= z₂
    #@. z₃ = zero(u)
    @. tmp = uprev + a31*z₁ + a32*z₂
  end

  @. u = tmp + γ*z₃
  f(tstep,u,k)
  @. b = dt*k - z₃
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₃ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₃
    f(tstep,u,k)
    @. b = dt*k - z₃
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₃ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  iter = 1
  tstep = t + dt

  if typeof(integrator.f) <: SplitFunction
    @. u = tmp + γ*z₃
    f2(t + c3*dt, u, k3); k3 .*= dt
    # z₄ is storage for the g1*chi2 from earlier
    #@. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + ea41*k1 + ea42*k2 + ea43*k3
    for i in eachindex(tmp)
      @inbounds tmp[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + ea41*k1[i] + ea42*k2[i] + ea43*k3[i] + nb043*z₄[i]
    end
    #@. z₄ = zero(u)
    z₄ .= z₂
  else
    @unpack α41,α42 = cache.tab
    # z₄ is storage for the g1*chi2
    @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃ + nb043*z₄
    #@. z₄ = zero(u)
    z₄ .= z₂
    #@. z₄ = α41*z₁ + α42*z₂
  end

  @. u = tmp + γ*z₄
  f(tstep,u,k)
  @. b = dt*k - z₄
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(b),false)
  end
  ndz = integrator.opts.internalnorm(dz)
  z₄ .+= dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  # Newton iteration
  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    @. u = tmp + γ*z₄
    f(tstep,u,k)
    @. b = dt*k - z₄
    if has_invW(f)
      A_mul_B!(vec(dz),W,vec(b)) # Here W is actually invW
    else
      cache.linsolve(vec(dz),W,vec(b),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κtol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κtol)
    z₄ .+= dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end


  g(t+dt,u,g4)

  if typeof(integrator.f) <: SplitFunction
    @. u = tmp + γ*z₄
    f2(t+dt, u, k4); k4 .*= dt
    if is_diagonal_noise(integrator.sol.prob)
      @. E₂ = chi2*(g1-g4)
      #@. u = uprev + a41*z₁ + a42*z₂ + a43*z₃ + γ*z₄ + eb1*k1 + eb2*k2 + eb3*k3 + eb4*k4 + integrator.W.dW*g4 + chi2*(g1-g4)
      for i in eachindex(u)
        @inbounds u[i] = uprev[i] + a41*z₁[i] + a42*z₂[i] + a43*z₃[i] + γ*z₄[i] + eb1*k1[i] + eb2*k2[i] + eb3*k3[i] + eb4*k4[i] + integrator.W.dW[i]*g4[i] + E₂[i]
      end
    else
      g1 .-= g4
      A_mul_B!(E₂,g1,chi2)
      A_mul_B!(tmp,g4,integrator.W.dW)
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
      A_mul_B!(E₂,g1,chi2)
      A_mul_B!(tmp,g4,integrator.W.dW)
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
        # dz already holds A_mul_B!(dz,g1,chi2)!
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
        # dz already holds A_mul_B!(dz,g1,chi2)!
        @. dz += btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
      end
    end

    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
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
      @inbounds tmp[i] = @muladd(δ*E₁[i]+E₂[i])/@muladd(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end

end
