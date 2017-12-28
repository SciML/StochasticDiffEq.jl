@muladd function perform_step!(integrator, cache::RackKenCarpConstantCache)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab

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

  # FSAL Step 1
  z₁ = dt.*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation for guess
  z₂ = z₁

  iter = 1
  tstep = t + 2*γdt
  tmp = @. uprev + γ*z₁
  u = @. tmp + γ*z₂
  b = dt.*f(tstep,u) .- z₂
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₂ = z₂ .+ dz

  η = max(cache.ηold,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = integrator.success_iter == 0 || η*ndz > κtol

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₂
    b = dt.*f(tstep,u) .- z₂
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
    z₂ = z₂ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 3

  # Guess is from Hermite derivative on z₁ and z₂
  z₃ = @. α31*z₁ + α32*z₂

  iter = 1
  tstep = t + c3*dt
  tmp = @. uprev + a31*z₁ + a32*z₂
  u = @. tmp + γ*z₃
  b = dt.*f(tstep,u) .- z₃
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₃ = z₃ .+ dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₃
    b = dt.*f(tstep,u) .- z₃
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
    z₃ = z₃ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  ################################## Solve Step 4

  @unpack α41,α42 = cache.tab
  z₄ = @. α41*z₁ + α42*z₂

  iter = 1
  tstep = t + dt
  tmp = @. uprev + a41*z₁ + a42*z₂ + a43*z₃
  u = @. tmp + γ*z₄
  b = dt.*f(tstep,u) .- z₄
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z₄ = z₄ .+ dz

  η = max(η,eps(eltype(integrator.opts.reltol)))^(0.8)
  do_newton = (η*ndz > κtol)

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    u = @. tmp + γ*z₄
    b = dt.*f(tstep,u) .- z₄
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
    z₄ = z₄ .+ dz
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  u = @. tmp + γ*z₄

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    tmp = @. btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
    if integrator.alg.smooth_est # From Shampine
      est = W\tmp
    else
      est = tmp
    end
    atmp = calculate_residuals(est, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  integrator.fsallast = z₄./dt
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

@muladd function perform_step!(integrator, cache::RackKenCarpCache)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uf,du1,dz,z₁,z₂,z₃,z₄,k,b,J,W,jac_config,tmp,atmp,κ,tol = cache
  @unpack γ,a31,a32,a41,a42,a43,btilde1,btilde2,btilde3,btilde4,c3,α31,α32 = cache.tab

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

  # FSAL Step 1
  @. z₁ = dt*integrator.fsalfirst

  ##### Step 2

  # TODO: Add extrapolation for guess
  @. z₂ = z₁

  # initial step of Newton iteration
  iter = 1
  tstep = t + 2*γdt
  @. tmp = uprev + γ*z₁
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

  # Guess is from Hermite derivative on z₁ and z₂
  @. z₃ = α31*z₁ + α32*z₂

  # initial step of Newton iteration
  iter = 1
  tstep = t + c3*dt
  @. tmp = uprev + a31*z₁ + a32*z₂
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
    u = @. tmp + γ*z₃
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

  @unpack α41,α42 = cache.tab
  @. z₄ = α41*z₁ + α42*z₂

  iter = 1
  tstep = t + dt
  @. tmp = uprev + a41*z₁ + a42*z₂ + a43*z₃
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

  @. u = tmp + γ*z₄

  ################################### Finalize

  cache.ηold = η
  cache.newton_iters = iter

  if integrator.opts.adaptive
    @. dz = btilde1*z₁ + btilde2*z₂ + btilde3*z₃ + btilde4*z₄
    if integrator.alg.smooth_est # From Shampine
      if has_invW(f)
        A_mul_B!(vec(tmp),W,vec(dz))
      else
        cache.linsolve(vec(tmp),W,vec(dz),false)
      end
    else
      tmp .= dz
    end
    calculate_residuals!(atmp, tmp, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end

  @. integrator.fsallast = z₄/dt
end
