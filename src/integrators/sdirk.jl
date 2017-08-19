@muladd function perform_step!(integrator,
                               cache::Union{ImplicitEMConstantCache,
                                            ImplicitEulerHeunConstantCache,
                                            ImplicitRKMilConstantCache},
                                            f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf = cache
  theta = integrator.alg.theta
  integrator.alg.symplectic ? a = dt/2 : a = dt
  uf.t = t

  # TODO: Stochastic extrapolants?
  u = uprev

  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    W = I - dt*theta*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - dt*theta*J
  end

  z = u - uprev
  iter = 0
  κ = cache.κ
  tol = cache.tol

  L = integrator.g(t,uprev)
  gtmp = L.*integrator.W.dW

  if typeof(cache) <: ImplicitRKMilConstantCache || integrator.alg.theta != 1
    ftmp = integrator.f(t,uprev)
  else
    ftmp = 0
  end

  if typeof(cache) <: ImplicitEulerHeunConstantCache
    utilde = @. uprev + gtmp
    gtmp = @. ((integrator.g(t,utilde) + L)/2)*integrator.W.dW
  end

  if typeof(cache) <: ImplicitRKMilConstantCache
    if alg_interpretation(integrator.alg) == :Ito
      K = @muladd uprev .+ dt.*ftmp
      utilde = @.  K + L*integrator.sqdt
      mil_correction = (integrator.g(t,utilde).-L)./(2 .* integrator.sqdt).*
                       (integrator.W.dW.^2 .- dt)
      gtmp += mil_correction
    elseif alg_interpretation(integrator.alg) == :Stratonovich
      utilde = @. uprev + L*integrator.sqdt
      mil_correction = (integrator.g(t,utilde).-L)./(2 .* integrator.sqdt).*
                       (integrator.W.dW.^2)
      gtmp += mil_correction
    end
  end

  iter += 1
  if integrator.alg.symplectic
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    u = @. uprev + z/2 + gtmp/2
  else
    u = @. uprev + dt*(1-theta)*ftmp + theta*z + gtmp
  end
  b = -z .+ dt.*f(t+a,u)
  dz = W\b
  ndz = integrator.opts.internalnorm(dz)
  z = z + dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    if integrator.alg.symplectic
      # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
      u = @. uprev + z/2 + gtmp/2
    else
      u = @. uprev + dt*(1-theta)*ftmp + theta*z + gtmp
    end
    b = -z .+ dt.*f(t+a,u)
    dz = W\b
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z = z + dz
  end

  if integrator.alg.symplectic
    u = @. uprev + z + gtmp
  else
    u = @. uprev + dt*(1-theta)*ftmp + theta*z + gtmp
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = @. uprev + dt*(1-theta)*ftmp + theta*z + gtmp

  #=
  if integrator.opts.adaptive && integrator.success_iter > 0
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    DD3 = ((u - uprev)/((dt)*(t+dt-tprev)) + (uprev-uprev2)/((t-tprev)*(t+dt-tprev)))
    dEst = (dt^2)*abs(DD3/6)
    integrator.EEst = @. dEst/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
  else
    integrator.EEst = 1
  end
  =#
  integrator.u = u
end

@muladd function perform_step!(integrator,
                               cache::Union{ImplicitEMCache,
                                            ImplicitEulerHeunCache,
                                            ImplicitRKMilCache},
                               f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config,gtmp,gtmp2,tmp = cache
  integrator.alg.symplectic ? a = dt/2 : a = dt
  dW = integrator.W.dW
  mass_matrix = integrator.sol.prob.mass_matrix
  theta = integrator.alg.theta

  if integrator.success_iter > 0 && !integrator.u_modified && integrator.alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif integrator.alg.extrapolant == :linear
    @. u = uprev + integrator.fsalfirst*dt
  else # :constant
    copy!(u,uprev)
  end

  uf.t = t

  if has_invW(f)
    f(Val{:invW},t,uprev,dt,W) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},t,uprev,J)
      else
        if alg_autodiff(integrator.alg)
          ForwardDiff.jacobian!(J,uf,vec(du1),vec(uprev),jac_config)
        else
          Calculus.finite_difference_jacobian!(uf,vec(uprev),vec(du1),J,integrator.alg.diff_type)
        end
      end
    end
    if integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps()
      new_W = true
      thetadt = dt*theta
      for j in 1:length(u), i in 1:length(u)
          @inbounds W[i,j] = mass_matrix[i,j]-thetadt*J[i,j]
      end
    else
      new_W = false
    end
  end

  ##############################################################################

  # Handle noise computations

  integrator.g(t,uprev,gtmp)

  if typeof(cache) <: ImplicitRKMilConstantCache || integrator.alg.theta != 1
    f(t,uprev,tmp)
  else
    tmp .= 0
  end

  if is_diagonal_noise(integrator.sol.prob)
    @tight_loop_macros for i in eachindex(u)
      @inbounds gtmp2[i]=gtmp[i]*dW[i]
    end
  else
    A_mul_B!(gtmp2,gtmp,dW)
  end

  if typeof(cache) <: ImplicitEulerHeunCache
    @. z = uprev + gtmp2
    integrator.g(t,z,gtmp2)
    @. gtmp = (gtmp2 + gtmp)/2
    if is_diagonal_noise(integrator.sol.prob)
      @tight_loop_macros for i in eachindex(u)
        @inbounds gtmp2[i]=gtmp[i]*dW[i]
      end
    else
      A_mul_B!(gtmp2,gtmp,dW)
    end
  end

  if typeof(cache) <: ImplicitRKMilCache
    gtmp3 = cache.gtmp3
    if alg_interpretation(integrator.alg) == :Ito
      @. z = @muladd uprev + dt*du1 + gtmp*integrator.sqdt
      integrator.g(t,z,gtmp3)
      @. gtmp2 += (gtmp3-gtmp)/(2integrator.sqdt)*(dW.^2 - dt)
    elseif alg_interpretation(integrator.alg) == :Stratonovich
      @. z = @muladd uprev + gtmp*integrator.sqdt
      integrator.g(t,z,gtmp3)
      @. gtmp2 += (gtmp3-gtmp)/(2integrator.sqdt)*(dW.^2)
    end
  end

  ##############################################################################

  @. z = u - uprev
  iter = 0
  κ = cache.κ
  tol = cache.tol
  if integrator.alg.symplectic
    @. u = uprev + z/2 + gtmp2/2
  else
    @. u = uprev + dt*(1-theta)*tmp + theta*z + gtmp2
  end
  iter += 1
  f(t+a,u,k)
  scale!(k,dt)
  if mass_matrix == I
    k .-= z
  else
    A_mul_B!(du1,mass_matrix,z)
    k .-= du1
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    integrator.alg.linsolve(vec(dz),W,vec(k),new_W)
  end
  ndz = integrator.opts.internalnorm(dz)
  z .+= dz

  η = max(cache.ηold,eps(first(u)))^(0.8)
  if integrator.success_iter > 0
    do_newton = (η*ndz > κ*tol)
  else
    do_newton = true
  end

  fail_convergence = false
  while (do_newton || iter < integrator.alg.min_newton_iter) && iter < integrator.alg.max_newton_iter
    iter += 1
    if integrator.alg.symplectic
      @. u = uprev + z/2 + gtmp2/2
    else
      @. u = uprev + dt*(1-theta)*tmp + theta*z + gtmp2
    end
    f(t+a,u,k)
    scale!(k,dt)
    if mass_matrix == I
      k .-= z
    else
      A_mul_B!(du1,mass_matrix,z)
      k .-= du1
    end
    if has_invW(f)
      A_mul_B!(dz,W,k) # Here W is actually invW
    else
      integrator.alg.linsolve(vec(dz),W,vec(k),false)
    end
    ndzprev = ndz
    ndz = integrator.opts.internalnorm(dz)
    θ = ndz/ndzprev
    if θ > 1 || ndz*(θ^(integrator.alg.max_newton_iter - iter)/(1-θ)) > κ*tol
      fail_convergence = true
      break
    end
    η = θ/(1-θ)
    do_newton = (η*ndz > κ*tol)
    z .+= dz
  end

  if integrator.alg.symplectic
    @. u = uprev + z + gtmp2
  else
    @. u = uprev + dt*(1-theta)*tmp + theta*z + gtmp2
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  #=
  if integrator.opts.adaptive && integrator.success_iter > 0
    # Use 2rd divided differences a la SPICE and Shampine
    uprev2 = integrator.uprev2
    tprev = integrator.tprev
    dt1 = (dt)*(t+dt-tprev)
    dt2 = (t-tprev)*(t+dt-tprev)
    @tight_loop_macros for (i,atol,rtol) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds DD3 = (u[i] - uprev[i])/dt1 + (uprev[i]-uprev2[i])/dt2
      dEst = (dt^2)*abs(DD3)/6
      @inbounds k[i] = dEst/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(k)
  else
    integrator.EEst = 1
  end
  =#
end
