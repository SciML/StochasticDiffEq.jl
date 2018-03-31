@muladd function perform_step!(integrator,
                               cache::Union{ImplicitEMConstantCache,
                                            ImplicitEulerHeunConstantCache,
                                            ImplicitRKMilConstantCache},
                                            f=integrator.f)
  @unpack t,dt,uprev,u,p = integrator
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

  iter = 0
  κ = cache.κ
  tol = cache.tol

  L = integrator.g(uprev,p,t)
  ftmp = integrator.f(uprev,p,t)
  gtmp = L.*integrator.W.dW

  if typeof(cache) <: ImplicitEulerHeunConstantCache
    utilde = uprev + gtmp
    gtmp = ((integrator.g(utilde,p,t) + L)/2)*integrator.W.dW
  end

  if typeof(cache) <: ImplicitRKMilConstantCache || integrator.opts.adaptive == true
    if alg_interpretation(integrator.alg) == :Ito ||
       typeof(cache) <: ImplicitEMConstantCache
      K = @muladd uprev .+ dt.*ftmp
      utilde =  K + L*integrator.sqdt
      ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
      mil_correction = ggprime .* (integrator.W.dW.^2 .- dt)./2
      gtmp += mil_correction
    elseif alg_interpretation(integrator.alg) == :Stratonovich ||
           typeof(cache) <: ImplicitEulerHeunConstantCache
      utilde = uprev + L*integrator.sqdt
      ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
      mil_correction = ggprime.*(integrator.W.dW.^2)./2
      gtmp += mil_correction
    end
  end

  if integrator.alg.symplectic
    z = zero(u) # constant extrapolation, justified by ODE IM
  else
    z = dt*ftmp # linear extrapolation
  end

  iter += 1
  if integrator.alg.symplectic
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    u = uprev + z/2 + gtmp/2
  else
    u = uprev + dt*(1-theta)*ftmp + theta*z + gtmp
  end
  b = -z .+ dt.*f(u,p,t+a)
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
      u = uprev + z/2 + gtmp/2
    else
      u = uprev + dt*(1-theta)*ftmp + theta*z + gtmp
    end
    b = -z .+ dt.*f(u,p,t+a)
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
    u = uprev + z + gtmp
  else
    u = uprev + dt*(1-theta)*ftmp + theta*z + gtmp
  end

  if (iter >= integrator.alg.max_newton_iter && do_newton) || fail_convergence
    integrator.force_stepfail = true
    return
  end

  cache.ηold = η
  cache.newton_iters = iter
  u = uprev + dt*(1-theta)*ftmp + theta*z + gtmp

  if integrator.opts.adaptive

    Ed = dt*J*ftmp/2
    if typeof(cache) <: Union{ImplicitEMConstantCache,ImplicitEulerHeunConstantCache}
        En = mil_correction
    else
        En = integrator.opts.internalnorm.(dW.^3) .*
             integrator.opts.internalnorm.(ggprime).^2 ./ 6
    end

    tmp = Ed+En
    integrator.EEst = integrator.opts.internalnorm(tmp./(integrator.opts.abstol + max.(integrator.opts.internalnorm.(uprev),integrator.opts.internalnorm.(u))*integrator.opts.reltol))
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,
                               cache::Union{ImplicitEMCache,
                                            ImplicitEulerHeunCache,
                                            ImplicitRKMilCache},
                               f=integrator.f)
  @unpack t,dt,uprev,u,p = integrator
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

  if has_invW(f)
    f(Val{:invW},W,uprev,p,dt,t) # W == inverse W
  else
    if !integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < integrator.alg.new_jac_conv_bound
      new_jac = false
    else # Compute a new Jacobian
      new_jac = true
      if has_jac(f)
        f(Val{:jac},J,uprev,p,t)
      else
        uf.t = t
        jacobian!(J, uf, uprev, du1, integrator, jac_config)
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

  integrator.g(gtmp,uprev,p,t)
  integrator.f(tmp,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @tight_loop_macros for i in eachindex(u)
      @inbounds gtmp2[i]=gtmp[i]*dW[i]
    end
  else
    A_mul_B!(gtmp2,gtmp,dW)
  end

  if typeof(cache) <: ImplicitEulerHeunCache
    gtmp3 = cache.gtmp3
    @. z = uprev + gtmp2
    integrator.g(gtmp3,z,p,t)
    @. gtmp = (gtmp3 + gtmp)/2
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
      @. z = @muladd uprev + dt*tmp + gtmp*integrator.sqdt
      integrator.g(gtmp3,z,p,t)
      @. gtmp3 = (gtmp3-gtmp)/(integrator.sqdt) # ggprime approximation
      @. gtmp2 += gtmp3*(dW.^2 - dt)/2
    elseif alg_interpretation(integrator.alg) == :Stratonovich
      @. z = @muladd uprev + gtmp*integrator.sqdt
      integrator.g(gtmp3,z,p,t)
      @. gtmp3 = (gtmp3-gtmp)/(integrator.sqdt) # ggprime approximation
      @. gtmp2 += gtmp3*(dW.^2)/2
    end
  end

  ##############################################################################

  if integrator.alg.symplectic
    @. z = zero(u) # Justified by ODE solvers, constrant extrapolation when IM
  else
    @. z = dt*tmp # linear extrapolation
  end

  iter = 0
  κ = cache.κ
  tol = cache.tol
  if integrator.alg.symplectic
    @. u = uprev + z/2 + gtmp2/2
  else
    @. u = uprev + dt*(1-theta)*tmp + theta*z + gtmp2
  end
  iter += 1
  f(k,u,p,t+a)
  scale!(k,dt)
  if mass_matrix == I
    k .-= z
  else
    A_mul_B!(vec(du1),mass_matrix,vec(z))
    k .-= du1
  end
  if has_invW(f)
    A_mul_B!(vec(dz),W,vec(k)) # Here W is actually invW
  else
    cache.linsolve(vec(dz),W,vec(k),new_W)
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
    f(k,u,p,t+a)
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
      cache.linsolve(vec(dz),W,vec(k),false)
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

  if integrator.opts.adaptive

    if has_invW(f)
      # This means the Jacobian was never computed!
      f(Val{:jac},J,uprev,p,t)
    end

    A_mul_B!(vec(z),J,vec(tmp))
    @. k = dt*dt*z/2

    # k is Ed
    # dz is En
    if typeof(cache) <: Union{ImplicitEMCache,ImplicitEulerHeunCache}

      if !is_diagonal_noise(integrator.sol.prob)
        g_sized = norm(gtmp,2)
      else
        g_sized = gtmp
      end

      if typeof(cache) <: ImplicitEMCache
        @. z = @muladd uprev + dt*tmp + g_sized*integrator.sqdt

        if !is_diagonal_noise(integrator.sol.prob)
          integrator.g(gtmp,z,p,t)
          g_sized2 = norm(gtmp,2)
        else
          integrator.g(gtmp2,z,p,t)
          g_sized2 = gtmp2
        end

        @. dz = (g_sized2-g_sized)/(2integrator.sqdt)*(dW.^2 - dt)
      elseif typeof(cache) <: ImplicitEulerHeunCache
        @. z = @muladd uprev + g_sized*integrator.sqdt

        if !is_diagonal_noise(integrator.sol.prob)
          integrator.g(gtmp,z,p,t)
          g_sized2 = norm(gtmp,2)
        else
          integrator.g(gtmp2,z,p,t)
          g_sized2 = gtmp2
        end

        @. dz = (g_sized2-g_sized)/(2integrator.sqdt)*(dW.^2)
      end

    elseif typeof(cache) <: ImplicitRKMilCache
      # gtmp3 is ggprime
      @. dz = abs(dW^3)*integrator.opts.internalnorm(gtmp3)^2 / 6
    end

    @tight_loop_macros for (i,atol,rtol,δ) in zip(eachindex(u),Iterators.cycle(integrator.opts.abstol),
                            Iterators.cycle(integrator.opts.reltol),Iterators.cycle(integrator.opts.delta))
      @inbounds tmp[i] = (k[i]+dz[i])/(atol + max(integrator.opts.internalnorm(uprev[i]),integrator.opts.internalnorm(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(tmp)

  end

end
