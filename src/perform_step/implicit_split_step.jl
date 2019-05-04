@muladd function perform_step!(integrator,
                               cache::Union{ISSEMConstantCache,
                                            ISSEulerHeunConstantCache},
                                            f=integrator.f)
  @unpack t,dt,uprev,u,p = integrator
  @unpack uf = cache
  nlsolve! = cache.nlsolve; nlcache = nlsolve!.cache
  alg = unwrap_alg(integrator, true)
  theta = alg.theta
  alg.symplectic ? a = dt/2 : a = dt

  # TODO: Stochastic extrapolants?
  u = uprev

  repeat_step = false
  if nlsolve! isa NLNewton
    uf.t = t
    J, nlcache.W = calc_W!(integrator, cache, dt*theta, repeat_step)
  end

  L = integrator.g(uprev,p,t)
  ftmp = integrator.f(uprev,p,t)

  if alg.symplectic
    z = zero(u) # constant extrapolation, justified by ODE IM
  else
    z = dt*ftmp # linear extrapolation
  end
  nlcache.z = z

  nlcache.c = a
  if alg.symplectic
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    #u = uprev + z/2
    tmp = uprev
  else
    tmp = uprev + dt*(1-theta)*ftmp
  end
  nlcache.tmp = tmp

  (z, η, iter, fail_convergence) = nlsolve!(integrator); fail_convergence && return nothing

  if alg.symplectic
    u = tmp + z
  else
    u = tmp + theta*z
  end

  gtmp = L.*integrator.W.dW

  if typeof(cache) <: ISSEulerHeunConstantCache
    utilde = u + gtmp
    gtmp = ((integrator.g(utilde,p,t) + L)/2)*integrator.W.dW
  end

  u += gtmp

  nlcache.ηold = η
  nlcache.nl_iters = iter

  if integrator.opts.adaptive

    Ed = dt*(J*ftmp)/2

    if typeof(cache) <: SplitStepEulerConstantCache
      K = @.. uprev + dt * ftmp
      utilde =  @.. K + integrator.sqdt * L
      ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
      En = ggprime .* (integrator.W.dW.^2 .- dt)./2
    elseif typeof(cache) <: ISSEulerHeunConstantCache
      utilde = uprev + L*integrator.sqdt
      ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
      En = ggprime.*(integrator.W.dW.^2)./2
    end

    resids = calculate_residuals(Ed, En, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(resids,t)
  end

  integrator.u = u
end

@muladd function perform_step!(integrator, cache::Union{ISSEMCache,
                                                        ISSEulerHeunCache},
                               f=integrator.f)
  @unpack t,dt,uprev,u,p = integrator
  @unpack uf,du1,dz,z,k,J,W,jac_config,gtmp,gtmp2,tmp,tmp,dW_cache = cache
  nlsolve! = cache.nlsolve; nlcache = nlsolve!.cache
  alg = unwrap_alg(integrator, true)
  alg.symplectic ? a = dt/2 : a = dt
  dW = integrator.W.dW
  mass_matrix = integrator.f.mass_matrix
  theta = alg.theta

  repeat_step = false

  if integrator.success_iter > 0 && !integrator.u_modified && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif alg.extrapolant == :linear
    @.. u = uprev + integrator.fsalfirst*dt
  else # :constant
    copyto!(u,uprev)
  end

  nlsolve! isa NLNewton && calc_W!(integrator, cache, dt, repeat_step)

  integrator.f(tmp,uprev,p,t)

  if alg.symplectic
    @.. z = zero(eltype(u)) # Justified by ODE solvers, constrant extrapolation when IM
  else
    @.. z = dt*tmp # linear extrapolation
  end

  if alg.symplectic
    #@.. u = uprev + z/2
    @.. tmp = uprev
  else
    #@.. u = uprev + dt*(1-theta)*tmp + theta*z
    @.. tmp = uprev + dt*(1-theta)*tmp
  end
  nlcache.c = a
  (z, η, iter, fail_convergence) = nlsolve!(integrator); fail_convergence && return nothing

  if alg.symplectic
    @.. u = uprev + z
  else
    #@.. u = uprev + dt*(1-theta)*tmp + theta*z
    @.. u = tmp + theta*z
  end

  nlcache.ηold = η
  nlcache.nl_iters = iter

  ##############################################################################

  # Handle noise computations

  integrator.g(gtmp,uprev,p,t)


  if is_diagonal_noise(integrator.sol.prob)
    @.. gtmp2 = gtmp*dW
  else
    mul!(gtmp2,gtmp,dW)
  end

  if typeof(cache) <: ISSEulerHeunCache
    gtmp3 = cache.gtmp3
    @.. z = uprev + gtmp2
    integrator.g(gtmp3,z,p,t)
    @.. gtmp = (gtmp3 + gtmp)/2
    if is_diagonal_noise(integrator.sol.prob)
      @.. gtmp2 = gtmp*dW
    else
      mul!(gtmp2,gtmp,dW)
    end
  end

  @.. u += gtmp2

  ##############################################################################

  if integrator.opts.adaptive

    if has_invW(f)
      # This means the Jacobian was never computed!
      f.jac(J,uprev,p,t)
    end

    mul!(vec(z),J,vec(tmp))
    @.. k = dt*dt*z/2

    # k is Ed
    # dz is En

      if !is_diagonal_noise(integrator.sol.prob)
        g_sized = norm(gtmp,2)
      else
        g_sized = gtmp
      end

      if typeof(cache) <: ISSEMCache
        @.. z = uprev + dt*tmp + integrator.sqdt * g_sized

        if !is_diagonal_noise(integrator.sol.prob)
          integrator.g(gtmp,z,p,t)
          g_sized2 = norm(gtmp,2)
          @.. dW_cache = dW.^2 - dt
          diff_tmp = integrator.opts.internalnorm(dW_cache,t)
          En = (g_sized2-g_sized)/(2integrator.sqdt)*diff_tmp
          @.. dz = En
        else
          integrator.g(gtmp2,z,p,t)
          g_sized2 = gtmp2
          @.. dz = (g_sized2-g_sized)/(2integrator.sqdt)*(dW.^2 - dt)
        end

      elseif typeof(cache) <: ISSEulerHeunCache
        @.. z = uprev + integrator.sqdt * g_sized

        if !is_diagonal_noise(integrator.sol.prob)
          integrator.g(gtmp,z,p,t)
          g_sized2 = norm(gtmp,2)
          @.. dW_cache = dW.^2
          diff_tmp = integrator.opts.internalnorm(dW_cache,t)
          En = (g_sized2-g_sized)/(2integrator.sqdt)*diff_tmp
          @.. dz = En
        else
          integrator.g(gtmp2,z,p,t)
          g_sized2 = gtmp2
          @.. dz = (g_sized2-g_sized)/(2integrator.sqdt)*(dW.^2)
        end


    end

    calculate_residuals!(tmp, k, dz, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(tmp,t)

  end
end
