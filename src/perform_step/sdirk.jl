@muladd function perform_step!(integrator,
                               cache::Union{ImplicitEMConstantCache,
                                            ImplicitEulerHeunConstantCache,
                                            ImplicitRKMilConstantCache})
  @unpack t,dt,uprev,u,p,P,c,f = integrator
  @unpack nlsolver = cache
  alg = unwrap_alg(integrator, true)
  OrdinaryDiffEq.markfirststage!(nlsolver)

  theta = alg.theta
  alg.symplectic ? a = dt/2 : a = theta*dt
  repeat_step = false

  # TODO: Stochastic extrapolants?
  u = uprev

  L = integrator.g(uprev,p,t)
  ftmp = integrator.f(uprev,p,t)
  gtmp = L.*integrator.W.dW

  if cache isa ImplicitEulerHeunConstantCache
    utilde = uprev + gtmp
    gtmp = ((integrator.g(utilde,p,t) + L)/2)*integrator.W.dW
  end

  if cache isa ImplicitRKMilConstantCache || integrator.opts.adaptive == true
    if alg_interpretation(alg) == :Ito ||
       cache isa ImplicitEMConstantCache
      K = @.. uprev + dt * ftmp
      utilde =  K + L*integrator.sqdt
      ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
      mil_correction = ggprime .* (integrator.W.dW.^2 .- abs(dt))./2
      gtmp += mil_correction
    elseif alg_interpretation(alg) == :Stratonovich ||
           cache isa ImplicitEulerHeunConstantCache
      utilde = uprev + L*integrator.sqdt
      ggprime = (integrator.g(utilde,p,t).-L)./(integrator.sqdt)
      mil_correction = ggprime.*(integrator.W.dW.^2)./2
      gtmp += mil_correction
    else
      error("Alg interpretation invalid. Use either :Ito or :Stratonovich")
    end
  end

  if alg.symplectic
    z = zero(u) # constant extrapolation, justified by ODE IM
  else
    z = dt*ftmp # linear extrapolation
  end
  nlsolver.z = z
  nlsolver.c = a

  if alg.symplectic
    # u = uprev + z then  u = (uprev+u)/2 = (uprev+uprev+z)/2 = uprev + z/2
    #u = uprev + z/2 + gtmp/2
    tmp = uprev + gtmp/2
  else
    #u = uprev + dt*(1-theta)*ftmp + theta*z + gtmp
    tmp = uprev + dt*(1-theta)*ftmp + gtmp
  end

  if P !== nothing
      ctmp = c(uprev, p, t, P.dW, nothing)
      tmp += ctmp
  end

  nlsolver.tmp = tmp
  z = OrdinaryDiffEq.nlsolve!(nlsolver, integrator, cache, repeat_step)
  OrdinaryDiffEq.nlsolvefail(nlsolver) && return nothing

  if alg.symplectic
    u = uprev + z + gtmp
  else
    #u = uprev + dt*(1-theta)*ftmp + theta*z + gtmp
    u = tmp + theta*z
  end

  if integrator.opts.adaptive

    if !OrdinaryDiffEq.isnewton(nlsolver)
      is_compos = isa(integrator.alg, StochasticDiffEqCompositeAlgorithm)
      J = OrdinaryDiffEq.calc_J(integrator, nlsolver.cache)
    end

    Ed = _reshape(dt*(nlsolver.cache.J*_vec(ftmp))/2, axes(ftmp))
    if cache isa Union{ImplicitEMConstantCache,ImplicitEulerHeunConstantCache}
        En = mil_correction
    else
        En = integrator.opts.internalnorm.(integrator.W.dW.^3, t) .*
             integrator.opts.internalnorm.(ggprime, t).^2 ./ 6
    end

    resids = calculate_residuals(Ed, En, uprev, u, integrator.opts.abstol,
                                 integrator.opts.reltol, integrator.opts.delta,
                                 integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(resids, t)
  end

  integrator.u = u
  return nothing
end

@muladd function perform_step!(integrator,
                               cache::Union{ImplicitEMCache,
                                            ImplicitEulerHeunCache,
                                            ImplicitRKMilCache})
  @unpack t,dt,uprev,u,p,P,c,f = integrator
  @unpack gtmp,gtmp2,nlsolver = cache
  @unpack z,tmp = nlsolver
  @unpack k,dz = nlsolver.cache # alias to reduce memory
  J = (OrdinaryDiffEq.isnewton(nlsolver) ? nlsolver.cache.J : nothing)
  alg = unwrap_alg(integrator, true)
  alg.symplectic ? a = dt/2 : a = alg.theta*dt
  dW = integrator.W.dW
  mass_matrix = integrator.f.mass_matrix
  theta = alg.theta

  OrdinaryDiffEq.markfirststage!(nlsolver)

  repeat_step = false

  if integrator.success_iter > 0 && !integrator.u_modified && alg.extrapolant == :interpolant
    current_extrapolant!(u,t+dt,integrator)
  elseif alg.extrapolant == :linear
    @.. u = uprev + integrator.fsalfirst*dt
  else # :constant
    copyto!(u,uprev)
  end

  ##############################################################################

  # Handle noise computations

  integrator.g(gtmp,uprev,p,t)
  integrator.f(tmp,uprev,p,t)

  if is_diagonal_noise(integrator.sol.prob)
    @.. gtmp2 = gtmp*dW
  else
    mul!(gtmp2,gtmp,dW)
  end

  if cache isa ImplicitEulerHeunCache
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

  if cache isa ImplicitRKMilCache
    gtmp3 = cache.gtmp3
    if alg_interpretation(alg) == :Ito
      @.. z = uprev + dt * tmp + integrator.sqdt * gtmp
      integrator.g(gtmp3,z,p,t)
      @.. gtmp3 = (gtmp3-gtmp)/(integrator.sqdt) # ggprime approximation
      @.. gtmp2 += gtmp3*(dW.^2 - abs(dt))/2
    elseif alg_interpretation(alg) == :Stratonovich
      @.. z = uprev + integrator.sqdt * gtmp
      integrator.g(gtmp3,z,p,t)
      @.. gtmp3 = (gtmp3-gtmp)/(integrator.sqdt) # ggprime approximation
      @.. gtmp2 += gtmp3*(dW.^2)/2
    else
      error("Alg interpretation invalid. Use either :Ito or :Stratonovich")
    end
  end

  if P !== nothing
    c(k, uprev, p, t, P.dW, nothing)
  end

  ##############################################################################

  if alg.symplectic
    @.. z = zero(eltype(u)) # Justified by ODE solvers, constant extrapolation when IM
  else
    @.. z = dt*tmp # linear extrapolation
  end

  nlsolver.c = a
  if alg.symplectic
    #@.. u = uprev + z/2 + gtmp2/2
    if P !== nothing
      @.. tmp = uprev + gtmp2/2 + k
    else
      @.. tmp = uprev + gtmp2/2
    end
  else
    #@.. u = uprev + dt*(1-theta)*tmp + theta*z + gtmp2
    if P !== nothing
      @.. tmp = uprev + dt*(1-theta)*tmp + gtmp2 + k
    else
      @.. tmp = uprev + dt*(1-theta)*tmp + gtmp2
    end
  end
  z = OrdinaryDiffEq.nlsolve!(nlsolver, integrator, cache, repeat_step)
  OrdinaryDiffEq.nlsolvefail(nlsolver) && return

  if alg.symplectic
    @.. u = uprev + z + gtmp2
  else
    @.. u = tmp + theta*z
  end

  if integrator.opts.adaptive

    if has_Wfact(f)
      # This means the Jacobian was never computed!
      f.jac(J,uprev,p,t)
    end

    mul!(vec(z),J,vec(tmp))
    @.. k = dt*dt*z/2

    # k is Ed
    # dz is En
    if cache isa Union{ImplicitEMCache,ImplicitEulerHeunCache}
      dW_cache = cache.dW_cache
      if !is_diagonal_noise(integrator.sol.prob)
        g_sized = norm(gtmp,2)
      else
        g_sized = gtmp
      end

      if cache isa ImplicitEMCache
        @.. z = uprev + dt * tmp + integrator.sqdt * g_sized

        if !is_diagonal_noise(integrator.sol.prob)
          integrator.g(gtmp,z,p,t)
          g_sized2 = norm(gtmp,2)
          @.. dW_cache = dW.^2 - dt
          diff_tmp = integrator.opts.internalnorm(dW_cache, t)
          En = (g_sized2-g_sized)/(2integrator.sqdt)*diff_tmp
          @.. dz = En
        else
          integrator.g(gtmp2,z,p,t)
          g_sized2 = gtmp2
          @.. dz = (g_sized2-g_sized)/(2integrator.sqdt)*(dW.^2 - dt)
        end

      elseif cache isa ImplicitEulerHeunCache
        @.. z = uprev + integrator.sqdt * g_sized

        if !is_diagonal_noise(integrator.sol.prob)
          integrator.g(gtmp,z,p,t)
          g_sized2 = norm(gtmp,2)
          @.. dW_cache = dW.^2
          diff_tmp = integrator.opts.internalnorm(dW_cache, t)
          En = (g_sized2-g_sized)/(2integrator.sqdt)*diff_tmp
          @.. dz = En
        else
          integrator.g(gtmp2,z,p,t)
          g_sized2 = gtmp2
          @.. dz = (g_sized2-g_sized)/(2integrator.sqdt)*(dW.^2)
        end

      end

    elseif cache isa ImplicitRKMilCache
      # gtmp3 is ggprime
      @.. dz = integrator.opts.internalnorm(dW^3, t)*integrator.opts.internalnorm(gtmp3, t)^2 / 6
    end

    calculate_residuals!(tmp, k, dz, uprev, u, integrator.opts.abstol,
                         integrator.opts.reltol, integrator.opts.delta,
                         integrator.opts.internalnorm, t)
    integrator.EEst = integrator.opts.internalnorm(tmp, t)

  end
end
