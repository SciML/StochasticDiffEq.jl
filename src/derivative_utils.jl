function calc_J!(integrator, cache, is_compos)
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du1,uf,J,jac_config = cache
    if has_jac(f)
      f(Val{:jac}, J, uprev, p, t)
    else
      uf.t = t
      uf.p = p
      jacobian!(J, uf, uprev, du1, integrator, jac_config)
      if is_compos
        integrator.eigen_est = norm(J, Inf)
      end
    end
end

function calc_W!(integrator, cache::StochasticDiffEqMutableCache, γdt, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p, alg = integrator
    @unpack J,W,jac_config = cache
    is_compos = is_composite(alg)
    alg = unwrap_alg(integrator, true)
    mass_matrix = integrator.sol.prob.mass_matrix

    new_W = true
    if has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && f(Val{:invW},W,uprev,p,γdt,t) # W == inverse W
      is_compos && calc_J!(integrator, cache, true)
    else
      # skip calculation of J if step is repeated
      if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < alg.new_jac_conv_bound)
        new_jac = false
      else # Compute a new Jacobian
        new_jac = true
        calc_J!(integrator, cache, is_compos)
      end
      # skip calculation of W if step is repeated
      if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t)))
        for j in 1:length(u), i in 1:length(u)
            @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
        end
      else
        new_W = false
      end
    end
    return new_W
  end
end

function calc_W!(integrator, cache::StochasticDiffEqConstantCache, γdt, repeat_step)
  uprev = integrator.uprev
  uf = cache.uf
  is_compos = is_composite(integrator.alg)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    is_compos && ( integrator.eigen_est = norm(J, Inf) )
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
    is_compos && ( integrator.eigen_est = J )
  end
  J, W
end
