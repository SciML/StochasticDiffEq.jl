
function stepsize_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.q11 = DiffEqBase.value(FastPower.fastpower(integrator.EEst,controller.beta1))
    integrator.q = DiffEqBase.value(integrator.q11/FastPower.fastpower(integrator.qold,controller.beta2))
    @fastmath integrator.q = DiffEqBase.value(max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma)))
end

@inline function step_accept_controller!(integrator::SDEIntegrator, alg)
    step_accept_controller!(integrator, integrator.opts.controller, alg)
end

function step_accept_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.dtnew = DiffEqBase.value(integrator.dt/integrator.q) * oneunit(integrator.dt)
end

function step_reject_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
end


function stepsize_controller!(integrator::SDEIntegrator, alg::TauLeaping)
  nothing  # Post-leap adjustment happens in perform_step!
end

function step_accept_controller!(integrator::SDEIntegrator, alg::TauLeaping)
  if alg.adaptive
    integrator.q = min(integrator.opts.gamma / integrator.EEst, integrator.opts.qmax)
    return integrator.dt * integrator.q
  else
    return integrator.dt
  end
end

function step_reject_controller!(integrator::SDEIntegrator, alg::TauLeaping)
  if alg.adaptive
    integrator.dt = integrator.opts.gamma * integrator.dt / integrator.EEst
  end
end

# CaoTauLeaping: Pre-leap τ computation
function stepsize_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
  if !alg.adaptive
    return
  end
  
  @unpack u, p, t, P, opts, c = integrator
  cache = integrator.cache

  P === nothing && error("CaoTauLeaping requires a JumpProblem with a RegularJump")
  
  # Handle both constant and mutable caches
  if isa(cache, CaoTauLeapingConstantCache)
    rate = P.cache.rate(u, p, t)  # Compute propensities directly
    mu = zero(u)
    sigma2 = zero(u)
  else  # CaoTauLeapingCache
    @unpack mu, sigma2, rate = cache
    P.cache.rate(rate, u, p, t)  # Compute propensities into cache
    fill!(mu, zero(eltype(mu)))
    fill!(sigma2, zero(eltype(sigma2)))
  end

  # Infer ν_ij using c by applying unit counts for each reaction
  num_reactions = length(rate)
  ν = zeros(eltype(u), length(u), num_reactions)
  unit_counts = zeros(eltype(rate), num_reactions)
  for j in 1:num_reactions
    unit_counts[j] = 1
    c(ν[:, j], u, p, t, unit_counts, nothing)  # ν[:, j] is the change vector for reaction j
    unit_counts[j] = 0  # Reset
  end

  # Compute μ_i and σ_i^2
  for i in eachindex(u)
    for j in 1:num_reactions
      ν_ij = ν[i, j]
      mu[i] += ν_ij * rate[j]
      sigma2[i] += ν_ij^2 * rate[j]
    end
  end

  # Compute τ per species
  ϵ = alg.epsilon
  τ_vals = similar(u, Float64)
  for i in eachindex(u)
    max_term = max(ϵ * u[i], 1.0)
    τ1 = abs(mu[i]) > 0 ? max_term / abs(mu[i]) : Inf
    τ2 = sigma2[i] > 0 ? max_term^2 / sigma2[i] : Inf
    τ_vals[i] = min(τ1, τ2)
  end

  τ = min(minimum(τ_vals), opts.dtmax)
  integrator.dt = max(τ, opts.dtmin)
  integrator.EEst = 1.0
end

function step_accept_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
  return integrator.dt
end

function step_reject_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
  error("CaoTauLeaping should never reject steps")
end
