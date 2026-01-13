struct TauLeapingConstantCache <: StochasticDiffEqConstantCache end

@cache struct TauLeapingCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    newrate::rateType
    EEstcache::rateType
end

function alg_cache(
        alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return TauLeapingConstantCache()
end

function alg_cache(
        alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    newrate = zero(jump_rate_prototype)
    EEstcache = zero(jump_rate_prototype)
    return TauLeapingCache(u, uprev, tmp, newrate, EEstcache)
end

function alg_cache(
        alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return TauLeapingConstantCache()
end

function alg_cache(
        alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    return TauLeapingCache(u, uprev, tmp, nothing, nothing)
end

# ThetaTrapezoidalTauLeaping cache
# Uses NLFunctional-style fixed-point iteration adapted for tau-leaping
# (The standard nlsolver infrastructure assumes ODE mass_matrix which DiscreteProblem lacks)

struct ThetaTrapezoidalTauLeapingConstantCache{T, N} <: StochasticDiffEqConstantCache
    theta::T
    nlalg::N  # NLFunctional algorithm with κ, max_iter settings
end

@cache struct ThetaTrapezoidalTauLeapingCache{uType, rateType, T, N} <:
              StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType               # Explicit contribution: X_n + ν*k - θ*dt*drift(X_n)
    z::uType                 # Current iteration value
    ztmp::uType              # New iteration result
    k::uType                 # Function evaluation storage
    drift_at_uprev::uType    # drift(X_n) = ν*a(X_n)
    atmp::uType              # Residual storage for convergence check
    poisson_counts::rateType # k ~ Poisson(dt * a(X_n))
    rate_at_uprev::rateType  # a(X_n)
    rate_tmp::rateType       # Temporary storage for rates
    theta::T
    nlalg::N
end

function alg_cache(alg::ThetaTrapezoidalTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ThetaTrapezoidalTauLeapingConstantCache(alg.theta, alg.nlsolve)
end

function alg_cache(alg::ThetaTrapezoidalTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    z = zero(u)
    ztmp = zero(u)
    k = zero(u)
    drift_at_uprev = zero(u)
    atmp = zero(u)
    poisson_counts = zero(jump_rate_prototype)
    rate_at_uprev = zero(jump_rate_prototype)
    rate_tmp = zero(jump_rate_prototype)
    return ThetaTrapezoidalTauLeapingCache(u, uprev, tmp, z, ztmp, k, drift_at_uprev, atmp,
        poisson_counts, rate_at_uprev, rate_tmp, alg.theta, alg.nlsolve)
end
