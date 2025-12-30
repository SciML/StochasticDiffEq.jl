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

# ThetaTrapezoidalTauLeaping cache with theta parameter and extra storage for predictor step
struct ThetaTrapezoidalTauLeapingConstantCache{T} <: StochasticDiffEqConstantCache
    theta::T
    alpha1::T
    alpha2::T
end

@cache struct ThetaTrapezoidalTauLeapingCache{uType, rateType, T} <:
              StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    predictor::uType
    predictor_counts::rateType
    corrector_counts::rateType
    rate_at_predictor::rateType
    rate_at_uprev::rateType
    corrector_rate::rateType
    theta::T
    alpha1::T
    alpha2::T
end

function alg_cache(alg::ThetaTrapezoidalTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    theta = alg.theta
    alpha1 = one(theta) / (2 * (one(theta) - theta) * theta)
    alpha2 = ((one(theta) - theta)^2 + theta^2) / (2 * (one(theta) - theta) * theta)
    ThetaTrapezoidalTauLeapingConstantCache(theta, alpha1, alpha2)
end

function alg_cache(alg::ThetaTrapezoidalTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    theta = alg.theta
    alpha1 = one(theta) / (2 * (one(theta) - theta) * theta)
    alpha2 = ((one(theta) - theta)^2 + theta^2) / (2 * (one(theta) - theta) * theta)
    tmp = zero(u)
    predictor = zero(u)
    predictor_counts = zero(jump_rate_prototype)
    corrector_counts = zero(jump_rate_prototype)
    rate_at_predictor = zero(jump_rate_prototype)
    rate_at_uprev = zero(jump_rate_prototype)
    corrector_rate = zero(jump_rate_prototype)
    ThetaTrapezoidalTauLeapingCache(u, uprev, tmp, predictor, predictor_counts,
        corrector_counts, rate_at_predictor, rate_at_uprev,
        corrector_rate, theta, alpha1, alpha2)
end
