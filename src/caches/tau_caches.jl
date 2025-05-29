@cache mutable struct TauLeapingCache{uType, rateType, jumpRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rate::rateType
  newrate::rateType
  EEstcache::jumpRateType
end

function alg_cache(alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt, ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
  tmp = zero(u)
  rate = zero(jump_rate_prototype)
  newrate = zero(jump_rate_prototype)
  EEstcache = zero(jump_rate_prototype)
  TauLeapingCache(u, uprev, tmp, rate, newrate, EEstcache)
end

@cache mutable struct CaoTauLeapingCache{uType, rateType, muType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rate::rateType
  mu::muType
  sigma2::muType
end

function alg_cache(alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt, ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
  tmp = zero(u)
  rate = zero(jump_rate_prototype)
  mu = zero(u)
  sigma2 = zero(u)
  CaoTauLeapingCache(u, uprev, tmp, rate, mu, sigma2)
end

struct TauLeapingConstantCache <: StochasticDiffEqConstantCache end
struct CaoTauLeapingConstantCache <: StochasticDiffEqConstantCache end

alg_cache(alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt, ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits} = TauLeapingConstantCache()
alg_cache(alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype, noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt, ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits} = CaoTauLeapingConstantCache()
