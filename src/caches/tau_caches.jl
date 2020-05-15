struct TauLeapingConstantCache <: StochasticDiffEqConstantCache end

@cache struct TauLeapingCache{uType,rateType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
end

alg_cache(alg::TauLeaping,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = TauLeapingConstantCache()

function alg_cache(alg::TauLeaping,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u)
  TauLeapingCache(u,uprev,tmp)
end
