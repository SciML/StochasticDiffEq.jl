struct TauLeapingConstantCache <: StochasticDiffEqConstantCache end

@cache struct TauLeapingCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  newrate::rateType
  EEstcache::rateType
end

alg_cache(alg::TauLeaping,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = TauLeapingConstantCache()

function alg_cache(alg::TauLeaping,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u)
  newrate = zero(jump_rate_prototype)
  EEstcache = zero(jump_rate_prototype)
  TauLeapingCache(u,uprev,tmp,newrate,EEstcache)
end

alg_cache(alg::CaoTauLeaping,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = TauLeapingConstantCache()

function alg_cache(alg::CaoTauLeaping,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u)
  TauLeapingCache(u,uprev,tmp,nothing,nothing)
end
