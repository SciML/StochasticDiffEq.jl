struct SROCK_1ConstantCache <: StochasticDiffEqConstantCache end
@cache struct SROCK_1Cache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::SROCK_1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = SROCK_1ConstantCache()

function alg_cache(alg::SROCK_1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  SROCK_1Cache(u,uprev,k,k₁,tmp)
end
