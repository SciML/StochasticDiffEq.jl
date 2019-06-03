mutable struct SROCK_1ConstantCache{zType,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  ms::SVector{10,Int}
  mη::SVector{10,uEltypeNoUnits}
  zprev::zType
  mdeg::Int
  optimal_η::uEltypeNoUnits
end
@cache struct SROCK_1Cache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uᵢ₋₁::uType
  k::rateType
  gₘ₋₂::rateType
  gₘ₋₁::rateType
  tmp::uType
  constantcache::SROCK_1ConstantCache
end

function alg_cache(alg::SROCK_1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  SROCK_1ConstantCache(uEltypeNoUnits,u)
end

function alg_cache(alg::SROCK_1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  gₘ₋₁ = zero(rate_prototype)
  gₘ₋₂ = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  tmp  = zero(u)
  constantcache = SROCK_1ConstantCache(uEltypeNoUnits,u)
  SROCK_1Cache(u,uprev,uᵢ₋₁,k,gₘ₋₁,gₘ₋₂,tmp,constantcache)
end
