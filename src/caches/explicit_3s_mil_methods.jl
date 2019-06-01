struct WangLi3SMil_AConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_ACache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::WangLi3SMil_A,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = WangLi3SMil_AConstantCache()

function alg_cache(alg::WangLi3SMil_A,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  WangLi3SMil_ACache(u,uprev,k,k₁,tmp)
end

struct WangLi3SMil_BConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_BCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::WangLi3SMil_B,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = WangLi3SMil_BConstantCache()

function alg_cache(alg::WangLi3SMil_B,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  WangLi3SMil_BCache(u,uprev,k,k₁,tmp)
end

struct WangLi3SMil_CConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_CCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::WangLi3SMil_C,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = WangLi3SMil_CConstantCache()

function alg_cache(alg::WangLi3SMil_C,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  WangLi3SMil_CCache(u,uprev,k,k₁,tmp)
end

struct WangLi3SMil_DConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_DCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::WangLi3SMil_D,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = WangLi3SMil_DConstantCache()

function alg_cache(alg::WangLi3SMil_D,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  WangLi3SMil_DCache(u,uprev,k,k₁,tmp)
end

struct WangLi3SMil_EConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_ECache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::WangLi3SMil_E,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = WangLi3SMil_EConstantCache()

function alg_cache(alg::WangLi3SMil_E,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  WangLi3SMil_ECache(u,uprev,k,k₁,tmp)
end

struct WangLi3SMil_FConstantCache <: StochasticDiffEqConstantCache end
@cache struct WangLi3SMil_FCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  k::rateType
  k₁::rateType
  tmp::uType
end

alg_cache(alg::WangLi3SMil_F,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = WangLi3SMil_FConstantCache()

function alg_cache(alg::WangLi3SMil_F,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  k₁ = zero(rate_prototype)
  tmp = zero(u)
  WangLi3SMil_FCache(u,uprev,k,k₁,tmp)
end
