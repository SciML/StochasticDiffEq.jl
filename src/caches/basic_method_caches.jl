struct EMConstantCache <: StochasticDiffEqConstantCache end
@cache struct EMCache{uType,rateType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp1::rateType
  rtmp2::rateNoiseType
end

alg_cache(alg::EM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = EMConstantCache()

function alg_cache(alg::EM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u); rtmp1 = zero(rate_prototype);
  rtmp2 = zero(noise_rate_prototype)
  EMCache(u,uprev,tmp,rtmp1,rtmp2)
end

struct SplitEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct SplitEMCache{uType,rateType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp1::rateType
  rtmp2::rateNoiseType
end

alg_cache(alg::SplitEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = SplitEMConstantCache()

function alg_cache(alg::SplitEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u); rtmp1 = zero(rate_prototype);
  rtmp2 = zero(noise_rate_prototype)
  SplitEMCache(u,uprev,tmp,rtmp1,rtmp2)
end

struct EulerHeunConstantCache <: StochasticDiffEqConstantCache end
@cache struct EulerHeunCache{uType,rateType,rateNoiseType,rateNoiseCollectionType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  ftmp1::rateType
  ftmp2::rateType
  nrtmp::rateNoiseCollectionType
  gtmp1::rateNoiseType
  gtmp2::rateNoiseType
end

alg_cache(alg::EulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = EulerHeunConstantCache()

function alg_cache(alg::EulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u); ftmp1 = zero(rate_prototype); ftmp2 = zero(rate_prototype)
  nrtmp = zero(rate_prototype)
  gtmp1 = zero(noise_rate_prototype); gtmp2 = zero(noise_rate_prototype)
  EulerHeunCache(u,uprev,tmp,ftmp1,ftmp2,nrtmp,gtmp1,gtmp2)
end

struct RandomEMConstantCache <: StochasticDiffEqConstantCache end
@cache struct RandomEMCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  tmp::uType
  rtmp::rateType
end

alg_cache(alg::RandomEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = RandomEMConstantCache()

function alg_cache(alg::RandomEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tmp = zero(u); rtmp = zero(rate_prototype)
  RandomEMCache(u,uprev,tmp,rtmp)
end

struct RKMilConstantCache <: StochasticDiffEqConstantCache end
@cache struct RKMilCache{uType,rateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  du2::rateType
  K::rateType
  tmp::uType
  L::rateType
end

alg_cache(alg::RKMil,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = RKMilConstantCache()

function alg_cache(alg::RKMil,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  du1 = zero(rate_prototype); du2 = zero(rate_prototype)
  K = zero(rate_prototype); tmp = zero(u); L = zero(rate_prototype)
  RKMilCache(u,uprev,du1,du2,K,tmp,L)
end

struct RKMilCommuteConstantCache <: StochasticDiffEqConstantCache end
@cache struct RKMilCommuteCache{uType,rateType,rateNoiseType,WikType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  du2::rateType
  K::rateType
  gtmp::rateNoiseType
  L::rateNoiseType
  WikJ::WikType
  Dg::WikType
  mil_correction::rateType
  Kj::uType
  Dgj::rateNoiseType
  tmp::uType
end

alg_cache(alg::RKMilCommute,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = RKMilCommuteConstantCache()

function alg_cache(alg::RKMilCommute,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  du1 = zero(rate_prototype); du2 = zero(rate_prototype)
  K = zero(rate_prototype); gtmp = zero(noise_rate_prototype);
  L = zero(noise_rate_prototype); tmp = zero(rate_prototype)
  WikJ = false .* vec(ΔW) .* vec(ΔW)'
  Dg = false .* vec(ΔW) .* vec(ΔW)'
  mil_correction = zero(rate_prototype)
  Kj = zero(u); Dgj = zero(noise_rate_prototype)
  RKMilCommuteCache(u,uprev,du1,du2,K,gtmp,L,WikJ,Dg,mil_correction,Kj,Dgj,tmp)
end
