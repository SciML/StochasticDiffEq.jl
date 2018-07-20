struct LambaEMConstantCache <: StochasticDiffEqConstantCache end
struct LambaEMCache{uType,rateType,rateNoiseType,dWType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  du2::rateType
  K::rateType
  tmp::uType
  L::rateNoiseType
  gtmp::rateNoiseType
  dW_cache::dWType
end

u_cache(c::LambaEMCache) = ()
du_cache(c::LambaEMCache) = (c.du1,c.du2,c.K,c.L)

alg_cache(alg::LambaEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = LambaEMConstantCache()

function alg_cache(alg::LambaEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  du1 = zero(rate_prototype); du2 = zero(rate_prototype)
  K = zero(rate_prototype); tmp = zero(u);
  L = zero(noise_rate_prototype)
  gtmp = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
    dW_cache = nothing
  else
    dW_cache = zero(ΔW)
  end
  LambaEMCache(u,uprev,du1,du2,K,tmp,L,gtmp,dW_cache)
end

struct LambaEulerHeunConstantCache <: StochasticDiffEqConstantCache end
struct LambaEulerHeunCache{uType,rateType,rateNoiseType,dWType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  du2::rateType
  K::rateType
  tmp::uType
  L::rateNoiseType
  gtmp::rateNoiseType
  dW_cache::dWType
end

u_cache(c::LambaEulerHeunCache) = ()
du_cache(c::LambaEulerHeunCache) = (c.du1,c.du2,c.K,c.L)

alg_cache(alg::LambaEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = LambaEulerHeunConstantCache()

function alg_cache(alg::LambaEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  du1 = zero(rate_prototype); du2 = zero(rate_prototype)
  K = zero(rate_prototype); tmp = zero(u);
  L = zero(noise_rate_prototype)
  gtmp = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
    dW_cache = nothing
  else
    dW_cache = zero(ΔW)
  end
  LambaEulerHeunCache(u,uprev,du1,du2,K,tmp,L,gtmp,dW_cache)
end
