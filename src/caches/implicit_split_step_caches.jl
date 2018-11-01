@cache mutable struct ISSEMCache{uType,rateType,JType,WType,JC,UF,
                          N,noiseRateType,F,randType} <:
                          StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  tmp::uType
  gtmp::noiseRateType
  gtmp2::rateType
  J::JType
  W::WType
  jac_config::JC
  linsolve::F
  uf::UF
  nlsolve::N
  dW_cache::randType
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  @iipnlcachefields
  gtmp = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
    gtmp2 = gtmp
    dW_cache = nothing
  else
    gtmp2 = zero(rate_prototype)
    dW_cache = zero(ΔW)
  end

  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,100000,new_W,z,W,alg.theta,zero(t),ηold,z₊,dz,tmp,b,k))
  ISSEMCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,J,W,jac_config,linsolve,uf,
                  nlsolve,dW_cache)
end

mutable struct ISSEMConstantCache{F,N} <: StochasticDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,100000,new_W,z,W,alg.theta,zero(t),ηold,z₊,dz,tmp,b,k))
  ISSEMConstantCache(uf,nlsolve)
end

@cache mutable struct ISSEulerHeunCache{uType,rateType,JType,WType,JC,UF,N,
                                 noiseRateType,F,randType} <:
                                 StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  tmp::uType
  gtmp::noiseRateType
  gtmp2::rateType
  gtmp3::noiseRateType
  J::JType
  W::WType
  jac_config::JC
  linsolve::F
  uf::UF
  nlsolve::N
  dW_cache::randType
end

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  @iipnlcachefields

  gtmp = zero(noise_rate_prototype)
  gtmp2 = zero(rate_prototype)

  if is_diagonal_noise(prob)
      gtmp3 = gtmp2
      dW_cache = nothing
  else
      gtmp3 = zero(noise_rate_prototype)
      dW_cache = zero(ΔW)
  end

  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,100000,new_W,z,W,alg.theta,zero(t),ηold,z₊,dz,tmp,b,k))
  ISSEulerHeunCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,gtmp3,
                         J,W,jac_config,linsolve,uf,nlsolve,dW_cache)
end

mutable struct ISSEulerHeunConstantCache{F,N} <: StochasticDiffEqConstantCache
  uf::F
  nlsolve::N
end

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  @oopnlcachefields
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,100000,new_W,z,W,alg.theta,zero(t),ηold,z₊,dz,tmp,b,k))
  ISSEulerHeunConstantCache(uf,nlsolve)
end
