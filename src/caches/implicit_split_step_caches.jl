@cache mutable struct ISSEMCache{uType,rateType,N,noiseRateType,randType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  gtmp::noiseRateType
  gtmp2::rateType
  nlsolver::N
  dW_cache::randType
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = alg.theta,zero(t)
  nlsolver = OrdinaryDiffEq.build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)
  gtmp = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
    gtmp2 = gtmp
    dW_cache = nothing
  else
    gtmp2 = zero(rate_prototype)
    dW_cache = zero(ΔW)
  end

  ISSEMCache(u,uprev,fsalfirst,gtmp,gtmp2,nlsolver,dW_cache)
end

mutable struct ISSEMConstantCache{N} <: StochasticDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = alg.theta,zero(t)
  nlsolver = OrdinaryDiffEq.build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  ISSEMConstantCache(nlsolver)
end

@cache mutable struct ISSEulerHeunCache{uType,rateType,N,noiseRateType,randType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  gtmp::noiseRateType
  gtmp2::rateType
  gtmp3::noiseRateType
  nlsolver::N
  dW_cache::randType
end

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = alg.theta,zero(t)
  nlsolver = OrdinaryDiffEq.build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(true))
  fsalfirst = zero(rate_prototype)

  gtmp = zero(noise_rate_prototype)
  gtmp2 = zero(rate_prototype)

  if is_diagonal_noise(prob)
      gtmp3 = gtmp2
      dW_cache = nothing
  else
      gtmp3 = zero(noise_rate_prototype)
      dW_cache = zero(ΔW)
  end

  ISSEulerHeunCache(u,uprev,fsalfirst,gtmp,gtmp2,gtmp3,nlsolver,dW_cache)
end

mutable struct ISSEulerHeunConstantCache{N} <: StochasticDiffEqConstantCache
  nlsolver::N
end

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  γ, c = alg.theta,zero(t)
  nlsolver = OrdinaryDiffEq.build_nlsolver(alg,u,uprev,p,t,dt,f,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,γ,c,Val(false))
  ISSEulerHeunConstantCache(nlsolver)
end
