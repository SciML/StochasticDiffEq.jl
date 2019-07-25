@cache mutable struct ISSEMCache{uType,rateType,N,noiseRateType,randType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  gtmp::noiseRateType
  gtmp2::rateType
  nlsolver::N
  dW_cache::randType
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  γ, c = alg.theta,zero(t)
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
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

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  γ, c = alg.theta,zero(t)
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
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

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  γ, c = alg.theta,zero(t)
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
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

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  γ, c = alg.theta,zero(t)
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  ISSEulerHeunConstantCache(nlsolver)
end
