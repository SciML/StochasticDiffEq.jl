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
  nlsolver::N
  dW_cache::randType
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  γ, c = alg.theta,zero(t)
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields
  gtmp = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
    gtmp2 = gtmp
    dW_cache = nothing
  else
    gtmp2 = zero(rate_prototype)
    dW_cache = zero(ΔW)
  end

  ISSEMCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,J,W,jac_config,linsolve,uf,
                  nlsolver,dW_cache)
end

mutable struct ISSEMConstantCache{F,N} <: StochasticDiffEqConstantCache
  uf::F
  nlsolver::N
end

function alg_cache(alg::ISSEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  γ, c = alg.theta,zero(t)
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  uf = nlsolver.uf
  ISSEMConstantCache(uf,nlsolver)
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
  nlsolver::N
  dW_cache::randType
end

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  γ, c = alg.theta,zero(t)
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @getiipnlsolvefields

  gtmp = zero(noise_rate_prototype)
  gtmp2 = zero(rate_prototype)

  if is_diagonal_noise(prob)
      gtmp3 = gtmp2
      dW_cache = nothing
  else
      gtmp3 = zero(noise_rate_prototype)
      dW_cache = zero(ΔW)
  end

  ISSEulerHeunCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,gtmp3,
                         J,W,jac_config,linsolve,uf,nlsolver,dW_cache)
end

mutable struct ISSEulerHeunConstantCache{F,N} <: StochasticDiffEqConstantCache
  uf::F
  nlsolver::N
end

function alg_cache(alg::ISSEulerHeun,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  γ, c = alg.theta,zero(t)
  W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  uf = nlsolver.uf
  ISSEulerHeunConstantCache(uf,nlsolver)
end
