mutable struct ImplicitEMCache{uType,rateType,J,JC,UF,uEltypeNoUnits,noiseRateType} <: StochasticDiffEqMutableCache
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
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::ImplicitEMCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitEMCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitEM,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u)
  dz = similar(u); tmp = similar(u); gtmp = similar(noise_rate_prototype)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)

  uf = DiffEqDiffTools.UJacobianWrapper(f,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-5 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  if is_diagonal_noise(prob)
    gtmp2 = gtmp
  else
    gtmp2 = similar(rate_prototype)
  end

  ImplicitEMCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,J,W,jac_config,uf,
                  ηold,κ,tol,10000)
end

mutable struct ImplicitEMConstantCache{F,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

function alg_cache(alg::ImplicitEM,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-5 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  ImplicitEMConstantCache(uf,ηold,κ,tol,100000)
end

mutable struct ImplicitEulerHeunCache{uType,rateType,J,JC,UF,uEltypeNoUnits,noiseRateType} <: StochasticDiffEqMutableCache
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
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::ImplicitEulerHeunCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitEulerHeunCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitEulerHeun,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u)
  dz = similar(u); tmp = similar(u); gtmp = similar(noise_rate_prototype)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)

  uf = DiffEqDiffTools.UJacobianWrapper(f,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-5 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  gtmp2 = similar(rate_prototype)

  ImplicitEulerHeunCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,J,W,jac_config,uf,
                  ηold,κ,tol,10000)
end

mutable struct ImplicitEulerHeunConstantCache{F,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

function alg_cache(alg::ImplicitEulerHeun,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-5 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  ImplicitEulerHeunConstantCache(uf,ηold,κ,tol,100000)
end

mutable struct ImplicitRKMilCache{uType,rateType,J,JC,UF,uEltypeNoUnits,noiseRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  tmp::uType
  gtmp::noiseRateType
  gtmp2::noiseRateType
  gtmp3::noiseRateType
  J::J
  W::J
  jac_config::JC
  uf::UF
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

u_cache(c::ImplicitRKMilCache)    = (c.uprev2,c.z,c.dz)
du_cache(c::ImplicitRKMilCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::ImplicitRKMil,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u)
  dz = similar(u); tmp = similar(u); gtmp = similar(noise_rate_prototype)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)

  uf = DiffEqDiffTools.UJacobianWrapper(f,t)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-5 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  gtmp2 = similar(rate_prototype)
  gtmp3 = similar(rate_prototype)

  ImplicitRKMilCache(u,uprev,du1,fsalfirst,k,z,dz,tmp,gtmp,gtmp2,gtmp3,
                   J,W,jac_config,uf,ηold,κ,tol,10000)
end

mutable struct ImplicitRKMilConstantCache{F,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

function alg_cache(alg::ImplicitRKMil,prob,u,ΔW,ΔZ,rate_prototype,noise_rate_prototype,
                   uEltypeNoUnits,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  uf = UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-5 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  ImplicitRKMilConstantCache(uf,ηold,κ,tol,100000)
end
