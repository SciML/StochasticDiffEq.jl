mutable struct SKenCarpConstantCache{F,uEltypeNoUnits,Tab} <: StochasticDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  if typeof(f) <: SplitFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-1 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = SKenCarpTableau(real(uBottomEltype),real(tTypeNoUnits))

  SKenCarpConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct SKenCarpCache{uType,rateType,uNoUnitsType,J,UF,JC,uEltypeNoUnits,Tab,F,kType,randType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  dz::uType
  b::uType
  tmp::uType
  atmp::uNoUnitsType
  J::J
  W::J
  uf::UF
  jac_config::JC
  linsolve::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
  chi2::randType
  g1::rateNoiseType
  g4::rateNoiseType
end

u_cache(c::SKenCarpCache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.dz)
du_cache(c::SKenCarpCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = zeros(J)
  z₁ = similar(u,indices(u)); z₂ = similar(u,indices(u))
  z₃ = similar(u,indices(u)); z₄ = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  tmp = zeros(ΔW); b = similar(u,indices(u));
  atmp = zeros(u,uEltypeNoUnits,indices(u))

  if typeof(f) <: SplitFunction
    k1 = zeros(u); k2 = zeros(u)
    k3 = zeros(u); k4 = zeros(u)
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end
  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    reltol = 1e-1 # TODO: generalize
    tol = min(0.03,first(reltol)^(0.5))
  end

  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zeros(ΔW)
  end

  g1 = zeros(noise_rate_prototype); g4 = zeros(noise_rate_prototype)

  tab = SKenCarpTableau(real(uBottomEltype),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  SKenCarpCache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uEltypeNoUnits,typeof(tab),typeof(linsolve),typeof(k1),
              typeof(chi2),typeof(g1)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,k1,k2,k3,k4,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab,chi2,g1,g4)
end
