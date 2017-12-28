mutable struct RackKenCarpCacheConstantCache{F,uEltypeNoUnits,Tab} <: OrdinaryDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
  tab::Tab
end

function alg_cache(alg::RackKenCarp,u,rate_prototype,uEltypeNoUnits,tTypeNoUnits,uBottomEltypeNoUnits,
                   uprev,uprev2,f,t,dt,reltol,::Type{Val{false}})
  uf = DiffEqDiffTools.UDerivativeWrapper(f,t)
  ηold = one(uEltypeNoUnits)

  if alg.κ != nothing
    κ = alg.κ
  else
    κ = uEltypeNoUnits(1//100)
  end
  if alg.tol != nothing
    tol = alg.tol
  else
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = KenCarp3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  RackKenCarpCacheConstantCache(uf,ηold,κ,tol,10000,tab)
end

mutable struct RackKenCarpCache{uType,rateType,uNoUnitsType,J,UF,JC,uEltypeNoUnits,Tab,F} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
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
end

u_cache(c::RackKenCarpCache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.dz)
du_cache(c::RackKenCarpCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::RackKenCarp,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z₁ = similar(u,indices(u)); z₂ = similar(u,indices(u));
  z₃ = similar(u,indices(u)); z₄ = similar(u,indices(u))
  dz = similar(u,indices(u))
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  tmp = similar(u); b = similar(u,indices(u));
  atmp = similar(u,uEltypeNoUnits,indices(u))

  uf = DiffEqDiffTools.UJacobianWrapper(f,t)
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
    tol = min(0.03,first(reltol)^(0.5))
  end

  tab = KenCarp3Tableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  RackKenCarpCache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(uf),
              typeof(jac_config),uEltypeNoUnits,typeof(tab),typeof(linsolve)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab)
end
