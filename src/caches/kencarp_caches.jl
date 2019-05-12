mutable struct SKenCarpConstantCache{F,N,Tab} <: StochasticDiffEqConstantCache
  uf::F
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  tab = SKenCarpTableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.c3
  @oopnlsolve
  if uf !== nothing && typeof(f) <: SplitSDEFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
  SKenCarpConstantCache(uf,nlsolver,tab)
end

@cache mutable struct SKenCarpCache{uType,rateType,uNoUnitsType,JType,WType,UF,JC,N,Tab,F,kType,randType,rateNoiseType} <: StochasticDiffEqMutableCache
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
  J::JType
  W::WType
  uf::UF
  jac_config::JC
  linsolve::F
  nlsolver::N
  tab::Tab
  chi2::randType
  g1::rateNoiseType
  g4::rateNoiseType
end

u_cache(c::SKenCarpCache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.dz)
du_cache(c::SKenCarpCache)   = (c.k,c.fsalfirst)

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tab = SKenCarpTableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.c3
  @iipnlsolve

  atmp = fill!(similar(u,uEltypeNoUnits),0)
  z₁ = similar(u); z₂ = similar(u)
  z₃ = similar(u); z₄ = z
  dz = similar(u)
  if typeof(f) <: SplitSDEFunction
    k1 = zero(u); k2 = zero(u)
    k3 = zero(u); k4 = zero(u)
    uf = DiffEqDiffTools.UJacobianWrapper(f.f1,t,p)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
    uf = DiffEqDiffTools.UJacobianWrapper(f,t,p)
  end
  linsolve = alg.linsolve(Val{:init},uf,u)
  jac_config = build_jac_config(alg,f,uf,du1,uprev,u,tmp,dz)

  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end

  g1 = zero(noise_rate_prototype); g4 = zero(noise_rate_prototype)

  SKenCarpCache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
                typeof(jac_config),typeof(nlsolver),typeof(tab),typeof(linsolve),typeof(k1),
              typeof(chi2),typeof(g1)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,k1,k2,k3,k4,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,nlsolver,tab,chi2,g1,g4)
end
