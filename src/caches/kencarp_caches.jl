mutable struct SKenCarpConstantCache{N,Tab} <: StochasticDiffEqConstantCache
  nlsolver::N
  tab::Tab
end

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  tab = SKenCarpTableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.c3
  J, W = oop_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  SKenCarpConstantCache(nlsolver,tab)
end

@cache mutable struct SKenCarpCache{uType,rateType,uNoUnitsType,N,Tab,kType,randType,rateNoiseType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  fsalfirst::rateType
  z₁::uType
  z₂::uType
  z₃::uType
  z₄::uType
  k1::kType
  k2::kType
  k3::kType
  k4::kType
  atmp::uNoUnitsType
  nlsolver::N
  tab::Tab
  chi2::randType
  g1::rateNoiseType
  g4::rateNoiseType
end

u_cache(c::SKenCarpCache)    = (c.z₁,c.z₂,c.z₃,c.z₄,c.nlsolver.dz)
du_cache(c::SKenCarpCache)   = (c.nlsolver.k,c.fsalfirst)

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  tab = SKenCarpTableau(real(uBottomEltypeNoUnits),real(tTypeNoUnits))
  γ, c = tab.γ,tab.c3
  J, W = iip_generate_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits)
  nlsolver = iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  fsalfirst = zero(rate_prototype)

  atmp = fill!(similar(u,uEltypeNoUnits),0)
  z₁ = similar(u); z₂ = similar(u)
  z₃ = similar(u); z₄ = nlsolver.z
  if typeof(f) <: SplitSDEFunction
    k1 = zero(u); k2 = zero(u)
    k3 = zero(u); k4 = zero(u)
  else
    k1 = nothing; k2 = nothing
    k3 = nothing; k4 = nothing
  end

  if typeof(ΔW) <: Union{SArray,Number}
    chi2 = copy(ΔW)
  else
    chi2 = zero(ΔW)
  end

  g1 = zero(noise_rate_prototype); g4 = zero(noise_rate_prototype)

  SKenCarpCache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(nlsolver),
              typeof(tab),typeof(k1),typeof(chi2),typeof(g1)}(
              u,uprev,fsalfirst,z₁,z₂,z₃,z₄,k1,k2,k3,k4,atmp,nlsolver,tab,chi2,g1,g4)
end
