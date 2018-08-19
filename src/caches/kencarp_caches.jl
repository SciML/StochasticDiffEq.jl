mutable struct SKenCarpConstantCache{F,N,Tab} <: StochasticDiffEqConstantCache
  uf::F
  nlsolve::N
  tab::Tab
end

function alg_cache(alg::SKenCarp,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}})
  @oopnlcachefields
  if uf != nothing && typeof(f) <: SplitSDEFunction
    uf = DiffEqDiffTools.UDerivativeWrapper(f.f1,t,p)
  else
    uf = DiffEqDiffTools.UDerivativeWrapper(f,t,p)
  end
  tab = SKenCarpTableau(real(uBottomEltype),real(tTypeNoUnits))
  nlsolve = typeof(_nlsolve)(NLSolverCache(κ,tol,min_iter,max_iter,10000,new_W,z,W,tab.γ,tab.c3,ηold,z₊,dz,tmp,b,k))
  SKenCarpConstantCache(uf,nlsolve,tab)
end

mutable struct SKenCarpCache{uType,rateType,uNoUnitsType,J,W,UF,JC,uEltypeNoUnits,Tab,F,kType,randType,rateNoiseType} <: StochasticDiffEqMutableCache
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
  W::W
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

  du1 = zero(rate_prototype)
  if has_jac(f) && !has_invW(f) && f.jac_prototype != nothing
    W = WOperator(f, zero(t))
    J = nothing
  else
    J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
    W = similar(J)
  end
  z₁ = similar(u,axes(u)); z₂ = similar(u,axes(u))
  z₃ = similar(u,axes(u)); z₄ = similar(u,axes(u))
  dz = similar(u,axes(u))
  fsalfirst = zero(rate_prototype)
  k = zero(rate_prototype)
  tmp = zero(u); b = similar(u,axes(u));
  atmp = fill!(similar(u,uEltypeNoUnits,axes(u)),0)

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
    chi2 = zero(ΔW)
  end

  g1 = zero(noise_rate_prototype); g4 = zero(noise_rate_prototype)

  tab = SKenCarpTableau(real(uBottomEltype),real(tTypeNoUnits))

  ηold = one(uEltypeNoUnits)

  SKenCarpCache{typeof(u),typeof(rate_prototype),typeof(atmp),typeof(J),typeof(W),typeof(uf),
              typeof(jac_config),uEltypeNoUnits,typeof(tab),typeof(linsolve),typeof(k1),
              typeof(chi2),typeof(g1)}(
              u,uprev,du1,fsalfirst,k,z₁,z₂,z₃,z₄,k1,k2,k3,k4,dz,b,tmp,atmp,J,
              W,uf,jac_config,linsolve,ηold,κ,tol,10000,tab,chi2,g1,g4)
end
