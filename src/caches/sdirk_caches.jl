mutable struct ImplicitEMCache{uType,rateType,J,JC,UF,uEltypeNoUnits} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uprev2::uType
  du1::rateType
  fsalfirst::rateType
  k::rateType
  z::uType
  dz::uType
  tmp::uType
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

function alg_cache(alg::ImplicitEM,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{true}})

  du1 = zeros(rate_prototype)
  J = zeros(uEltypeNoUnits,length(u),length(u)) # uEltype?
  W = similar(J)
  z = similar(u)
  dz = similar(u); tmp = similar(u)
  fsalfirst = zeros(rate_prototype)
  k = zeros(rate_prototype)
  vfr = VectorFReturn(f,size(u))
  uf = UJacobianWrapper(vfr,t)
  if alg_autodiff(alg)
    jac_config = ForwardDiff.JacobianConfig(uf,vec(du1),vec(uprev),
                    ForwardDiff.Chunk{determine_chunksize(u,alg)}())
  else
    jac_config = nothing
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
    tol = min(0.03,first(reltol)^(0.5))
  end
  ImplicitEMCache(u,uprev,uprev2,du1,fsalfirst,k,z,dz,tmp,J,W,jac_config,uf,ηold,κ,tol,10000)
end

mutable struct ImplicitEMConstantCache{F,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  uf::F
  ηold::uEltypeNoUnits
  κ::uEltypeNoUnits
  tol::uEltypeNoUnits
  newton_iters::Int
end

function alg_cache(alg::ImplicitEM,u,rate_prototype,uEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,reltol,::Type{Val{false}})
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
    tol = min(0.03,first(reltol)^(0.5))
  end

  ImplicitEMConstantCache(uf,ηold,κ,tol,100000)
end
