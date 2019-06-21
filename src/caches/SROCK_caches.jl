mutable struct SROCK1ConstantCache{zType,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  ms::SVector{10,Int}
  mη::SVector{10,uEltypeNoUnits}
  zprev::zType
  mdeg::Int
  optimal_η::uEltypeNoUnits
end
@cache struct SROCK1Cache{uType,rateType,noiseRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  gₘ₋₁::noiseRateType
  gₘ₋₂::noiseRateType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  atmp::rateType
  constantcache::SROCK1ConstantCache
end

function alg_cache(alg::SROCK1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  SROCK1ConstantCache{uEltypeNoUnits}(u)
end

function alg_cache(alg::SROCK1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  gₘ₋₁ = zero(noise_rate_prototype)
  gₘ₋₂ = zero(noise_rate_prototype)
  tmp  = uᵢ₋₂             # these 3 variables are dummied to use same memory
  fsalfirst = k
  atmp = zero(rate_prototype)
  constantcache = SROCK1ConstantCache{uEltypeNoUnits}(u)
  SROCK1Cache(u,uprev,uᵢ₋₁,uᵢ₋₂,gₘ₋₁,gₘ₋₂,tmp,k,fsalfirst,atmp,constantcache)
end

mutable struct SROCK2ConstantCache{zType,T} <: StochasticDiffEqConstantCache
  ms::SVector{46,Int}
  recf::Vector{T}
  mσ::SVector{46,T}
  mτ::SVector{46,T}
  recf2::Vector{T}
  mα::SVector{46,T}
  zprev::zType
  mdeg::Int
  deg_index::Int
  start::Int
end

@cache struct SROCK2Cache{uType,rateType,noiseRateType,T} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uᵢ::uType
  uₓ::uType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  Gₛ::noiseRateType
  Gₛ₁::noiseRateType
  vec_χ::T
  tmp::uType
  k::rateType
  fsalfirst::rateType
  atmp::rateType
  constantcache::SROCK2ConstantCache
end

function alg_cache(alg::SROCK2,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  SROCK2ConstantCache{uEltypeNoUnits}(u)
end

function alg_cache(alg::SROCK2,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  uᵢ = zero(u)
  uₓ = zero(u)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  Gₛ = zero(noise_rate_prototype)
  Gₛ₁ = zero(noise_rate_prototype)
  vec_χ = false .* vec(ΔW)
  tmp  = uᵢ₋₂            # these 2 variables are dummied to use same memory
  fsalfirst = k
  atmp = zero(rate_prototype)
  constantcache = SROCK2ConstantCache{uEltypeNoUnits}(u)
  SROCK2Cache{typeof(u),typeof(k),typeof(noise_rate_prototype),typeof(vec_χ)}(u,uprev,uᵢ,uₓ,uᵢ₋₁,uᵢ₋₂,Gₛ,Gₛ₁,vec_χ,tmp,k,fsalfirst,atmp,constantcache)
end

mutable struct SROCKEMConstantCache{zType,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  ms::SVector{10,Int}
  mη::SVector{10,uEltypeNoUnits}
  zprev::zType
  mdeg::Int
  optimal_η::uEltypeNoUnits
end

@cache struct SROCKEMCache{uType,rateType,noiseRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  Gₛ::noiseRateType
  Gₛ₁::noiseRateType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  atmp::rateType
  constantcache::SROCKEMConstantCache
end

function alg_cache(alg::SROCKEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  SROCKEMConstantCache{uEltypeNoUnits}(u)
end

function alg_cache(alg::SROCKEM,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  Gₛ = zero(noise_rate_prototype)
  if (!alg.strong_order_1 || is_diagonal_noise(prob) || typeof(ΔW) <: Number || length(ΔW) == 1)
    Gₛ₁ = Gₛ
  else
    Gₛ₁ = zero(noise_rate_prototype)
  end
  tmp  = zero(u)             # these 3 variables are dummied to use same memory
  fsalfirst = k
  atmp = zero(rate_prototype)
  constantcache = SROCKEMConstantCache{uEltypeNoUnits}(u)
  SROCKEMCache{typeof(u),typeof(k),typeof(Gₛ)}(u,uprev,uᵢ₋₁,uᵢ₋₂,Gₛ,Gₛ₁,tmp,k,fsalfirst,atmp,constantcache)
end

mutable struct SKSROCKConstantCache{zType,T} <: StochasticDiffEqConstantCache
  mc::Vector{T}
  mα::Vector{T}
  zprev::zType
end
@cache struct SKSROCKCache{uType,rateType,noise_rate_prototype} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  Gₛ::noise_rate_prototype
  tmp::uType
  k::rateType
  fsalfirst::rateType
  atmp::rateType
  constantcache::SKSROCKConstantCache
end

function alg_cache(alg::SKSROCK,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  SKSROCKConstantCache{uEltypeNoUnits}(u)
end

function alg_cache(alg::SKSROCK,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  Gₛ = zero(noise_rate_prototype)
  tmp  = uᵢ₋₂             # Dummmy variables
  fsalfirst = k
  atmp = zero(rate_prototype)
  constantcache = SKSROCKConstantCache{uEltypeNoUnits}(u)
  SKSROCKCache(u,uprev,uᵢ₋₁,uᵢ₋₂,Gₛ,tmp,k,fsalfirst,atmp,constantcache)
end

mutable struct TangXiaoSROCK2ConstantCache{zType,T} <: StochasticDiffEqConstantCache
  ms::SVector{46,Int}
  recf::Vector{T}
  mσ::SVector{46,T}
  mτ::SVector{46,T}
  recf2::Vector{T}
  mα::SVector{5, T}
  mn̂::SVector{5, Int}
  c1::SVector{13, T}
  c2::SVector{13, T}
  zprev::zType
  mdeg::Int
  deg_index::Int
  start::Int
  start_mcs::Int
end

@cache struct TangXiaoSROCK2Cache{uType,rateType,noiseRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType
  uᵢ::uType
  uₓ::uType
  Û₁::uType
  Û₂::uType
  uᵢ₋₁::uType
  uᵢ₋₂::uType
  Gₛ::noiseRateType
  Gₛ₁::noiseRateType
  tmp::uType
  k::rateType
  fsalfirst::rateType
  atmp::rateType
  constantcache::TangXiaoSROCK2ConstantCache
end

function alg_cache(alg::TangXiaoSROCK2,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  TangXiaoSROCK2ConstantCache{uEltypeNoUnits}(u)
end

function alg_cache(alg::TangXiaoSROCK2,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  k = zero(rate_prototype)
  uᵢ = zero(u)
  uₓ = zero(u)
  Û₁ = zero(u)
  Û₂ = zero(u)
  uᵢ₋₁ = zero(u)
  uᵢ₋₂ = zero(u)
  Gₛ = zero(noise_rate_prototype)
  if typeof(ΔW) <: Number || length(ΔW) == 1 || is_diagonal_noise(prob)
    Gₛ₁ = Gₛ
  else
    Gₛ₁ = zero(noise_rate_prototype)
  end
  tmp  = uᵢ₋₂            # these 2 variables are dummied to use same memory
  fsalfirst = k
  atmp = zero(rate_prototype)
  constantcache = TangXiaoSROCK2ConstantCache{uEltypeNoUnits}(u)
  TangXiaoSROCK2Cache{typeof(u),typeof(k),typeof(noise_rate_prototype)}(u,uprev,uᵢ,uₓ,Û₁,Û₂,uᵢ₋₁,uᵢ₋₂,Gₛ,Gₛ₁,tmp,k,fsalfirst,atmp,constantcache)
end
