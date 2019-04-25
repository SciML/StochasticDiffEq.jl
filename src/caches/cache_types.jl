abstract type StochasticDiffEqCache <: DECache end
abstract type StochasticDiffEqConstantCache <: StochasticDiffEqCache end
abstract type StochasticDiffEqMutableCache <: StochasticDiffEqCache end

mutable struct StochasticCompositeCache{T,F} <: StochasticDiffEqCache
  caches::T
  choice_function::F
  current::Int
end

function alg_cache(alg::algType,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{T}}) where {T,algType<:StochasticCompositeAlgorithm}
  caches = map((x)->alg_cache(x,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,Val{T}),alg.algs)
  StochasticCompositeCache(caches,alg.choice_function,1)
end
