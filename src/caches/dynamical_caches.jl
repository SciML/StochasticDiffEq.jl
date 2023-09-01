abstract type StochasticDynamicalEqConstantCache <: StochasticDiffEqConstantCache end # Pourquoi faire ça, Si c'est pour avoir une seul function de check dans initialize!
abstract type StochasticDynamicalEqMutableCache <: StochasticDiffEqMutableCache end


mutable struct BAOABConstantCache{uType,uEltypeNoUnits,uCoeffType, uCoeffMType} <: StochasticDynamicalEqConstantCache
  k::uType
  half::uEltypeNoUnits
  c1::uCoeffType
  c2::uCoeffMType
end
@cache struct BAOABCache{uType,uEltypeNoUnits,rateNoiseType,uCoeffType, uCoeffMType,uTypeCombined} <: StochasticDynamicalEqMutableCache
  utmp::uType
  dumid::uType
  dutmp::uType
  dunoise::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c1::uCoeffType
  c2::uCoeffMType
  tmp::uTypeCombined
end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  if typeof(alg.gamma) <: AbstractMatrix
      c1 = exp(-alg.gamma*dt)
      c2 = cholesky(I - alg.scale_noise*c1*transpose(c1)).U# if scale_noise == false, c2 = 1
  else
      c1 = exp.(-alg.gamma*dt)
      c2 = sqrt.(1 .- alg.scale_noise*c1.^2)# if scale_noise == false, c2 = 1
  end
  BAOABConstantCache(k, uEltypeNoUnits(1//2),c1, c2)
end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dumid = zero(u.x[1])
  dutmp = zero(u.x[1])
  dunoise = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  gtmp = zero(noise_rate_prototype)
  noise = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)

  if typeof(alg.gamma) <: AbstractMatrix
      c1 = exp(-alg.gamma*dt)
      c2 = cholesky(I - alg.scale_noise*c1*transpose(c1)).U# if scale_noise == false, c2 = 1
  else
      c1 = exp.(-alg.gamma*dt)
      c2 = sqrt.(1 .- alg.scale_noise*c1.^2)# if scale_noise == false, c2 = 1
  end

  tmp = zero(u)

  BAOABCache(utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c1, c2, tmp)
end



mutable struct ABOBAConstantCache{uType,uEltypeNoUnits, uCoeffType, uCoeffMType} <: StochasticDynamicalEqConstantCache
  k::uType
  half::uEltypeNoUnits
  c₂::uCoeffType
  σ::uCoeffMType
end
@cache struct ABOBACache{uType,uEltypeNoUnits,rateNoiseType,uCoeffType, uCoeffMType,uTypeCombined} <: StochasticDynamicalEqMutableCache
  utmp::uType
  dumid::uType
  dutmp::uType
  dunoise::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::uCoeffType
  σ::uCoeffMType
  tmp::uTypeCombined
end

function alg_cache(alg::ABOBA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])

  if typeof(alg.gamma) <: AbstractMatrix
      c₂ = exp(-alg.gamma*dt)
      σ = cholesky(I - alg.scale_noise*c₂*transpose(c₂)).U
  else
      c₂ = exp.(-alg.gamma*dt)
      σ = sqrt.(1 .- alg.scale_noise*c₂.^2)
  end
   # if scale_noise == false, c2 = 1
  ABOBAConstantCache(k, uEltypeNoUnits(1//2), c₂, σ)
end

function alg_cache(alg::ABOBA,prob,u,ΔW,ΔZ,p,rate_prototype, noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  dumid = zero(u.x[1])
  dunoise = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  gtmp = zero(noise_rate_prototype)
  noise = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)

  if typeof(alg.gamma) <: AbstractMatrix
      c₂ = exp(-alg.gamma*dt)
      σ = cholesky(I - alg.scale_noise*c₂*transpose(c₂)).U
  else
      c₂ = exp.(-alg.gamma*dt)
      σ = sqrt.(1 .- alg.scale_noise*c₂.^2)
  end

  tmp = zero(u)

  ABOBACache(utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c₂, σ, tmp)
end




mutable struct OBABOConstantCache{uType,rateNoiseType, uEltypeNoUnits, uCoeffType, uCoeffMType} <: StochasticDynamicalEqConstantCache
  k::uType
  sig::rateNoiseType
  half::uEltypeNoUnits
  c₂::uCoeffType
  σ::uCoeffMType
end

@cache struct OBABOCache{uType,uEltypeNoUnits,rateNoiseType,uCoeffType, uCoeffMType,uTypeCombined} <: StochasticDynamicalEqMutableCache
  utmp::uType
  dumid::uType
  dutmp::uType
  dunoise::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::uCoeffType
  σ::uCoeffMType
  tmp::uTypeCombined
end

function alg_cache(alg::OBABO,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  sig = zero(noise_rate_prototype)
  half=uEltypeNoUnits(1//2)

  if typeof(alg.gamma) <: AbstractMatrix
      c₂ = exp(-alg.gamma*half*dt)
      σ = cholesky(I - alg.scale_noise*c₂*transpose(c₂)).U
  else
      c₂ = exp.(-alg.gamma*half*dt)
      σ = sqrt.(1 .- alg.scale_noise*c₂.^2)
  end
   # if scale_noise == false, c2 = 1
  OBABOConstantCache(k, sig, half, c₂, σ)
end

function alg_cache(alg::OBABO,prob,u,ΔW,ΔZ,p,rate_prototype, noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  dumid = zero(u.x[1])
  dunoise = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  gtmp = zero(noise_rate_prototype)
  noise = zero(rate_prototype.x[1])


  half = uEltypeNoUnits(1//2)

  if typeof(alg.gamma) <: AbstractMatrix
      c₂ = exp(-alg.gamma*half*dt)
      σ = cholesky(I - alg.scale_noise*c₂*transpose(c₂)).U
  else
      c₂ = exp.(-alg.gamma*half*dt)
      σ = sqrt.(1 .- alg.scale_noise*c₂.^2)
  end

  tmp = zero(u)

  OBABOCache(utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c₂, σ, tmp)
end
