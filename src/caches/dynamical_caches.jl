abstract type StochasticDynamicalEqConstantCache <: StochasticDiffEqConstantCache end # Pourquoi faire ça, Si c'est pour avoir une seul function de check dans initialize!
abstract type StochasticDynamicalEqMutableCache <: StochasticDiffEqMutableCache end


mutable struct BAOABConstantCache{uType,uEltypeNoUnits,uCoeffType, uCoeffMType} <: StochasticDynamicalEqConstantCache
  k::uType
  half::uEltypeNoUnits
  c2::uCoeffType
  σ::uCoeffMType
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
  c2::uCoeffType
  σ::uCoeffMType
  tmp::uTypeCombined
end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  if typeof(alg.gamma) <: AbstractMatrix
      c2 = exp(-alg.gamma*dt)
      σ = cholesky(I - alg.scale_noise*c2*transpose(c2)).U# if scale_noise == false, c2 = 1
  else
      c2 = exp.(-alg.gamma*dt)
      σ = sqrt.(1 .- alg.scale_noise*c2.^2)# if scale_noise == false, c2 = 1
  end
  BAOABConstantCache(k, uEltypeNoUnits(1//2),c2, c2)
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
      c2 = exp(-alg.gamma*dt)
      σ = cholesky(I - alg.scale_noise*c2*transpose(c2)).U# if scale_noise == false, c2 = 1
  else
      c2 = exp.(-alg.gamma*dt)
      σ = sqrt.(1 .- alg.scale_noise*c2.^2)# if scale_noise == false, c2 = 1
  end

  tmp = zero(u)

  BAOABCache(utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c2, σ, tmp)
end



mutable struct ABOBAConstantCache{uType,uEltypeNoUnits, uCoeffType, uCoeffMType} <: StochasticDynamicalEqConstantCache
  k::uType
  half::uEltypeNoUnits
  c2::uCoeffType
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
  c2::uCoeffType
  σ::uCoeffMType
  tmp::uTypeCombined
end

function alg_cache(alg::ABOBA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])

  if typeof(alg.gamma) <: AbstractMatrix
      c2 = exp(-alg.gamma*dt)
      σ = cholesky(I - alg.scale_noise*c2*transpose(c2)).U
  else
      c2 = exp.(-alg.gamma*dt)
      σ = sqrt.(1 .- alg.scale_noise*c2.^2)
  end
   # if scale_noise == false, c2 = 1
  ABOBAConstantCache(k, uEltypeNoUnits(1//2), c2, σ)
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
      c2 = exp(-alg.gamma*dt)
      σ = cholesky(I - alg.scale_noise*c2*transpose(c2)).U
  else
      c2 = exp.(-alg.gamma*dt)
      σ = sqrt.(1 .- alg.scale_noise*c2.^2)
  end

  tmp = zero(u)

  ABOBACache(utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c2, σ, tmp)
end




mutable struct OBABOConstantCache{uType,rateNoiseType, uEltypeNoUnits, uCoeffType, uCoeffMType} <: StochasticDynamicalEqConstantCache
  k::uType
  gt::rateNoiseType
  half::uEltypeNoUnits
  c2::uCoeffType
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
  c2::uCoeffType
  σ::uCoeffMType
  tmp::uTypeCombined
end

function alg_cache(alg::OBABO,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  gt = zero(noise_rate_prototype)
  half=uEltypeNoUnits(1//2)

  if typeof(alg.gamma) <: AbstractMatrix
      c2 = exp(-alg.gamma*half*dt)
      σ = cholesky(I - alg.scale_noise*c2*transpose(c2)).U
  else
      c2 = exp.(-alg.gamma*half*dt)
      σ = sqrt.(1 .- alg.scale_noise*c2.^2)
  end
   # if scale_noise == false, c2 = 1
  OBABOConstantCache(k, gt, half, c2, σ)
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
      c2 = exp(-alg.gamma*half*dt)
      σ = cholesky(I - alg.scale_noise*c2*transpose(c2)).U
  else
      c2 = exp.(-alg.gamma*half*dt)
      σ = sqrt.(1 .- alg.scale_noise*c2.^2)
  end

  tmp = zero(u)

  OBABOCache(utmp, dumid, dutmp, dunoise, k, gtmp, noise, half, c2, σ, tmp)
end
