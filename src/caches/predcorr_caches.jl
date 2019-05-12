struct PCEulerConstantCache <: StochasticDiffEqConstantCache end

alg_cache(alg::PCEuler,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}}) = PCEulerConstantCache()

@cache struct PCEulerCache{uType,rateType,rateNoiseType,rateNoiseCollectionType} <: StochasticDiffEqMutableCache
  utmp::uType
  ftmp::rateType
  gtmp::rateNoiseType
  gdWtmp::rateNoiseCollectionType
  bbprimetmp::rateType
end

function alg_cache(alg::PCEuler,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  utmp = zero(u); ftmp = zero(rate_prototype);
  gtmp = zero(noise_rate_prototype)
  bbprimetmp = zero(ftmp)
  if is_diagonal_noise(prob)
    gdWtmp = gtmp
  else
    gdWtmp = zero(rate_prototype)
  end
  PCEulerCache(utmp,ftmp,gtmp,gdWtmp,bbprimetmp)
end
