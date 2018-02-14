struct PCEulerConstantCache <: StochasticDiffEqConstantCache end

alg_cache(alg::PCEuler,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = PCEulerConstantCache()

struct PCEulerCache{uType,rateType,rateNoiseType,rateNoiseCollectionType} <: StochasticDiffEqMutableCache
  utmp::uType
  ftmp::rateType
  gtmp::rateNoiseType
  gdWtmp::rateNoiseCollectionType
  bbprimetmp::rateType
end

function alg_cache(alg::PCEuler,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{true}})
  utmp = similar(u); ftmp = zeros(rate_prototype);
  gtmp = zeros(noise_rate_prototype)
  bbprimetmp = similar(ftmp)
  if is_diagonal_noise(prob)
    gdWtmp = gtmp
  else
    gdWtmp = zeros(rate_prototype)
  end
  PCEulerCache(utmp,ftmp,gtmp,gdWtmp,bbprimetmp)
end
