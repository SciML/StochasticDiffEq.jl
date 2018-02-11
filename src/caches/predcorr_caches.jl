struct PCEulerConstantCache <: StochasticDiffEqConstantCache end

alg_cache(alg::PCEuler,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,uEltypeNoUnits,uBottomEltype,tTypeNoUnits,uprev,f,t,::Type{Val{false}}) = PCEulerConstantCache()
