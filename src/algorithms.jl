abstract StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm
abstract StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm

immutable EM <: StochasticDiffEqAlgorithm end
immutable RKMil <: StochasticDiffEqAlgorithm end
@with_kw immutable SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRA1()
end
@with_kw immutable SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRIW1()
end
immutable SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end
immutable SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

isadaptive(alg::StochasticDiffEqAlgorithm) = false
isadaptive(alg::StochasticDiffEqAdaptiveAlgorithm) = true
