abstract StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm
abstract StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm

immutable EM <: StochasticDiffEqAlgorithm end
immutable RKMil <: StochasticDiffEqAlgorithm end
immutable SRA <: StochasticDiffEqAdaptiveAlgorithm end
immutable SRI <: StochasticDiffEqAdaptiveAlgorithm end
immutable SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end
immutable SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

isadaptive(alg::StochasticDiffEqAlgorithm) = false
isadaptive(alg::StochasticDiffEqAdaptiveAlgorithm) = true
