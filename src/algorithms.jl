abstract StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm

immutable EM <: StochasticDiffEqAlgorithm end
immutable RKMil <: StochasticDiffEqAlgorithm end
immutable SRA <: StochasticDiffEqAlgorithm end
immutable SRI <: StochasticDiffEqAlgorithm end
immutable SRIW1 <: StochasticDiffEqAlgorithm end
immutable SRA1 <: StochasticDiffEqAlgorithm end
immutable SRAVectorized <: StochasticDiffEqAlgorithm end
immutable SRIVectorized <: StochasticDiffEqAlgorithm end
