@compat abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
@compat abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
@compat abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

@compat abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
@compat abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
@compat abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

immutable EM <: StochasticDiffEqAlgorithm end
immutable EulerHeun <: StochasticDiffEqAlgorithm end
immutable RKMil{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMil(;interpretation=:Ito) = RKMil{interpretation}()

@with_kw immutable SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRA1()
end
@with_kw immutable SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRIW1()
  error_terms = 4
end

immutable SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end
immutable SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

immutable StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

immutable RandomEM <: StochasticDiffEqRODEAlgorithm end
