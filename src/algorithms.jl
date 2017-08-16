@compat abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
@compat abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
@compat abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

@compat abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
@compat abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
@compat abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

immutable EM <: StochasticDiffEqAlgorithm end
immutable SplitEM <: StochasticDiffEqAlgorithm end

immutable IIF1M{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF1M(;nlsolve=NLSOLVEJL_SETUP()) = IIF1M{typeof(nlsolve)}(nlsolve)

immutable IIF2M{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF2M(;nlsolve=NLSOLVEJL_SETUP()) = IIF2M{typeof(nlsolve)}(nlsolve)

immutable IIF1Mil{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF1Mil(;nlsolve=NLSOLVEJL_SETUP()) = IIF1Mil{typeof(nlsolve)}(nlsolve)

immutable EulerHeun <: StochasticDiffEqAlgorithm end
immutable RKMil{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMil(;interpretation=:Ito) = RKMil{interpretation}()

immutable RKMilCommute{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMilCommute(;interpretation=:Ito) = RKMilCommute{interpretation}()

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
