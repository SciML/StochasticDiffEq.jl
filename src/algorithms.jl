@compat abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
@compat abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
@compat abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

@compat abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
@compat abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
@compat abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

struct EM <: StochasticDiffEqAlgorithm end
struct SplitEM <: StochasticDiffEqAlgorithm end

struct IIF1M{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF1M(;nlsolve=NLSOLVEJL_SETUP()) = IIF1M{typeof(nlsolve)}(nlsolve)

struct IIF2M{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF2M(;nlsolve=NLSOLVEJL_SETUP()) = IIF2M{typeof(nlsolve)}(nlsolve)

struct IIF1Mil{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
Base.@pure IIF1Mil(;nlsolve=NLSOLVEJL_SETUP()) = IIF1Mil{typeof(nlsolve)}(nlsolve)

struct EulerHeun <: StochasticDiffEqAlgorithm end
struct RKMil{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMil(;interpretation=:Ito) = RKMil{interpretation}()

struct RKMilCommute{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMilCommute(;interpretation=:Ito) = RKMilCommute{interpretation}()

@with_kw struct SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRA1()
end
@with_kw struct SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRIW1()
  error_terms = 4
end

struct SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end
struct SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

struct StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

struct RandomEM <: StochasticDiffEqRODEAlgorithm end
