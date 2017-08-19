@compat abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
@compat abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
@compat abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

@compat abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
@compat abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
@compat abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

@compat abstract type StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller} <: StochasticDiffEqAdaptiveAlgorithm end
@compat abstract type StochasticDiffEqNewtonAlgorithm{CS,AD,Controller} <: StochasticDiffEqAlgorithm end

################################################################################

# Basics

struct EM <: StochasticDiffEqAlgorithm end
struct SplitEM <: StochasticDiffEqAlgorithm end
struct EulerHeun <: StochasticDiffEqAlgorithm end
struct RKMil{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMil(;interpretation=:Ito) = RKMil{interpretation}()

struct RKMilCommute{interpretation} <: StochasticDiffEqAlgorithm end
Base.@pure RKMilCommute(;interpretation=:Ito) = RKMilCommute{interpretation}()

################################################################################

# Rossler

@with_kw struct SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRA1()
end
@with_kw struct SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRIW1()
  error_terms = 4
end

struct SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end
struct SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

################################################################################

# IIF

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

################################################################################

# SDIRK

struct ImplicitEM{CS,AD,F,K,T,T2,Controller} <: StochasticDiffEqNewtonAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::Symbol
  κ::K
  tol::T
  extrapolant::Symbol
  min_newton_iter::Int
  max_newton_iter::Int
  new_jac_conv_bound::T2
end
Base.@pure ImplicitEM(;chunk_size=0,autodiff=true,diff_type=:central,
                          linsolve=DEFAULT_LINSOLVE,κ=nothing,tol=nothing,
                          extrapolant=:constant,min_newton_iter=1,
                          max_newton_iter=7,new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEM{chunk_size,autodiff,typeof(linsolve),
                          typeof(κ),typeof(tol),typeof(new_jac_conv_bound),controller}(
                          linsolve,diff_type,κ,tol,extrapolant,min_newton_iter,
                          max_newton_iter,new_jac_conv_bound)

################################################################################

# Etc.

struct StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

struct RandomEM <: StochasticDiffEqRODEAlgorithm end
