abstract StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm
abstract StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm
abstract StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm

immutable EM <: StochasticDiffEqAlgorithm end
immutable RKMil <: StochasticDiffEqAlgorithm end
@with_kw immutable SRA{TabType,RSType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRA1()
  rswm::RSType = RSWM()
end
@with_kw immutable SRI{TabType,RSType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRIW1()
  rswm::RSType = RSWM()
end
@with_kw immutable SRIW1{RSType} <: StochasticDiffEqAdaptiveAlgorithm
  rswm::RSType = RSWM()
end
@with_kw immutable SRA1{RSType} <: StochasticDiffEqAdaptiveAlgorithm
  rswm::RSType = RSWM()
end

immutable StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

isadaptive(alg::StochasticDiffEqAlgorithm) = false
isadaptive(alg::StochasticDiffEqAdaptiveAlgorithm) = true
