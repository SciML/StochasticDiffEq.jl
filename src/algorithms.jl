abstract StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm
abstract StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm
abstract StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm

@with_kw immutable EM{RSType} <: StochasticDiffEqAlgorithm
  rswm::RSType = RSWM(adaptivealg=:RSwM1)
end
@with_kw immutable RKMil{RSType} <: StochasticDiffEqAlgorithm
  rswm::RSType = RSWM(adaptivealg=:RSwM1)
end
@with_kw immutable SRA{TabType,RSType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRA1()
  rswm::RSType = RSWM()
end
@with_kw immutable SRI{TabType,RSType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType = constructSRIW1()
  rswm::RSType = RSWM()
  error_terms = 4
end
@with_kw immutable SRIW1{RSType} <: StochasticDiffEqAdaptiveAlgorithm
  rswm::RSType = RSWM()
end
@with_kw immutable SRA1{RSType} <: StochasticDiffEqAdaptiveAlgorithm
  rswm::RSType = RSWM()
end

immutable StochasticCompositeAlgorithm{T,F,RSType} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
  rswm::RSType
end

Base.@pure StochasticCompositeAlgorithm(algs,choice_function;rswm=RSWM()) = StochasticCompositeAlgorithm(algs,choice_function,rswm)
isadaptive(alg::StochasticDiffEqAlgorithm) = false
isadaptive(alg::StochasticDiffEqAdaptiveAlgorithm) = true
isadaptive(alg::StochasticDiffEqCompositeAlgorithm) = isadaptive(alg.algs[1])
