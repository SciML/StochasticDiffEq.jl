@compat abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
@compat abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
@compat abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

@compat abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
@compat abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
@compat abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

@with_kw immutable EM{RSType} <: StochasticDiffEqAlgorithm
  rswm::RSType = RSWM(adaptivealg=:RSwM1)
end
@with_kw immutable EulerHeun{RSType} <: StochasticDiffEqAlgorithm
  rswm::RSType = RSWM(adaptivealg=:RSwM1)
end

immutable RKMil{RSType,interpretation} <: StochasticDiffEqAlgorithm
  rswm::RSType
end
Base.@pure function RKMil(;rswm=RSWM(adaptivealg=:RSwM1),interpretation=:Ito)
  RKMil{typeof(rswm),interpretation}(rswm)
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

Base.@pure StochasticCompositeAlgorithm(algs,choice_function;rswm = isadaptive(algs[1]) ? RSWM() : RSWM(adaptivealg=:RSwM1)) = StochasticCompositeAlgorithm(algs,choice_function,rswm)
isadaptive(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = false
isadaptive(alg::Union{StochasticDiffEqAdaptiveAlgorithm,StochasticDiffEqRODEAdaptiveAlgorithm}) = true
isadaptive(alg::Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}) = isadaptive(alg.algs[1])

@with_kw immutable RandomEM{RSType} <: StochasticDiffEqRODEAlgorithm
  rswm::RSType = RSWM(adaptivealg=:RSwM1)
end
