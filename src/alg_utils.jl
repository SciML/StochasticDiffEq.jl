qmax_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = 9//8
qmin_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = 1//5

alg_order(alg::EM) = 1//2
alg_order(alg::EulerHeun) = 1//2
alg_order(alg::RandomEM) = 1//2
alg_order(alg::RKMil) = 1//1
alg_order(alg::SRI) = alg.tableau.order
alg_order(alg::SRIW1) = 3//2
alg_order(alg::SRA) = alg.tableau.order
alg_order(alg::SRA1) = 2//1
alg_order(alg::Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}) = alg_order(alg.algs[1])

beta2_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = 2//(5alg_order(alg))
beta1_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm},beta2) = 7//(10alg_order(alg))

isdtchangeable(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = true

alg_interpretation(alg::StochasticDiffEqAlgorithm) = :Ito
alg_interpretation(alg::EulerHeun) = :Stratonovich
alg_interpretation{RSWM,interpretation}(alg::RKMil{RSWM,interpretation}) = interpretation

alg_compatible(prob,alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = true

alg_compatible(prob,alg::SRI) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SRIW1) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SRA) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SRA1) = is_diagonal_noise(prob)
