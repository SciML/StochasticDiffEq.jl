qmax_default(alg::StochasticDiffEqAlgorithm) = 9//8
qmin_default(alg::StochasticDiffEqAlgorithm) = 1//5

alg_order(alg::EM) = 1//2
alg_order(alg::RKMil) = 1//1
alg_order(alg::SRI) = alg.tableau.order
alg_order(alg::SRIW1) = 3//2
alg_order(alg::SRA) = alg.tableau.order
alg_order(alg::SRA1) = 2//1

beta2_default(alg::StochasticDiffEqAlgorithm) = 2//(5alg_order(alg))
beta1_default(alg::StochasticDiffEqAlgorithm,beta2) = 7//(10alg_order(alg))
