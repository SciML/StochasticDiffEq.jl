qmax_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = 9//8
qmin_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = 1//5

isadaptive(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = false
isadaptive(alg::Union{StochasticDiffEqAdaptiveAlgorithm,StochasticDiffEqRODEAdaptiveAlgorithm}) = true
isadaptive(alg::Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}) = isadaptive(alg.algs[1])

alg_order(alg::EM) = 1//2
alg_order(alg::ImplicitEM) = 1//2
alg_order(alg::ImplicitEulerHeun) = 1//2
alg_order(alg::ImplicitRKMil) = 1//1
alg_order(alg::SplitEM) = 1//2
alg_order(alg::PCEuler) = 1//2
alg_order(alg::IIF1M) = 1//2
alg_order(alg::IIF2M) = 1//2
alg_order(alg::IIF1Mil) = 1//1
alg_order(alg::EulerHeun) = 1//2
alg_order(alg::RandomEM) = 1//2
alg_order(alg::RKMil) = 1//1
alg_order(alg::RKMilCommute) = 1//1
alg_order(alg::SRI) = alg.tableau.order
alg_order(alg::SRIW1) = 3//2
alg_order(alg::SRIW2) = 3//2
alg_order(alg::SOSRI) = 3//2
alg_order(alg::SOSRI2) = 3//2
alg_order(alg::SRA) = alg.tableau.order
alg_order(alg::SRA1) = 2//1
alg_order(alg::SRA2) = 2//1
alg_order(alg::SRA3) = 2//1
alg_order(alg::SOSRA) = 2//1
alg_order(alg::SOSRA2) = 2//1
alg_order(alg::RackKenCarp) = 2//1
alg_order(alg::Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}) = alg_order(alg.algs[1])

beta2_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = 2//(5alg_order(alg))
beta1_default(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm},beta2) = 7//(10alg_order(alg))

isdtchangeable(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = true

alg_interpretation(alg::StochasticDiffEqAlgorithm) = :Ito
alg_interpretation(alg::EulerHeun) = :Stratonovich
alg_interpretation(alg::RKMil{interpretation}) where {interpretation} = interpretation
alg_interpretation(alg::RKMilCommute{interpretation}) where {interpretation} = interpretation
alg_interpretation(alg::ImplicitRKMil{CS,AD,F,S,K,T,T2,Controller,interpretation}) where {CS,AD,F,S,K,T,T2,Controller,interpretation} = interpretation

alg_compatible(prob,alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = true

alg_compatible(prob,alg::StochasticDiffEqAlgorithm) = false
alg_compatible(prob,alg::SRI) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SRIW1) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SRIW2) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SOSRI) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SOSRI2) = is_diagonal_noise(prob)
alg_compatible(prob,alg::SRA) = true
alg_compatible(prob,alg::SRA1) = true
alg_compatible(prob,alg::SRA2) = true
alg_compatible(prob,alg::SRA3) = true
alg_compatible(prob,alg::SOSRA) = true
alg_compatible(prob,alg::SOSRA2) = true
alg_compatible(prob,alg::RackKenCarp) = true
alg_compatible(prob,alg::EM) = true
alg_compatible(prob,alg::EulerHeun) = true
alg_compatible(prob,alg::SplitEM) = true
alg_compatible(prob,alg::PCEuler) = true
alg_compatible(prob,alg::ImplicitEM) = true
alg_compatible(prob,alg::ImplicitEulerHeun) = true
alg_compatible(prob,alg::RKMil) = is_diagonal_noise(prob)
alg_compatible(prob,alg::ImplicitRKMil) = is_diagonal_noise(prob)
alg_compatible(prob,alg::RKMilCommute) = true # No good check for commutative noise
alg_compatible(prob,alg::IIF1M) = true
alg_compatible(prob,alg::IIF2M) = true
alg_compatible(prob,alg::Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}) = max((alg_compatible(prob,a) for a in alg.algs)...)

alg_needs_extra_process(alg::Union{StochasticDiffEqAlgorithm,StochasticDiffEqRODEAlgorithm}) = false
alg_needs_extra_process(alg::Union{StochasticDiffEqCompositeAlgorithm,StochasticDiffEqRODECompositeAlgorithm}) = max((alg_needs_extra_process(a) for a in alg.algs)...)
alg_needs_extra_process(alg::SRI) = true
alg_needs_extra_process(alg::SRIW1) = true
alg_needs_extra_process(alg::SRIW2) = true
alg_needs_extra_process(alg::SOSRI) = true
alg_needs_extra_process(alg::SOSRI2) = true
alg_needs_extra_process(alg::SRA) = true
alg_needs_extra_process(alg::SRA1) = true
alg_needs_extra_process(alg::SRA2) = true
alg_needs_extra_process(alg::SRA3) = true
alg_needs_extra_process(alg::SOSRA) = true
alg_needs_extra_process(alg::SOSRA2) = true
alg_needs_extra_process(alg::RackKenCarp) = true

alg_autodiff(alg::StochasticDiffEqNewtonAlgorithm{CS,AD,Controller}) where {CS,AD,Controller} = AD
alg_autodiff(alg::StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}) where {CS,AD,Controller} = AD
get_chunksize(alg::StochasticDiffEqNewtonAlgorithm{CS,AD,Controller}) where {CS,AD,Controller} = CS
get_chunksize(alg::StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}) where {CS,AD,Controller} = CS

alg_mass_matrix_compatible(alg::StochasticDiffEqAlgorithm) = false

function alg_mass_matrix_compatible(alg::StochasticDiffEqNewtonAlgorithm)
    if alg.symplectic
        return true
    elseif alg.theta == 1
        return true
    else
        error("Algorithm must be set as symplectic or theta=1 for mass matrices")
    end
end
