struct SDEDefaultInit <: DiffEqBase.DAEInitializationAlgorithm end

function DiffEqBase.initialize_dae!(integrator::AbstractSDEIntegrator, initializealg = integrator.initializealg)
    OrdinaryDiffEqCore._initialize_dae!(integrator, integrator.sol.prob, initializealg, Val(DiffEqBase.isinplace(integrator.sol.prob)))
end

function OrdinaryDiffEqCore._initialize_dae!(integrator::AbstractSDEIntegrator, prob, ::SDEDefaultInit, isinplace)
    if SciMLBase.has_initializeprob(prob.f)
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.OverrideInit(), isinplace)
    elseif SciMLBase.__has_mass_matrix(prob.f)
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.CheckInit(), isinplace)
    end
end

