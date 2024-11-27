struct SDEDefaultInit <: DiffEqBase.DAEInitializationAlgorithm end

function DiffEqBase.initialize_dae!(integrator::SDEIntegrator, initializealg = integrator.initializealg)
    OrdinaryDiffEqCore._initialize_dae!(integrator, integrator.sol.prob, initializealg, Val(DiffEqBase.isinplace(integrator.sol.prob)))
end

function OrdinaryDiffEqCore._initialize_dae!(integrator::SDEIntegrator, prob, ::SDEDefaultInit, isinplace)
    if SciMLBase.has_initializeprob(prob.f)
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.OverrideInit(), isinplace)
    else
        OrdinaryDiffEqCore._initialize_dae!(integrator, prob, SciMLBase.CheckInit(), isinplace)
    end
end

