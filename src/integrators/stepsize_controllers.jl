
function stepsize_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.q11 = DiffEqBase.value(DiffEqBase.fastpow(integrator.EEst,controller.beta1))
    integrator.q = DiffEqBase.value(integrator.q11/DiffEqBase.fastpow(integrator.qold,controller.beta2))
    @fastmath integrator.q = DiffEqBase.value(max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma)))
    nothing
end

@inline function step_accept_controller!(integrator::SDEIntegrator, alg)
    step_accept_controller!(integrator, integrator.opts.controller, alg)
    nothing
end

function step_accept_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.dtnew = DiffEqBase.value(integrator.dt/integrator.q) * oneunit(integrator.dt)
    nothing
end

function step_reject_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
    nothing
end


function stepsize_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    nothing
end

function step_accept_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    integrator.q = min(integrator.opts.gamma / integrator.EEst, integrator.opts.qmax)
    return integrator.dt * integrator.q
end

function step_reject_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    integrator.dt = integrator.opts.gamma * integrator.dt / integrator.EEst
    nothing
end


function stepsize_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    nothing
end

function step_accept_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    return integrator.EEst # use EEst for the τ
end

function step_reject_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    error("CaoTauLeaping should never reject steps")
end
