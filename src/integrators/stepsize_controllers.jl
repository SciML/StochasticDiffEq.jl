function stepsize_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    integrator.q11 = DiffEqBase.value(FastPower.fastpower(integrator.EEst, controller.beta1))
    q = DiffEqBase.value(integrator.q11 / FastPower.fastpower(integrator.qold, controller.beta2))
    return @fastmath DiffEqBase.value(
        max(
            inv(integrator.opts.qmax), min(inv(integrator.opts.qmin), q / integrator.opts.gamma)
        )
    )
end

# Dispatch to controller-specific method, passing q from stepsize_controller!
@inline function step_accept_controller!(integrator::SDEIntegrator, alg, q)
    return step_accept_controller!(integrator, integrator.opts.controller, alg, q)
end

function step_accept_controller!(integrator::SDEIntegrator, controller::PIController, alg, q)
    integrator.qold = max(integrator.EEst, integrator.opts.qoldinit)
    return DiffEqBase.value(integrator.dt / q) * oneunit(integrator.dt)
end

function step_reject_controller!(integrator::SDEIntegrator, controller::PIController, alg)
    return integrator.dt = integrator.dt / min(inv(integrator.opts.qmin), integrator.q11 / integrator.opts.gamma)
end

function stepsize_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    return nothing
end

function step_accept_controller!(integrator::SDEIntegrator, alg::TauLeaping, q)
    q_tau = min(integrator.opts.gamma / integrator.EEst, integrator.opts.qmax)
    return integrator.dt * q_tau
end

function step_reject_controller!(integrator::SDEIntegrator, alg::TauLeaping)
    return integrator.dt = integrator.opts.gamma * integrator.dt / integrator.EEst
end

function stepsize_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    return nothing
end

function step_accept_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping, q)
    return integrator.EEst # use EEst for the τ
end

function step_reject_controller!(integrator::SDEIntegrator, alg::CaoTauLeaping)
    error("CaoTauLeaping should never reject steps")
end
