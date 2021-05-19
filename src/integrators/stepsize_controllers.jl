
# TODO: HR
# function stepsize_controller!(integrator, alg)
#     integrator.q11 = DiffEqBase.value(DiffEqBase.fastpow(integrator.EEst,integrator.opts.beta1))
#     integrator.q = DiffEqBase.value(integrator.q11/DiffEqBase.fastpow(integrator.qold,integrator.opts.beta2))
#     @fastmath integrator.q = DiffEqBase.value(max(inv(integrator.opts.qmax),min(inv(integrator.opts.qmin),integrator.q/integrator.opts.gamma)))
# end

# function step_accept_controller!(integrator, alg)
#     integrator.dtnew = DiffEqBase.value(integrator.dt/integrator.q) * oneunit(integrator.dt)
# end

# function step_reject_controller!(integrator, alg)
#     integrator.dtnew = integrator.dt/min(inv(integrator.opts.qmin),integrator.q11/integrator.opts.gamma)
# end


function stepsize_controller!(integrator, alg::TauLeaping)
    nothing
end

function step_accept_controller!(integrator, alg::TauLeaping, dummy)
    integrator.q = min(integrator.opts.gamma / integrator.EEst, integrator.opts.qmax)
    return integrator.dt * integrator.q
end

function step_reject_controller!(integrator, alg::TauLeaping)
    integrator.dt = integrator.opts.gamma * integrator.dt / integrator.EEst
end


function stepsize_controller!(integrator, alg::CaoTauLeaping)
    nothing
end

function step_accept_controller!(integrator, alg::CaoTauLeaping, dummy)
    return integrator.EEst # use EEst for the τ
end

function step_reject_controller!(integrator, alg::CaoTauLeaping)
    error("CaoTauLeaping should never reject steps")
end
