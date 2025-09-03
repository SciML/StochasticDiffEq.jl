@muladd function perform_step!(integrator, cache::TauLeapingConstantCache)
    @unpack t, dt, uprev, u, W, p, P, c = integrator
    # Handle case where P is Nothing (for pure discrete problems)
    dW = P === nothing ? nothing : P.dW
    tmp = c(uprev, p, t, dW, nothing)
    integrator.u = uprev .+ tmp

    if integrator.opts.adaptive && P !== nothing
        if integrator.alg isa TauLeaping
            oldrate = P.cache.currate
            newrate = P.cache.rate(integrator.u, p, t+dt)
            EEstcache = @. abs(newrate - oldrate) /
                           max(50integrator.opts.reltol*oldrate, integrator.rate_constants/integrator.dt)
            integrator.EEst = maximum(EEstcache)
            if integrator.EEst <= 1
                P.cache.currate = newrate
            end
        elseif integrator.alg isa CaoTauLeaping
            # Calculate τ as EEst
        end
    end
end

@muladd function perform_step!(integrator, cache::TauLeapingCache)
    @unpack t, dt, uprev, u, W, p, P, c = integrator
    @unpack tmp, newrate, EEstcache = cache
    # Handle case where P is Nothing (for pure discrete problems)
    dW = P === nothing ? nothing : P.dW
    c(tmp, uprev, p, t, dW, nothing)
    @.. u = uprev + tmp

    if integrator.opts.adaptive && P !== nothing
        if integrator.alg isa TauLeaping
            oldrate = P.cache.currate
            P.cache.rate(newrate, u, p, t+dt)
            @.. EEstcache = abs(newrate - oldrate) /
                            max(50integrator.opts.reltol*oldrate, integrator.rate_constants/integrator.dt)
            integrator.EEst = maximum(EEstcache)
            if integrator.EEst <= 1
                P.cache.currate .= newrate
            end
        elseif integrator.alg isa CaoTauLeaping
            # Calculate τ as EEst
        end
    end
end
