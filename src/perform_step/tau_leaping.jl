@muladd function perform_step!(integrator, cache::TauLeapingConstantCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    tmp = c(uprev, p, t, P.dW, nothing)
    integrator.u = uprev .+ tmp

    if integrator.opts.adaptive
        if integrator.alg isa TauLeaping
            oldrate = P.cache.currate
            newrate = P.cache.rate(integrator.u, p, t + dt)
            EEstcache = @. abs(newrate - oldrate) /
                max(50integrator.opts.reltol * oldrate, integrator.rate_constants / integrator.dt)
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
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; tmp, newrate, EEstcache) = cache
    c(tmp, uprev, p, t, P.dW, nothing)
    @.. u = uprev + tmp

    if integrator.opts.adaptive
        if integrator.alg isa TauLeaping
            oldrate = P.cache.currate
            P.cache.rate(newrate, u, p, t + dt)
            @.. EEstcache = abs(newrate - oldrate) /
                max(50integrator.opts.reltol * oldrate, integrator.rate_constants / integrator.dt)
            integrator.EEst = maximum(EEstcache)
            if integrator.EEst <= 1
                P.cache.currate .= newrate
            end
        elseif integrator.alg isa CaoTauLeaping
            # Calculate τ as EEst
        end
    end
end

# ThetaTrapezoidalTauLeaping: Weak second order tau-leaping (Hu, Li, Min 2011)
# Algorithm:
# 1. Generate Poisson jumps k'_j = Poisson(a_j(X_n) * θ * τ)
# 2. Predictor: X' = X_n + ν * k'
# 3. Compute corrector rates: l_j = max(α₁ * a_j(X') - α₂ * a_j(X_n), 0)
# 4. Generate Poisson jumps k_j = Poisson(l_j * (1-θ) * τ)
# 5. Final update: X_{n+1} = X' + ν * k

@muladd function perform_step!(integrator, cache::ThetaTrapezoidalTauLeapingConstantCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; theta, alpha1, alpha2) = cache
    rng = P.rng

    # Get rates at current state
    rate_at_uprev = P.cache.rate(uprev, p, t)

    # Step 1: Generate Poisson random numbers for predictor step with rate θ*dt*a(X_n)
    predictor_counts = JumpProcesses.pois_rand.(Ref(rng), theta * dt .* rate_at_uprev)

    # Step 2: Compute predictor state X' = X_n + c(k')
    predictor = uprev .+ c(uprev, p, t, predictor_counts, nothing)

    # Step 3: Compute rates at predictor state
    rate_at_predictor = P.cache.rate(predictor, p, t + theta * dt)

    # Step 4: Compute corrector rates l_j = max(α₁*a_j(X') - α₂*a_j(X_n), 0)
    corrector_rate = @. max(alpha1 * rate_at_predictor - alpha2 * rate_at_uprev, zero(eltype(rate_at_uprev)))

    # Step 5: Generate Poisson random numbers for corrector step with rate (1-θ)*dt*l_j
    corrector_counts = JumpProcesses.pois_rand.(Ref(rng), (one(theta) - theta) * dt .* corrector_rate)

    # Step 6: Final update X_{n+1} = X' + c(k)
    tmp = c(predictor, p, t + theta * dt, corrector_counts, nothing)
    integrator.u = predictor .+ tmp

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.currate = P.cache.rate(integrator.u, p, t + dt)
    end
end

@muladd function perform_step!(integrator, cache::ThetaTrapezoidalTauLeapingCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; tmp, predictor, predictor_counts, corrector_counts,
        rate_at_predictor, rate_at_uprev, corrector_rate,
        theta, alpha1, alpha2) = cache
    rng = P.rng

    # Get rates at current state (in-place)
    P.cache.rate(rate_at_uprev, uprev, p, t)

    # Step 1: Generate Poisson random numbers for predictor step with rate θ*dt*a(X_n)
    @. predictor_counts = JumpProcesses.pois_rand(rng, theta * dt * rate_at_uprev)

    # Step 2: Compute predictor state X' = X_n + c(k')
    c(tmp, uprev, p, t, predictor_counts, nothing)
    @.. predictor = uprev + tmp

    # Step 3: Compute rates at predictor state (in-place)
    P.cache.rate(rate_at_predictor, predictor, p, t + theta * dt)

    # Step 4: Compute corrector rates l_j = max(α₁*a_j(X') - α₂*a_j(X_n), 0)
    @. corrector_rate = max(alpha1 * rate_at_predictor - alpha2 * rate_at_uprev,
        zero(eltype(rate_at_uprev)))

    # Step 5: Generate Poisson random numbers for corrector step with rate (1-θ)*dt*l_j
    @. corrector_counts = JumpProcesses.pois_rand(rng, (one(theta) - theta) * dt * corrector_rate)

    # Step 6: Final update X_{n+1} = X' + c(k)
    c(tmp, predictor, p, t + theta * dt, corrector_counts, nothing)
    @.. u = predictor + tmp

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.rate(P.cache.currate, u, p, t + dt)
    end
end
