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

# ThetaTrapezoidalTauLeaping: Implicit weak second order tau-leaping
# Based on Hu, Li, Min (2011) and Anderson, Mattingly (2011)
#
# Implicit equation:
#   X_{n+1} = X_n + ν*k + θ*dt*(drift(X_{n+1}) - drift(X_n))
# where k ~ Poisson(dt * a(X_n)) and drift(u) = ν*a(u) = c(u, p, t, rate(u, p, t), nothing)
#
# Rearranged for fixed-point iteration:
#   X_{n+1} = tmp + θ*z
# where:
#   tmp = X_n + ν*k - θ*dt*drift(X_n)  (explicit contribution)
#   z = dt*drift(X_{n+1})              (iteration variable)
#
# Uses NLFunctional-style fixed-point iteration (same convergence logic as OrdinaryDiffEqNonlinearSolve)

@muladd function perform_step!(integrator, cache::ThetaTrapezoidalTauLeapingConstantCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; theta, nlalg) = cache
    rng = P.rng

    # Step 1: Get rates at current state and generate Poisson counts
    rate_at_uprev = P.cache.rate(uprev, p, t)
    poisson_counts = JumpProcesses.pois_rand.(Ref(rng), dt .* rate_at_uprev)

    # Step 2: Compute explicit contributions
    # jump_contribution = ν*k
    jump_contribution = c(uprev, p, t, poisson_counts, nothing)
    # drift_at_uprev = ν*a(X_n)
    drift_at_uprev = c(uprev, p, t, rate_at_uprev, nothing)

    # Step 3: Set up iteration
    # tmp = X_n + ν*k - θ*dt*drift(X_n)
    tmp = uprev .+ jump_contribution .- theta .* dt .* drift_at_uprev

    # Initial guess: z = dt * drift(uprev) (linear extrapolation)
    z = dt .* drift_at_uprev

    # Step 4: Fixed-point iteration using NLFunctional parameters
    # Solve: z = dt * drift(tmp + θ*z)
    κ = nlalg.κ
    maxiters = nlalg.max_iter
    tstep = t + dt  # Evaluate drift at t + dt

    η = one(eltype(z))
    ndz = one(eltype(z))

    for iter in 1:maxiters
        # Compute current guess: ustep = tmp + θ*z
        ustep = tmp .+ theta .* z

        # Evaluate drift at ustep
        rate_at_ustep = P.cache.rate(ustep, p, tstep)
        drift_at_ustep = c(ustep, p, tstep, rate_at_ustep, nothing)

        # New iteration: ztmp = dt * drift(ustep)
        ztmp = dt .* drift_at_ustep

        # Compute residual
        dz = ztmp .- z

        # Compute norm of residuals (relative to state scale)
        atmp = calculate_residuals(
            dz, uprev, ustep, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        ndzprev = ndz
        ndz = integrator.opts.internalnorm(atmp, t)

        # Check for divergence and convergence
        if iter > 1
            θ_rate = ndz / ndzprev
            if θ_rate > 2  # Diverging
                integrator.force_stepfail = true
                return nothing
            end
            η = DiffEqBase.value(θ_rate / (1 - θ_rate))
        end

        z = ztmp

        # Check convergence
        if η * ndz < κ || ndz < 1e-12
            break
        end
    end

    # Step 5: Final update: X_{n+1} = tmp + θ*z
    integrator.u = tmp .+ theta .* z

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.currate = P.cache.rate(integrator.u, p, t + dt)
    end

    return nothing
end

@muladd function perform_step!(integrator, cache::ThetaTrapezoidalTauLeapingCache)
    (; t, dt, uprev, u, W, p, P, c) = integrator
    (; tmp, z, ztmp, k, drift_at_uprev, atmp, poisson_counts, rate_at_uprev,
        rate_tmp, theta, nlalg) = cache
    rng = P.rng

    # Step 1: Get rates at current state (in-place) and generate Poisson counts
    P.cache.rate(rate_at_uprev, uprev, p, t)
    @. poisson_counts = JumpProcesses.pois_rand(rng, dt * rate_at_uprev)

    # Step 2: Compute explicit contributions
    # Use u as scratch space for jump contribution: ν*k (will be overwritten at end)
    c(u, uprev, p, t, poisson_counts, nothing)  # u now holds ν*k temporarily
    # Compute drift at uprev: drift(X_n) = ν*a(X_n)
    c(drift_at_uprev, uprev, p, t, rate_at_uprev, nothing)

    # Step 3: Set up iteration
    # tmp = X_n + ν*k - θ*dt*drift(X_n)
    @.. tmp = uprev + u - theta * dt * drift_at_uprev

    # Initial guess: z = dt * drift(uprev)
    @.. z = dt * drift_at_uprev

    # Step 4: Fixed-point iteration using NLFunctional parameters
    # Solve: z = dt * drift(tmp + θ*z)
    κ = nlalg.κ
    maxiters = nlalg.max_iter
    tstep = t + dt

    η = one(eltype(z))
    ndz = one(eltype(z))

    for iter in 1:maxiters
        # Compute current guess: k = tmp + θ*z (reuse k as ustep)
        @.. k = tmp + theta * z

        # Evaluate rates at k
        P.cache.rate(rate_tmp, k, p, tstep)

        # Evaluate drift at k: ztmp = dt * drift(k)
        c(ztmp, k, p, tstep, rate_tmp, nothing)
        @.. ztmp = dt * ztmp

        # Compute residual: atmp = ztmp - z
        @.. atmp = ztmp - z

        # Compute norm of residuals
        ndzprev = ndz
        calculate_residuals!(atmp, atmp, uprev, k, integrator.opts.abstol,
            integrator.opts.reltol, integrator.opts.internalnorm, t)
        ndz = integrator.opts.internalnorm(atmp, t)

        # Check for divergence and convergence
        if iter > 1
            θ_rate = ndz / ndzprev
            if θ_rate > 2  # Diverging
                integrator.force_stepfail = true
                return nothing
            end
            η = DiffEqBase.value(θ_rate / (1 - θ_rate))
        end

        @.. z = ztmp

        # Check convergence
        if η * ndz < κ || ndz < 1e-12
            break
        end
    end

    # Step 5: Final update: X_{n+1} = tmp + θ*z
    @.. u = tmp + theta * z

    # Update currate for consistency
    if integrator.opts.adaptive
        P.cache.rate(P.cache.currate, u, p, t + dt)
    end

    return nothing
end
