function sde_determine_initdt(
        u0::uType, t::tType, tdir, dtmax, abstol, reltol,
        internalnorm, prob, order, integrator
    ) where {tType, uType}
    if integrator.P !== nothing
        # Don't have a good estimate with jumps
        return tdir * dtmax / 1.0e6
    end

    f = prob.f
    g = prob.g
    p = prob.p
    d₀ = internalnorm(
        ArrayInterface.aos_to_soa(u0) ./
            (abstol .+ internalnorm.(u0, t) .* reltol), t
    )
    dtmin = nextfloat(integrator.opts.dtmin)
    smalldt = tType(1 // 10^(6))
    if !isinplace(prob)
        f₀ = f(u0, p, t)
        if any(x -> any(isnan, x), f₀)
            @SciMLMessage(
                "First function call for f produced NaNs. Exiting.",
                integrator.opts.verbose, :init_NaN
            )
        end
        g₀ = 3g(u0, p, t)
        if any(x -> any(isnan, x), g₀)
            @SciMLMessage(
                "First function call for g produced NaNs. Exiting.",
                integrator.opts.verbose, :init_NaN
            )
        end
    else
        f₀ = zero(u0)
        if prob.noise_rate_prototype !== nothing
            g₀ = zero(prob.noise_rate_prototype)
        else
            g₀ = zero(u0)
        end
        f(f₀, u0, p, t)
        if any(x -> any(isnan, x), f₀)
            @SciMLMessage(
                "First function call for f produced NaNs. Exiting.",
                integrator.opts.verbose, :init_NaN
            )
        end
        g(g₀, u0, p, t)
        g₀ .*= 3
        if any(x -> any(isnan, x), g₀)
            @SciMLMessage(
                "First function call for g produced NaNs. Exiting.",
                integrator.opts.verbose, :init_NaN
            )
        end
    end

    d₁ = internalnorm(
        max.(internalnorm.(f₀ .+ g₀, t), internalnorm.(f₀ .- g₀, t)) ./
            (abstol .+ internalnorm.(u0, t) .* reltol),
        t
    )
    T0 = typeof(d₀)
    T1 = typeof(d₁)
    if d₀ < T0(1 // 10^(5)) || d₁ < T1(1 // 10^(5))
        dt₀ = smalldt
    else
        dt₀ = tType((d₀ / d₁) / 100)
    end
    dt₀ = min(dt₀, tdir * dtmax)
    u₁ = u0 .+ tdir .* dt₀ .* f₀
    if !isinplace(prob)
        f₁ = f(u₁, p, t + tdir * dt₀)
        g₁ = 3g(u₁, p, t + tdir * dt₀)
    else
        f₁ = zero(u0)
        if prob.noise_rate_prototype !== nothing
            g₁ = zero(prob.noise_rate_prototype)
        else
            g₁ = zero(u0)
        end
        f(f₁, u0, p, t)
        g(g₁, u0, p, t)
        g₁ .*= 3
    end
    ΔgMax = max.(internalnorm.(g₀ .- g₁, t), internalnorm.(g₀ .+ g₁, t))
    d₂ = internalnorm(
        max.(internalnorm.(f₁ .- f₀ .+ ΔgMax, t), internalnorm.(f₁ .- f₀ .- ΔgMax, t)) ./
            (abstol .+ internalnorm.(u0, t) .* reltol),
        t
    ) ./ dt₀
    if max(d₁, d₂) <= T1(1 // Int64(10)^(15))
        dt₁ = max(smalldt, dt₀ * 1 // 10^(3))
    else
        dt₁ = tType(10^(-(2 + log10(max(d₁, d₂))) / (order + 1 // 2)))
    end
    return dt = tdir * max(dtmin, min(100dt₀, dt₁, tdir * dtmax))
end
