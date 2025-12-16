# Regression test for issue #636: ImplicitRKMil error estimation bug with zero noise
# https://github.com/SciML/StochasticDiffEq.jl/issues/636
#
# Before fix: abstol=1e-9 gave 4.48 (4.2% error)
# After fix: abstol=1e-8 gives 4.303 (0.025% error)

using StochasticDiffEq, Test, Random

@testset "Issue 636 - ImplicitRKMil zero noise" begin
    # Simplified version of particle settling problem
    function nil!(dw, w, tp, t)
        for i in 1:6
            dw[i] = 0
        end
    end

    function finite_depth_u_crossing_nondim(x, y, z, t, p)
        (St, theta, W, K) = p
        k1 = [cos(theta / 2), sin(theta / 2)]
        k2 = [cos(theta / 2), -sin(theta / 2)]
        k1x = k1[1] * x + k1[2] * y
        k2x = k2[1] * x + k2[2] * y
        coshfac = cosh(z + K) / sinh(K)
        sinhfac = sinh(z + K) / sinh(K)
        u_x = -coshfac * (sin(k1x - t) * k1[1] + sin(k2x - t) * k2[1])
        u_y = -coshfac * (sin(k1x - t) * k1[2] + sin(k2x - t) * k2[2])
        u_z = sinhfac * (cos(k1x - t) + cos(k2x - t))
        return [u_x, u_y, u_z]
    end

    function drag_fd_3d_crossing_nondim!(dw, w, p, t)
        (St, theta, W, K) = p
        u = finite_depth_u_crossing_nondim(w[1], w[2], w[3], t, p)
        dw[1] = w[4]
        dw[2] = w[5]
        dw[3] = w[6]
        dw[4] = 1 / St * (u[1] - w[4])
        dw[5] = 1 / St * (u[2] - w[5])
        dw[6] = 1 / St * (u[3] - w[6] - W)
    end

    St = 0.1
    theta = 0
    tspan = (0.0, 10000.0)
    u0 = [0.0, 0.0, -2.0, 0.0, 0.0, 0.0]
    W = 0.0072
    K = 10
    p = (St, theta, W, K)

    condition_bottom(u, t, integrator) = u[3] + K
    cad = ContinuousCallback(condition_bottom, terminate!)

    prob_sde = SDEProblem(drag_fd_3d_crossing_nondim!, nil!, u0, tspan, p)

    # Test with medium tolerance - should give ~4.30 with the fix
    sol = solve(prob_sde, ImplicitRKMil(), callback = cad, abstol = 1e-8, reltol = 1e-8,
        maxiters = 1e9)

    d = sqrt(last(sol.u)[1]^2 + last(sol.u)[2]^2)

    # Expected value is ~4.302 (ODE reference)
    # With fix: should be within 1% (was 4.2% error before fix)
    expected = 4.302
    rel_error = abs(d - expected) / expected

    @test rel_error < 0.01  # Within 1% (before fix: 4.2% error)
    @test 4.0 < d < 4.5     # Sanity check on result
end
