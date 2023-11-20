using StochasticDiffEq, DiffEqNoiseProcess, Test, DiffEqDevTools, Random
Random.seed!(1)

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
g(u,p,t) = 1 .+zero(u)
γ = 1


@testset "Vector u" begin

    u0 = zeros(2)
    v0 = ones(2)

    f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
    f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
    g_iip(du,u,p,t) = du .= g(u,p,t)

    ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
    prob1 = DynamicalSDEProblem(ff_harmonic,v0,u0,(0.0,5.0))
    sol1 = solve(prob1,BAOAB(gamma=[γ,γ]);dt=1/10,save_noise=true)

    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,BAOAB(gamma=[γ,γ]);dt=1/10)

    @test sol1[:] ≈ sol2[:]

    dts = (1/2) .^ (8:-1:4)

    # # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    # sim1  = analyticless_test_convergence(dts,prob1,BAOAB(gamma=γ),(1/2)^10;trajectories=Int(1e2),use_noise_grid=false)
    # @test abs(sim1.𝒪est[:weak_final]-1.5) < 0.3


    sol1 = solve(prob1,ABOBA(gamma=[γ,γ]);dt=1/10,save_noise=true)
    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,ABOBA(gamma=[γ,γ]);dt=1/10)

    @test sol1[:] ≈ sol2[:]

    # # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    # sim1  = analyticless_test_convergence(dts,prob1,ABOBA(gamma=γ),(1/2)^10;trajectories=Int(1e4),use_noise_grid=false)
    # @test abs(sim1.𝒪est[:weak_final]-2) < 0.3


    sol1 = solve(prob1,OBABO(gamma=[γ,γ]);dt=1/10,save_noise=true)
    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,OBABO(gamma=[γ,γ]);dt=1/10)

    @test sol1[:] ≈ sol2[:]

    # # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    # sim1  = analyticless_test_convergence(dts,prob1,OBABO(gamma=γ),(1/2)^10;trajectories=Int(1e4),use_noise_grid=false)
    # @test abs(sim1.𝒪est[:weak_final]-1.5) < 0.3
end



@testset "Matricial gamma" begin

    u0 = zeros(2)
    v0 = ones(2)

    gamma_mat = [γ -0.1*γ;0.5*γ γ]

    f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
    f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
    g_iip(du,u,p,t) = du .= g(u,p,t)

    ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
    prob1 = DynamicalSDEProblem(ff_harmonic,v0,u0,(0.0,5.0))
    sol1 = solve(prob1,BAOAB(gamma=gamma_mat);dt=1/10,save_noise=true)

    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,BAOAB(gamma=gamma_mat);dt=1/10)

    @test sol1[:] ≈ sol2[:]


    sol1 = solve(prob1,ABOBA(gamma=gamma_mat);dt=1/10,save_noise=true)
    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,ABOBA(gamma=gamma_mat);dt=1/10)

    @test sol1[:] ≈ sol2[:]


    sol1 = solve(prob1,OBABO(gamma=gamma_mat);dt=1/10,save_noise=true)
    prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0); noise=NoiseWrapper(sol1.W))
    sol2 = solve(prob2,OBABO(gamma=gamma_mat);dt=1/10)

    @test sol1[:] ≈ sol2[:]
end

@testset "IIP and OOP consistency" begin
    f1_harmonic_iip(dv, v, u, p, t) = dv .= f1_harmonic(v, u, p, t)
    f2_harmonic_iip(du, v, u, p, t) = du .= f2_harmonic(v, u, p, t)
    g_iip(du, u, p, t) = du .= g(u, p, t)

    u0 = zeros(1)
    v0 = ones(1)

    dt = 0.01
    T = 1.5
    t = 0:dt:T

    brownian_values = cumsum([[zeros(length(u0))];
        [sqrt(dt) * randn(length(u0)) for i in 1:(length(t)-1)]])
    brownian_values2 = cumsum([[zeros(length(u0))];
        [sqrt(dt) * randn(length(u0)) for i in 1:(length(t)-1)]])
    W = NoiseGrid(t, brownian_values, brownian_values2)

    ff_harmonic = DynamicalSDEFunction(f1_harmonic, f2_harmonic, g)
    prob = DynamicalSDEProblem(ff_harmonic, v0, u0, (0.0, T), noise=W)

    ff_harmonic_iip = DynamicalSDEFunction(f1_harmonic_iip, f2_harmonic_iip, g_iip)
    prob_iip = DynamicalSDEProblem(ff_harmonic_iip, v0, u0, (0.0, T), noise=W)

    sol = solve(prob, BAOAB(gamma=γ), dt=dt)
    sol_iip = solve(prob_iip, BAOAB(gamma=γ), dt=dt)

    @test sol.u ≈ sol_iip.u

    sol = solve(prob, ABOBA(gamma=γ), dt=dt)
    sol_iip = solve(prob_iip, ABOBA(gamma=γ), dt=dt)

    @test sol.u ≈ sol_iip.u

    sol = solve(prob, OBABO(gamma=γ), dt=dt)
    sol_iip = solve(prob_iip, OBABO(gamma=γ), dt=dt)

    @test sol.u ≈ sol_iip.u
end
