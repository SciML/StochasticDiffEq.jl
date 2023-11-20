using StochasticDiffEq, DiffEqNoiseProcess, Test, DiffEqDevTools, Random
Random.seed!(1)

f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
g(u,p,t) = 1 .+zero(u)
γ = 1

@testset "Dynamical convergence" begin
    u0 = 0
    v0 = 1

    ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
    prob1 = DynamicalSDEProblem(ff_harmonic,v0,u0,(0.0,1.5))

    dts = (1/2) .^ (6:-1:3)

    # Can't use NoiseGrid as noise is not generated with the correct size in convergence.jl. We require noise with shape of v.
    sim1  = analyticless_test_convergence(dts,prob1,BAOAB(gamma=γ),(1/2)^10;trajectories=Int(1e5),use_noise_grid=false)
    display(sim1.𝒪est)
    @test abs(sim1.𝒪est[:weak_final]-2) < 0.5

    sim1  = analyticless_test_convergence(dts,prob1,ABOBA(gamma=γ),(1/2)^10;trajectories=Int(1e5),use_noise_grid=false)
    display(sim1.𝒪est)
    @test abs(sim1.𝒪est[:weak_final]-2) < 0.5


    sim1  = analyticless_test_convergence(dts,prob1,OBABO(gamma=γ),(1/2)^10;trajectories=Int(1e5),use_noise_grid=false)
    display(sim1.𝒪est)
    @test abs(sim1.𝒪est[:weak_final]-2) < 0.8
end
