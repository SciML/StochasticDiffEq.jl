"""
 Tests for  https://link.springer.com/article/10.1007/s10543-007-0130-3 with test problems as in the paper.
 RS1, RS2
 and for https://www.sciencedirect.com/science/article/pii/S0377042706003906
 NON
"""


import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
#using DiffEqGPU

function generate_weak_solutions(prob, alg, dts, numtraj; ensemblealg=EnsembleThreads())
  sols = []
  for i in 1:length(dts)
    sol = solve(prob,alg;ensemblealg=ensemblealg,dt=dts[i],adaptive=false,save_start=false,save_everystep=false,weak_timeseries_errors=false,weak_dense_errors=false,trajectories=Int(numtraj))
    println(i)
    push!(sols,sol)
  end
  return sols
end

function strat_weak_errors(prob, alg, dts, numtraj; ensemblealg=EnsembleThreads())
  errors = Array{typeof(prob.prob.u0)}(undef, length((dts)), 2)
  for i in 1:length(dts)
    sol = solve(prob, alg, dt=dts[i], #save_start=false,save_everystep=true,weak_timeseries_errors=false,weak_dense_errors=false,
    trajectories=numtraj, adaptive=false)

    computed_exp = Statistics.mean([_sol.u[end] for _sol in sol])
    true_exp1 = prob.prob.u0*exp(1.0*(prob.prob.p[1]+0.5*prob.prob.p[2]^2))
    #true_exp1b = prob.u0*exp(1.0*(prob.p[1][1]))
    true_exp2 = Statistics.mean([_sol.u_analytic[end] for _sol in sol])
    println(i)
    errors[i,1] = LinearAlgebra.norm(computed_exp-true_exp1)
    errors[i,2] = LinearAlgebra.norm(computed_exp-true_exp2)
  end
  return errors
end

function prob_func(prob, i, repeat)
    remake(prob,seed=seeds[i])
end

"""
 Test Scalar SDEs (oop)
"""

@info "Scalar oop noise"

numtraj = Int(1e5) # in the paper they use 1e9
u₀ = 0.1
p = [1.5, 0.1]
f(u,p,t) = p[1]*u
g(u,p,t) = p[2]*u
dts = 1 .//2 .^(5:-1:1)
tspan = (0.0,1.0)

h1(z) = z

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)



prob = SDEProblem(f,g,u₀,tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())



errors = [LinearAlgebra.norm(Statistics.mean(sol.u) - u₀*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3
println("RS1:", m)

#using Plots; convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "RS1-final-"*string(numtraj)*".pdf")

# f_strat(u0,p,t,W) = @.(u0*exp(p[1]*t+p[2]*W))
# seed = 100
# Random.seed!(seed)
# seeds = rand(UInt, numtraj)
# prob_test = SDEProblem(SDEFunction(f,g,analytic=f_strat),g,u₀,tspan,p)
# ensemble_test = EnsembleProblem(prob_test;
#         #output_func = (sol,i) -> (h1(sol[end]),false),
#         prob_func = prob_func
#         )
# errors2 = @time strat_weak_errors(ensemble_test, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())
#
# @test m ≈ log(errors2[end,1]/errors2[1,1])/log(dts[end]/dts[1])
# @test -(log(errors2[end,2]/errors2[1,2])/log(dts[end]/dts[1])-2) < 0.3


seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) - u₀*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS2:", m)

dts = 1 .//2 .^(5:-1:0)
numtraj = Int(1e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f,g,u₀,tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )


_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("NON:", m)

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

f!(du,u,p,t) = du[1] = p[1]*u[1]
g!(du,u,p,t) = du[1] = p[2]*u[1]

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f!,g!,[u₀],tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts,
  numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3
println("RS1:", m)


seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f!,g!,[u₀],tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts,
  numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS1:", m)

dts = 1 .//2 .^(5:-1:0)
numtraj = Int(1e5)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f!,g!,[u₀],tspan, p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts,
  numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]+0.5*p[2]^2))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])


"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"

u₀ = [0.1,0.1]
function f2!(du,u,p,t)
  du[1] = 5//4*u[2]-5//4*u[1]
  du[2] = 1//4*u[1]-1//4*u[2]
end
function g2!(du,u,p,t)
  du[1,1] = sqrt(3)/2*(u[1]-u[2])
  du[1,2] = 1//2*(u[1]+u[2])
  #du[2,1] = 0
  du[2,2] = u[1]
end
dts = 1 .//2 .^(4:-1:1)
tspan = (0.0,1.0)

h2(z) = z^2 # E(x_i) = 1/10 exp(1/2t) or E(x_1* x_2) = 1/100 exp(2t)

prob = SDEProblem(f2!,g2!,u₀,tspan,noise_rate_prototype=zeros(2,2))
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h2(sol[end][1]),false),
        prob_func = prob_func
        )

numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS1:", m)


numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS2:", m)

numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(2)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("NON:", m)

"""
 Test Diagonal noise SDEs (iip)
"""

@info "Diagonal noise"


u₀ = [0.1,0.1]
function f3!(du,u,p,t)
  du[1] = 299//200*u[1]
  du[2] = 299//200*u[2]
end
function g3!(du,u,p,t)
  du[1] = 1//10*u[1]
  du[2] = 1//10*u[2]
end
dts = 1 .//2 .^(5:-1:1)
tspan = (0.0,1.0)

h3(z) = z^2 # == 1//10**exp(3//2*t) if h3(z) = z and  == 1//100**exp(301//100*t) if h3(z) = z^2 )

prob = SDEProblem(f3!,g3!,u₀,tspan)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h3(sol[end][1]),false),
        prob_func = prob_func
        )

numtraj = Int(4e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS1:", m)

numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS2:", m)


numtraj = Int(5e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, NON(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("NON:", m)

using Plots; convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
