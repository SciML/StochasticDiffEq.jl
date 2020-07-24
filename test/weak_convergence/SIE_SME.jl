"""
 Tests for SIE and SME methods
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
    sol = solve(prob,alg;ensemblealg=ensemblealg,dt=dts[i],save_start=false,save_everystep=false,weak_timeseries_errors=false,weak_dense_errors=false,trajectories=Int(numtraj))
    println(i)
    push!(sols,sol)
  end
  return sols
end

function prob_func(prob, i, repeat)
    remake(prob,seed=seeds[i])
end

"""
 Test Scalar SDEs (oop)
"""

@info "Scalar oop noise"

# PC exercise 14.2.2
numtraj = Int(7e3)
u₀ = 0.1
f(u,p,t) = p[1]*u
g(u,p,t) = p[2]*u
dts = 1 .//2 .^(6:-1:3)
tspan = (0.0,1.0)
p = [3//2,1//100]

h1(z) = z
#analytical_sol(t) = E(f(X(t))) =

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f,g,u₀,tspan,p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, SIEA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

#using Plots; convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, " SIEA.png")
println("SIEA:", m)


_solutions = @time generate_weak_solutions(ensemble_prob, SMEA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SMEA:", m)


_solutions = @time generate_weak_solutions(ensemble_prob, SIEB(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SIEB:", m)

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"

u₀ = [0.1]
f1!(du,u,p,t) =  (du[1] = p[1]*u[1])
g1!(du,u,p,t) =  (du[1] = p[2]*u[1])
dts = 1 .//2 .^(6:-1:3)
tspan = (0.0,1.0)
p = [3//2,1//100]


prob = SDEProblem(f1!,g1!,u₀,tspan,p)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end][1]),false),
        prob_func = prob_func
        )


numtraj = Int(7e3)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)


_solutions = @time generate_weak_solutions(ensemble_prob, SIEA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SIEA:", m)


_solutions = @time generate_weak_solutions(ensemble_prob, SMEA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SMEA:", m)


_solutions = @time generate_weak_solutions(ensemble_prob, SIEB(), dts, numtraj, ensemblealg=EnsembleThreads())

errors =  [LinearAlgebra.norm(Statistics.mean(sol.u) .- u₀.*exp(1.0*(p[1]))) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SIEB:", m)

"""
 Test Diagonal noise SDEs (iip), SIAM Journal on Numerical Analysis, 47 (2009), pp. 1713–1738
"""

@info "Diagonal noise"

u₀ = [0.1,0.1]
function f3!(du,u,p,t)
  du[1] = 3//2*u[1]
  du[2] = 3//2*u[2]
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

numtraj = Int(5e4)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, SIEA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SIEA:", m)


_solutions = @time generate_weak_solutions(ensemble_prob, SMEA(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SMEA:", m)

_solutions = @time generate_weak_solutions(ensemble_prob, SIEB(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)-1//100*exp(301//100)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("SIEB:", m)
