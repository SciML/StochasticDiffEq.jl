"""
 Tests for  https://link.springer.com/article/10.1007/s10543-007-0130-3 with test problems as in the paper.
 RS1, RS2
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

numtraj = Int(1e6) # in the paper they use 1e9
u₀ = 0.1
f(u,p,t) = 1.5*u
g(u,p,t) = 0.1*u
dts = 1 .//2 .^(5:-1:1)
tspan = (0.0,1.0) # 2.0 in paper


h1(z) = z

seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

prob = SDEProblem(f,g,u₀,tspan)
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(sol[end]),false),
        prob_func = prob_func
        )
_solutions = @time generate_weak_solutions(ensemble_prob, RS1(), dts, numtraj, ensemblealg=EnsembleThreads())



errors = [LinearAlgebra.norm(Statistics.mean(sol.u) - u₀*exp(1.0*1.5) ) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

using Plots
convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "DRI1-CPU-final-"*string(numtraj)*".pdf")
println("RS1:", m)

solve(prob, RS2(), dt=dts[end])

numtraj = Int(1e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

_solutions = @time generate_weak_solutions(ensemble_prob, RS2(), dts, numtraj, ensemblealg=EnsembleThreads())

errors = [abs(LinearAlgebra.norm(Statistics.mean(sol.u)) - u₀*exp(1.0*1.5)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test -(m-2) < 0.3

println("RS2:", m)

"""
 Test Scalar SDEs (iip)
"""

@info "Scalar iip noise"


"""
 Test non-commutative noise SDEs (iip)
"""

@info "Non-commutative noise"




"""
 Test Diagonal noise SDEs (iip), SIAM Journal on Numerical Analysis, 47 (2009), pp. 1713–1738
"""

@info "Diagonal noise"
