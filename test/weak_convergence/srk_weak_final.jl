"""
 oop tests for https://arxiv.org/abs/1303.5103 with test problems as in the paper.
 DRI1 (and RI1)
"""


import Statistics # for mean values of trajectories
import LinearAlgebra # for the normn
using StochasticDiffEq
using Test
using Random
Random.seed!(100)
#using DiffEqGPU
#using Plots

"""
 Test Scalar SDEs
"""

numtraj = 5e5 # in the paper they use 1e9
u₀ = 0.0
f(u,p,t) = 1//2*u+sqrt(u^2+1)
g(u,p,t) = sqrt(u^2+1)
dts = 1 .//2 .^(4:-1:1)
tspan = (0.0,2.0)
prob = SDEProblem(f,g,u₀,tspan)

h1(z) = z^3-6*z^2+8*z

#analytical_sol(t) = E(f(X(t))) = E(h1(arsinh(X(t))) = t^3-3*t^2+2*t
#analytical_sol(2) = 0 and analytical_sol(1)=0
ensemble_prob = EnsembleProblem(prob;
        output_func = (sol,i) -> (h1(asinh(sol[end])),false)
        )
N = length(dts)
_solutions = @time [solve(ensemble_prob,
        DRI1();
        #ensemblealg = EnsembleGPUArray(); # EnsembleGPUArray or EnsembleDistributed()
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.3

#convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "DRI1-CPU-final-"*string(numtraj)*".pdf")
#println(m)



_solutions = @time [solve(ensemble_prob,
        RI1();
        #ensemblealg = EnsembleGPUArray(); # EnsembleGPUArray or EnsembleDistributed()
        dt=dts[i],
        save_start=false,
        save_everystep=false,
        weak_timeseries_errors=false,
        weak_dense_errors=false,
        trajectories=numtraj) for i in 1:N]

errors = [LinearAlgebra.norm(Statistics.mean(sol.u)) for sol in _solutions]
m = log(errors[end]/errors[1])/log(dts[end]/dts[1])
@test abs(m-2) < 0.31

#convergence_plot = plot(dts, errors, xaxis=:log, yaxis=:log)
#savefig(convergence_plot, "RI1-CPU-final-"*string(numtraj)*".pdf")
#println(m)
