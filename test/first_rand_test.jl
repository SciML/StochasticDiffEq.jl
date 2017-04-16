using DiffEqBase, StochasticDiffEq, Base.Test

f1 = (t,u) -> 0.
g1 = (t,u) -> 1.
dt = 1//2^(4)
prob1 = SDEProblem(f1,g1,0.,(0.0,1.0))
integrator = init(prob1,EM(),dt=dt)

k = integrator.ΔW
@test integrator.ΔW != 0
solve!(integrator)
@test integrator.sol.W[2] == k

prob1 = SDEProblem(f1,g1,zeros(4),(0.0,1.0))
sol = solve(prob1,EM(),dt=dt)
@test sol.W[2] != zeros(4)


f1 = (t,u,du) -> du.=0.
g1 = (t,u,du) -> du.=1.
dt = 1//2^(4)
prob1 = SDEProblem(f1,g1,zeros(4),(0.0,1.0))
sol = solve(prob1,EM(),dt=dt)
@test sol.W[2] != zeros(4)
