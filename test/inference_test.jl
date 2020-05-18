using StochasticDiffEq, Test
u0=1/2
f(u,p,t) = u
g(u,p,t) = u
dt = 1//2^(4)
tspan = (0.0,1.0)
prob = SDEProblem(f,g,u0,(0.0,1.0))
sol = solve(prob,EM(),dt=dt)
@inferred solve(prob,EM(),dt=dt)
