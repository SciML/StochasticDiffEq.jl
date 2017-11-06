using StochasticDiffEq, DiffEqNoiseProcess

f(t,u,du) = (du .= u)
g(t,u,du) = (du .= u)
u0 = rand(4,2)
W = WienerProcess(0.0,0.0,0.0)
prob = SDEProblem(f,g,u0,(0.0,1.0),noise=W)
sol = solve(prob,SRIW1())

similar(0.0)
