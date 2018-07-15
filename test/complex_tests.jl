using DiffEqBase, StochasticDiffEq

f(u,p,t)=1.01u
g(u,p,t)=1.01u
u0 = 1.0+1.0im
tspan = (0.0,1.0)
prob = SDEProblem(f,g,u0,tspan)
sol = solve(prob,SRIW1())
typeof(sol[end]) == ComplexF64

f(du,u,p,t) = (du.=1.01u)
g(du,u,p,t) = (du.=1.01u)
u0 = ones(2,4) + im*ones(2,4)
prob = SDEProblem(f,g,u0,tspan)
sol = solve(prob,SRIW1())
eltype(sol[end]) == ComplexF64
