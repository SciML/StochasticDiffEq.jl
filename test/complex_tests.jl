using DiffEqBase, StochasticDiffEq

f(t,u)=1.01u
g(t,u)=1.01u
u0 = 1.0+1.0im
tspan = (0.0,1.0)
prob = SDEProblem(f,g,u0,tspan)
sol = solve(prob,SRIW1())
typeof(sol[end]) == Complex128

f(t,u,du) = (du.=1.01u)
g(t,u,du) = (du.=1.01u)
u0 = ones(2,4) + im*ones(2,4)
prob = SDEProblem(f,g,u0,tspan)
sol = solve(prob,SRIW1())
eltype(sol[end]) == Complex128
