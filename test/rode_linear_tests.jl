using DiffEqBase, StochasticDiffEq

f = (t,u,W) -> 1.01u.+0.87u.*W
(p::typeof(f))(::Type{Val{:analytic}},t,u0,W) = u0.*exp.(0.63155*t+0.87*W)
u0 = 1.00
tspan = (0.0,1.0)
prob = RODEProblem(f,u0,tspan)

sol = solve(prob,RandomEM(),dt=1/100)

f = (t,u,W,du) -> du.=1.01u.+0.87u.*W
(p::typeof(f))(::Type{Val{:analytic}},t,u0,W) = u0.*exp.(0.63155*t+0.87*W)
u0 = ones(4)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)
