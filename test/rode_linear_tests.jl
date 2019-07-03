using StochasticDiffEq

f(u,p,t,W) = 1.01u.+0.87u.*W
u0 = 1.00
tspan = (0.0,1.0)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)

f(u,p,t,W,du) = (du.=1.01u.+0.87u.*W)
u0 = ones(4)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)

f(u,p,t,W) = 2u*sin(W)
u0 = 1.00
tspan = (0.0,5.0)
prob = RODEProblem{false}(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)

function f(du,u,p,t,W)
  du[1] = 2u[1]*sin(W[1] - W[2])
  du[2] = -2u[2]*cos(W[1] + W[2])
end
u0 = [1.00;1.00]
tspan = (0.0,5.0)
prob = RODEProblem(f,u0,tspan)
sol = solve(prob,RandomEM(),dt=1/100)

function f(du,u,p,t,W)
  du[1] = -2W[3]*u[1]*sin(W[1] - W[2])
  du[2] = -2u[2]*cos(W[1] + W[2])
end
u0 = [1.00;1.00]
tspan = (0.0,5.0)
prob = RODEProblem(f,u0,tspan,rand_prototype=zeros(3))
sol = solve(prob,RandomEM(),dt=1/100)
