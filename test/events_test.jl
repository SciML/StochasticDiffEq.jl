using StochasticDiffEq, RecursiveArrayTools, DiffEqBase, Base.Test#, ParameterizedFunctions
function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -9.81
end

function g(du,u,p,t)
  nothing
end

function condtion(u,t,integrator) # Event when event_f(u,p,t,k) == 0
  u[1]
end

affect! = nothing
function affect_neg!(integrator)
  integrator.u[2] = -integrator.u[2]
end

callback = ContinuousCallback(condtion,affect!,affect_neg!)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
prob = SDEProblem(f,g,u0,tspan)

sol = solve(prob,SRIW1(),callback=callback,adaptive=false,dt=3/4)

@test minimum(sol[1,:]) > -1e-12 && minimum(sol[1,:]) < 1e-12

function g(du,u,p,t)
  du[2] = .125*u[2]
end

prob = SDEProblem(f,g,u0,tspan)

sol = solve(prob,SRIW1(),callback=callback)

sol = solve(prob,EM(),callback=callback,dt=1/4)
