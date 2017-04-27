using StochasticDiffEq, RecursiveArrayTools, DiffEqBase, Base.Test#, ParameterizedFunctions
f = function (t,u,du)
  du[1] = u[2]
  du[2] = -9.81
end

g = function (t,u,du)
  nothing
end

condtion= function (t,u,integrator) # Event when event_f(t,u,k) == 0
  u[1]
end

affect! = nothing
affect_neg! = function (integrator)
  integrator.u[2] = -integrator.u[2]
end

callback = ContinuousCallback(condtion,affect!,affect_neg!)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
prob = SDEProblem(f,g,u0,tspan)

sol = solve(prob,SRIW1(),callback=callback,adaptive=false,dt=3/4)

@test sol[1,6] < 1e-14

g = function (t,u,du)
  du[2] = .125*u[2]
end

prob = SDEProblem(f,g,u0,tspan)

sol = solve(prob,SRIW1(),callback=callback)

sol = solve(prob,EM(),callback=callback,dt=1/4)
