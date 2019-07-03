using StochasticDiffEq, Test
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

# Continuous callback
callback = ContinuousCallback(condtion,affect!,affect_neg!)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
prob = SDEProblem(f,g,u0,tspan)

sol = solve(prob,SRIW1(),callback=callback,adaptive=false,dt=3/4)

@test minimum(sol[1,:]) > -1e-12 && minimum(sol[1,:]) < 1e-12

sol = solve(prob,SRIW1(),callback=callback,save_everystep=false)
t = sol.t[end÷2] # this is the callback time point
sol = solve(prob,SRIW1(),callback=callback,saveat=t)
@test count(x->x==t, sol.t) == 2
sol = solve(prob,SRIW1(),callback=callback,saveat=t-eps(t))
@test count(x->x==t, sol.t) == 2

function g(du,u,p,t)
  du[2] = .125*u[2]
end

prob = SDEProblem(f,g,u0,tspan)

sol = solve(prob,SRIW1(),callback=callback)

sol = solve(prob,EM(),callback=callback,dt=1/4)

# Discrete callback
tstop = [5.;8.]
condition_dc = (u,t,integrator) -> t in tstop
affect!_dc = (integrator) -> integrator.u .= 1.0
save_positions = (true,true)
callback_dc = DiscreteCallback(condition_dc, affect!_dc, save_positions=save_positions)
sol = solve(prob, SRIW1(), callback=callback_dc, tstops=tstop, saveat=tstop)
@test count(x->x==tstop[1], sol.t) == 2
@test count(x->x==tstop[2], sol.t) == 2
sol = solve(prob, SRIW1(), callback=callback_dc, tstops=tstop, saveat=prevfloat.(tstop))
@test count(x->x==tstop[1], sol.t) == 2
@test count(x->x==tstop[2], sol.t) == 2
