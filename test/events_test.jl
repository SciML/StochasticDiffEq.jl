using StochasticDiffEq, Test
function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -9.81
end

function g(du,u,p,t)
  nothing
end

function condition(u,t,integrator) # Event when event_f(u,p,t,k) == 0
  u[1]
end

affect! = nothing
function affect_neg!(integrator)
  integrator.u[2] = -integrator.u[2]
end

# Continuous callback
callback = ContinuousCallback(condition,affect!,affect_neg!)

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
times_finalize_called = 0
callback_dc = DiscreteCallback(condition_dc, affect!_dc, save_positions=save_positions,
  finalize=(args...)->global times_finalize_called+=1)
sol = solve(prob, SRIW1(), callback=callback_dc, tstops=tstop, saveat=tstop)
@test count(x->x==tstop[1], sol.t) == 2
@test count(x->x==tstop[2], sol.t) == 2
@test times_finalize_called == 1
sol = solve(prob, SRIW1(), callback=callback_dc, tstops=tstop, saveat=prevfloat.(tstop))
@test count(x->x==tstop[1], sol.t) == 2
@test count(x->x==tstop[2], sol.t) == 2
@test times_finalize_called == 2



###
# https://github.com/SciML/DifferentialEquations.jl/issues/802
###

function HM_neuron!(du,u,Params,t)
    # Membrane voltage
    du[1] = -u[1]^3 + 3.0 * u[1]^2 + u[2] - u[3] + u[4]*u[1];

    # Fast channel
    du[2] = 1.0 - 5.0 * u[1]^2 - u[2];

    # Slow channel
    val = Params.s * (u[1] + 8.0/5.0);
    du[3] = Params.r * (val - u[3]);

    # Synapse
    du[4] = -u[4]/Params.syntau;
end

function HM_noise!(du,u,Params,t)
    du[1] = 0.1;
    du[2] = 0.0;
    du[3] = 0.0;
    du[4] = 0.0;
end

tvals = range(0.0,stop=1999.9,length=20000);
struct P
    r::Float64
    s::Float64
    I::Float64
    syntau::Float64
    gSyn::Float64
end

params = P(0.001,1,0.0,5,0.0);
x0 = Float64[-1.6,-11.8,0.0,0.0];

condition2(u,t,integrator) = u[1];
affect2!(integrator) = integrator.u[4] = integrator.u[4] - params.gSyn; # Inhibitory exponential decay synapse
cb = ContinuousCallback(condition2,affect2!,nothing);
prob = SDEProblem(HM_neuron!,HM_noise!,x0,(0.0,1999.9),params);
sol  = solve(prob,ImplicitEM(),reltol = 1e-4, abstol = 1e-6,dense = true,callback=cb);
sol  = solve(prob,SKenCarp(),reltol = 1e-4, abstol = 1e-6,dense = true,callback=cb);

using DiffEqCallbacks

function f(du,u,p,t)
  du[1] = p[1] -u[1]
end

function g(du,u,p,t)
  du[1] = p[2]
end

sprob = SDEProblem(f,g,[1.0],(0.0,10.0),[1.0,0.1])
sol = solve(sprob,ImplicitEM(),callback=PositiveDomain())
