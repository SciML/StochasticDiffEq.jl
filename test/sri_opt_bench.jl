using StochasticDiffEq, DiffEqProblemLibrary, DiffEqBase
srand(200)
prob = oval2ModelExample(largeFluctuations=true,useBigs=false)
quick_prob = deepcopy(prob)
quick_prob.tspan = (0.0,1.0)

using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.0001
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
@benchmark begin
  srand(100)
  sol = solve(quick_prob,SRI(),dt=(1/2)^(18),progress_steps=Int(1e5),
        qmax=1.125,save_timeseries=false,
        timeseries_steps=1000,abstol=1e-5,reltol=1e-4)
end

@benchmark begin
  srand(100)
  sol = solve(quick_prob,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),
        progress_steps=Int(1e5),
        qmax=1.125,save_timeseries=false,
        timeseries_steps=1000,abstol=1e-5,reltol=1e-4)
end


srand(200)
@time sol = solve(prob,SRI(tableau=StochasticDiffEq.constructSRIOpt1()),dt=(1/2)^(18),progress=true,qmax=1.125,
          timeseries_steps=1000,abstol=1e-5,reltol=1e-3);
using Plots; gr()
lw = 2
lw2 = 3
p1 = plot(sol.t,sol[:,16],top_margin=50px,
          title="(A) Timeseries of Ecad Concentration",xguide="Time (s)",
          yguide="Concentration",guidefont=font(16),tickfont=font(16),
          linewidth=lw,left_margin=85px,leg=false)

p2 = plot(sol.t,sol[:,17],top_margin=50px,
          title="(B) Timeseries of Vim Concentration",xguide="Time (s)",
          yguide="Concentration",guidefont=font(16),
          tickfont=font(16),linewidth=lw,leg=false)
