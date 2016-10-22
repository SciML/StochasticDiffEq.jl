
# Stochastic needs ΔW = s^(1/2).
#=
β = 0.6
σ = (t,y) -> β*y/(4.0s)
u₀ = 1.5Newton
prob = SDEProblem(f,σ,u₀)

sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:EM)
sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:SRIW1Optimized)

TEST_PLOT && plot(sol)

u₀ = [1.5Newton 2.0Newton
      3.0Newton 1.0Newton]

prob = SDEProblem(f,σ,u₀)

sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:EM)
sol =solve(prob::SDEProblem,[0,1],Δt=(1/2^4)Second,save_timeseries=true,alg=:SRIW1Optimized)
=#
