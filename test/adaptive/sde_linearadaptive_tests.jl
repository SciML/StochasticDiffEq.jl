using StochasticDiffEq, DiffEqProblemLibrary, DiffEqBase
srand(100)
prob = prob_sde_linear

f = (t,u) -> u
σ = (t,u) -> u
analytic = (t,u0,W) -> u0.*exp.(.5*t+W)

prob_sde_linear = SDETestProblem(f,σ,1/2,analytic)

#using Plots
#gr(); Test_PLOT = true

## Solve and plot
println("Solve and Plot")
sol =solve(prob,SRI(),dt=1//2^(4),abstol=1,reltol=0)
err1 = sol.errors[:final]
println("Final error for the first solution was $err1")
# p1 = plot(sol,plot_analytic=true)

sol2 =solve(prob,SRI(),dt=1//2^(4),abstol=1e-1,reltol=0)
err2 = sol2.errors[:final]
println("Final error for the second solution was $err2")
#TEST_PLOT && p2 = plot(sol2,plot_analytic=true)

sol3 =solve(prob,SRI(),dt=BigInt(1)//BigInt(2)^(4),abstol=1e-2,reltol=0)
err3 = sol3.errors[:final]
println("Final error for the third solution was $err3")
#TEST_PLOT && p3 = plot(sol3,plot_analytic=true)

sol4 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-3,reltol=0)
err4 = sol4.errors[:final]
println("Final error for the fourth solution was $err4")
#TEST_PLOT && p4 = plot(sol4,plot_analytic=true)

#plot(p1,p2,p3,p4,title="Solutions to Linear SDE at Different Tolerances")
#gui()

true
