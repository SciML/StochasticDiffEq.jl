using StochasticDiffEq
srand(100)
prob = prob_sde_2Dlinear

using Plots; gr()
## Solve and plot
println("Solve and Plot")
sol =solve(prob,SRI(),dt=1/2^(4),abstol=1,reltol=0)
err1 = sol.errors[:final]
p1 = plot(sol,plot_analytic=true,legend=false,title="tol = 1")

println("1e-1")
sol2 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-1,reltol=0)
err2 = sol2.errors[:final]
p2 = plot(sol2,plot_analytic=true,legend=false,title="tol = 1e-1")

println("1e-2")
sol3 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-2,reltol=0)
err3 = sol3.errors[:final]
p3 = plot(sol3,plot_analytic=true,legend=false,title="tol = 1e-2")

println("1e-3")
sol4 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-3,reltol=0)
err4 = sol4.errors[:final]
p4 = plot(sol4,plot_analytic=true,legend=false,title="tol = 1e-3")

plot(p1,p2,p3,p4,title="Solutions to Linear SDE at Different Tolerances",size=(1200,800))
#gui()

println("""
Final error for the solutions were:
          $err1
          $err2
          $err3
          $err4""")

err4 < err2 && err3 < err1
