using StochasticDiffEq
srand(100)
prob = prob_sde_2Dlinear

## Solve and plot
sol =solve(prob,SRI(),dt=1/2^(4),abstol=1,reltol=0)
err1 = sol.errors[:final]

println("1e-1")
sol2 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-1,reltol=0)
err2 = sol2.errors[:final]

println("1e-2")
sol3 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-2,reltol=0)
err3 = sol3.errors[:final]
@test err3 < err1
println("1e-3")
sol4 =solve(prob,SRI(),dt=1/2^(4),abstol=1e-3,reltol=0)
err4 = sol4.errors[:final]
@test err4 < err2
println("""
Final error for the solutions were:
          $err1
          $err2
          $err3
          $err4""")
