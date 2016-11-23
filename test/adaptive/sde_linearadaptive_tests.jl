using StochasticDiffEq, DiffEqProblemLibrary, DiffEqBase

probs = Vector{SDETestProblem}(2)
add_probs = Vector{SDETestProblem}(2)
probs[1] = prob_sde_2Dlinear
probs[2] = prob_sde_linear
add_probs[1] = prob_sde_additive
add_probs[2] = prob_sde_additivesystem


for i in 1:2

  ## SRIW1

  srand(100)
  sol =solve(probs[i],SRI(),dt=1//2^(4),abstol=1,reltol=0)
  err1 = sol.errors[:final]

  sol2 =solve(probs[i],SRI(),dt=1//2^(4),abstol=1e-1,reltol=0)
  err2 = sol2.errors[:final]

  sol3 =solve(probs[i],SRI(),dt=BigInt(1)//BigInt(2)^(4),abstol=1e-2,reltol=0)
  err3 = sol3.errors[:final]
  @test err1 > err3

  sol4 =solve(probs[i],SRI(),dt=1/2^(4),abstol=1e-3,reltol=0)
  err4 = sol4.errors[:final]
  @test err2 > err4

  srand(100)
  sol =solve(probs[i],SRIW1(),dt=1//2^(4),abstol=1,reltol=0)
  err21 = sol.errors[:final]
  @test err1 ≈ err21
  # p1 = plot(sol,plot_analytic=true)

  sol2 =solve(probs[i],SRIW1(),dt=1//2^(4),abstol=1e-1,reltol=0)
  err22 = sol2.errors[:final]
  @test err2 ≈ err22
  #TEST_PLOT && p2 = plot(sol2,plot_analytic=true)

  sol3 =solve(probs[i],SRIW1(),dt=BigInt(1)//BigInt(2)^(4),abstol=1e-2,reltol=0)
  err23 = sol3.errors[:final]
  @test err3 ≈ err23
  #TEST_PLOT && p3 = plot(sol3,plot_analytic=true)

  sol4 =solve(probs[i],SRIW1(),dt=1/2^(4),abstol=1e-3,reltol=0)
  err24 = sol4.errors[:final]
  @test err4 ≈ err24

  ## SRA1

  srand(100)
  sol =solve(add_probs[i],SRA(),dt=1//2^(4),abstol=1,reltol=0)
  err1 = sol.errors[:final]

  sol2 =solve(add_probs[i],SRA(),dt=1//2^(4),abstol=1e-1,reltol=0)
  err2 = sol2.errors[:final]

  sol3 =solve(add_probs[i],SRA(),dt=BigInt(1)//BigInt(2)^(4),abstol=1e-2,reltol=0)
  err3 = sol3.errors[:final]
  @test err1 > err3

  sol4 =solve(add_probs[i],SRA(),dt=1/2^(4),abstol=1e-4,reltol=0)
  err4 = sol4.errors[:final]
  @test err2 > err4

  srand(100)
  sol =solve(add_probs[i],SRA1(),dt=1//2^(4),abstol=1,reltol=0)
  err21 = sol.errors[:final]
  @test err1 ≈ err21
  # p1 = plot(sol,plot_analytic=true)

  sol2 =solve(add_probs[i],SRA1(),dt=1//2^(4),abstol=1e-1,reltol=0)
  err22 = sol2.errors[:final]
  @test err2 ≈ err22
  #TEST_PLOT && p2 = plot(sol2,plot_analytic=true)

  sol3 =solve(add_probs[i],SRA1(),dt=BigInt(1)//BigInt(2)^(4),abstol=1e-2,reltol=0)
  err23 = sol3.errors[:final]
  @test err3 ≈ err23
  #TEST_PLOT && p3 = plot(sol3,plot_analytic=true)

  sol4 =solve(add_probs[i],SRA1(),dt=1/2^(4),abstol=1e-4,reltol=0)
  err24 = sol4.errors[:final]
  @test isapprox(err4,err24;atol=1e-4)
end
