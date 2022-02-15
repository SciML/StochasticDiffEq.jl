using StochasticDiffEq, DiffEqNoiseProcess, Test, Random, LinearAlgebra
using LevyArea, Statistics

seed = 113213898
Random.seed!(seed)

m = 10
W = WienerProcess(0.0,zeros(m),nothing)

dt = 0.1
calculate_step!(W,dt,nothing,nothing)

for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

# LevyArea.jl algs: ITER_INT_ALGS

@testset "diagonal noise tests" begin
  true_diag = 1//2 .* W.dW .* W.dW

  Wikdiag = StochasticDiffEq.WikJDiagonal_iip(W.dW)
  StochasticDiffEq.get_iterated_I!(dt, W.dW, W.dZ, Wikdiag, nothing, 1)
  Wikdiagoop = StochasticDiffEq.WikJDiagonal_oop()

  @test StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, Wikdiagoop, nothing, 1) == Wikdiag.WikJ
  @test Wikdiag.WikJ == true_diag

  for alg ∈ LevyArea.ITER_INT_ALGS
    @test diag(StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, alg, nothing, 1)) == true_diag
  end   
end

# Commutative noise tests
@testset "Commutative noise tests" begin
  true_commute = 1//2 .* W.dW .* W.dW'

  Wikcommute = StochasticDiffEq.WikJCommute_iip(W.dW)
  StochasticDiffEq.get_iterated_I!(dt, W.dW, W.dZ, Wikcommute, nothing, 1)
  Wikcommuteoop = StochasticDiffEq.WikJCommute_oop()

  @test StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, Wikcommuteoop, nothing, 1) == Wikcommute.WikJ
  @test Wikcommute.WikJ == true_commute

  for alg ∈ LevyArea.ITER_INT_ALGS
    @test diag(StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, alg, nothing, 1)) == diag(true_commute)
  end   
end


# moment conditions
# E1(I_{j1} I_{j2}) = Δ δ_{j1,j2}
# E2(I_{j1, j2} I_{j1, j2}) = 1/2 Δ^2
# E3(I_{j1} I_{j2} I_{j3, j4}) = {Δ^2 if j1=j2=j3=j4, 1/2 Δ^2 if j3!=j4 and j1=j3, j2=j4 or j1=j4, j2=j3, 0 otherwise}
function test_moments(m, Wik, Δ, samples, p=nothing)
  # generate new random dW
  W = WienerProcess!(0.0,zeros(m),nothing)
  E1 = false .* vec(W.dW) .* vec(W.dW)'
  E2 = zero(E1)
  E3 = zero(E1)
  tmp = zero(E1)
  for _ in 1:samples
    calculate_step!(W,Δ,nothing,nothing)
    accept_step!(W,Δ,nothing,nothing)
    mul!(tmp,vec(W.dW),vec(W.dW)')

    I = StochasticDiffEq.get_iterated_I(Δ, W.dW, W.dZ, Wik, p, 1)
    @. E1 = E1 + tmp
    @. E2 = E2 + I * I
    @. E3 = E3 + tmp * I
  end
  @. E1 = E1/samples
  @. E2 = E2/samples
  @. E3 = E3/samples

  @test maximum(abs.(E1 - Δ*I)) < 2e-1
  @test maximum(abs.(E2- Δ^2*(zero(E2).+1)/2)) < 1e-2
  @test maximum(abs.(E3 - Δ^2*(one(E2).+1)/2)) < 1e-2
end


"""
Problem 2.3.3 from
Kloeden, P. E., Platen, E., & Schurz, H. Numerical solution of SDE through computer
experiments. Springer Science & Business Media. (2012)
"""
function test_path_convergence(Wik, dt = 1.0, ps = [Int(1e1),Int(1e3),Int(1e6)])  
  m = 2
  W = WienerProcess(0.0,zeros(m),nothing)
  calculate_step!(W,dt,nothing,nothing)
  for i in 1:10
    accept_step!(W,dt,nothing,nothing)
  end  
  sample_path = []
  for (i, p) in enumerate(ps)
    Random.seed!(seed)
    I = StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, Wik, p, 1)
    A = I - 1//2 .* W.dW .* W.dW'
    push!(sample_path,A[1,2])
  end
  v = abs.(sample_path .- sample_path[end])
  @test sort(v, rev = true) == v
end

"""
Exercise 2.3.1 from
Kloeden, P. E., Platen, E., & Schurz, H. Numerical solution of SDE through computer
experiments. Springer Science & Business Media. (2012)
"""

function test_compare_sample_mean_and_var(Wik, Δ, m, samples=Int(5e5), p=Int(1e2))
  W = WienerProcess(0.0,zeros(m),nothing)
  calculate_step!(W,dt,nothing,nothing)
  for i in 1:10
    accept_step!(W,dt,nothing,nothing)
  end  

  xs = [StochasticDiffEq.get_iterated_I(Δ, W.dW, W.dZ, Wik, p, 1)[1,2]]
  coms = [1//2*W.dW[1]*W.dW[2]]

  for i in 1:samples
    calculate_step!(W,Δ,nothing,nothing)
    accept_step!(W,Δ,nothing,nothing)
    com = 1//2*W.dW[1]*W.dW[2]
    x = StochasticDiffEq.get_iterated_I(Δ, W.dW, W.dZ, Wik, p, 1)[1,2]
    
    push!(xs,x)
    push!(coms,com)
  end

  @test mean(xs) ≈ 0 atol = 1e-2
  @test mean(coms) ≈ 0 atol = 1e-2
  @test var(xs) ≈ 1//2*Δ^2 rtol = 1e-2
  @test var(coms) ≈ 1//4*Δ^2 rtol = 1e-2
end

@testset "General noise tests" begin
  true_commute = 1//2 .* W.dW .* W.dW'
  samples = Int(1e4)
  for alg ∈ LevyArea.ITER_INT_ALGS
    # Test the relations given in Wiktorsson Eq.(2.1)
    Random.seed!(seed)
    I = StochasticDiffEq.get_iterated_I(dt, W.dW, W.dZ, alg, 1, 1)
    Random.seed!(seed)
    A = LevyArea.levyarea(W.dW/√dt, 1, alg) 
    @test dt*A + true_commute == I # because of correction with \mu term
    @test diag(A) == diag(zero(A))
    @test diag(true_commute) == diag(I)
    @test I + I' ≈ 2*true_commute atol=1e-12
    @test A ≈ -A' atol=1e-12

    # test moment conditions
    test_moments(m, alg, dt, samples)

    # test sample path convergence
    test_path_convergence(alg)

    # test other StatsJ
    @time test_compare_sample_mean_and_var(alg, 1.0, 2)
  end  
end
