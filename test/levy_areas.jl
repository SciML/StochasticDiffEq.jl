using StochasticDiffEq, DiffEqNoiseProcess, Test, Random, LinearAlgebra

function true_general_function(dt,dW,C,m)
  # Choose p
  sum_dW² = dot(dW,dW)
  M = div(m*(m-1),2)
  p = Int(floor(sqrt(M/(12*dt*C))*sqrt(m + 4*sum_dW²/dt)/π + 1))
  # Alternative choice of p which on average results in larger values:
  #p = Int(floor(sqrt(5*m*M/(12*dt*C))/π+1))
  #@show p
  # a_p below (22) for tail approx.
  ap = pi^2/6
  for r=1:p
    ap -= 1/r^2
  end

  # define permutation operator
  PermOp = zeros(Int,m^2,m^2)
  for i=1:m^2
    j = 1 + m*((i-1)%m) + div(i-1,m)
    PermOp[i,j] = 1
  end

  # define K operator
  KOp = zeros(Int,M,m^2)
  row = 1
  for i=1:m-1
    col = (i-1)*m + i + 1
    s = m-i
    KOp[row:(row+s)-1,col:(col+s)-1] = Matrix{eltype(dW)}(I, s, s)
    row += s
  end

  # define identiy matrices
  Idm = Matrix{eltype(dW)}(I, m, m)
  IdM = Matrix{eltype(dW)}(I, M, M)
  Idm2 = Matrix{eltype(dW)}(I, m^2, m^2)

  # compute Σinf
  op1 = 2/dt * KOp*(Idm2-PermOp)
  op2 = kron(Idm, dW .* dW') #kron(I+zeros(m,m),W.dW .* W.dW')
  op3 = (Idm2-PermOp)*KOp'
  Σinf = 2*IdM + op1*op2*op3

  # compute sqrt( Σinf ) from Cholesky decomposition
  α = sqrt(1 + sum_dW²/dt)
  SqΣinf = (Σinf + 2*α*IdM)/(sqrt(2)*(1+α))

  # initialze Gp from distribution N(0, I_M) to compute tail approx
  Gp₁ = randn(M)
  prefac = (dt/(2*π))*sqrt(ap)
  Atail = prefac*SqΣinf*Gp₁ # check influence of op3 here!!

  #@show Atail

  # compute Atilde
  Ap = zeros(M)

  # for r=1:p
  #   sd = sqrt(dt/(2*pi*r))
  #   ζ = randn(m)*sd
  #   η = randn(m)*sd
  #   aᵢ₀ -= (2/sqrt(π*r))*ζ
  #   for i=1:m
  #     for j=1:m
  #       Ap[i,j] += ζ[i]*η[j]-η[i]*ζ[j]
  #     end
  #   end
  # end
  # for i=1:m
  #   for j=1:m
  #     Jp[i,j] += Ap[i,j] - 1//2*(aᵢ₀[j]*dW[i]-dW[j]*aᵢ₀[i])
  #   end
  # end

  op4 = KOp*(PermOp-Idm2)

  for k in 1:p
    ζ = randn(eltype(dW),m)
    η = randn(eltype(dW),m)
    #@show ζ , η
    #@show ζ * sqrt(dt/(2*π*k)), η * sqrt(dt/(2*π*k))
    #kronprod1 = vec(ζ) .* vec(η+convert(eltype(dW),sqrt(2/dt))*dW)'
    kronprod2 = vec(η+convert(eltype(dW),sqrt(2/dt))*dW) .* vec(ζ)'

    #@show size(vec(kronprod2)), size(op4), size(Ap)
    Ap += op4*vec(kronprod2)/k #(kronprod1 - kronprod2)/k
  end

  Ap *= dt/(2*π)

  # compute J_{ij}^p
  Jp_com = vec(1//2 .* vec(dW) .* vec(dW)')

  Jp =  Jp_com + op3*(Ap + Atail)

  #@show op3*Atail

  return reshape(Jp, (m, m)), reshape((op3*Ap), (m, m)), reshape(op3*Atail, (m, m))
end

seed = 10
Random.seed!(seed)

m = 10
W = WienerProcess(0.0,zeros(m),nothing)

dt = 0.1
calculate_step!(W,dt,nothing,nothing)

for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

# Diagonal noise tests

true_diag = 1//2 .* W.dW .* W.dW

Wikdiag = StochasticDiffEq.WikJDiagonal_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikdiag, nothing, 1)
Wikdiagoop = StochasticDiffEq.WikJDiagonal_oop()

@test StochasticDiffEq.get_iterated_I(dt, W.dW, Wikdiagoop, nothing, 1) == Wikdiag.WikJ
@test Wikdiag.WikJ == true_diag


KPWdiagiip = StochasticDiffEq.KPWJ_iip(W.dW)
Random.seed!(seed)
StochasticDiffEq.get_iterated_I!(dt, W.dW, KPWdiagiip, Int(1e3), 1, 1//1)
KPWdiagoop = StochasticDiffEq.KPWJ_oop()
Random.seed!(seed)
@test isapprox(StochasticDiffEq.get_iterated_I(dt, W.dW, KPWdiagoop, Int(1e3), 1, 1//1), KPWdiagiip.WikJ, atol=1e-14)

KPWdiagonly = [KPWdiagiip.WikJ[i, i] for i in 1:m]

@test KPWdiagonly == true_diag
@test KPWdiagonly == Wikdiag.WikJ


# Commutative noise tests

true_commute = 1//2 .* W.dW .* W.dW'

Wikcommute = StochasticDiffEq.WikJCommute_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikcommute, nothing, 1)
Wikcommuteoop = StochasticDiffEq.WikJCommute_oop()

@test StochasticDiffEq.get_iterated_I(dt, W.dW, Wikcommuteoop, nothing, 1) == Wikcommute.WikJ
@test Wikcommute.WikJ == true_commute


# General noise tests

Random.seed!(seed)
Wikgeneral = StochasticDiffEq.WikJGeneral_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikgeneral, nothing, 1)

Random.seed!(seed)
Wikgeneraloop  = StochasticDiffEq.WikJGeneral_oop(W.dW)
@test StochasticDiffEq.get_iterated_I(dt, W.dW, Wikgeneraloop, nothing, 1) == Wikgeneral.WikJ

Random.seed!(seed)
true_noncom = true_general_function(dt, W.dW, 1.0, m)

@test isapprox(Wikgeneral.WikJ, true_noncom[1], atol=1e-15)
@test isapprox(Wikgeneral.WikJ2, true_noncom[3], atol=1e-15)


# Test the relations given in Wiktorsson Eq.(2.1)
Aij = true_noncom[2]+true_noncom[3]
Aii = [Aij[i, i] for i in 1:m]
Jii = [true_noncom[1][i, i] for i in 1:m]

@test isapprox(true_noncom[1] + true_noncom[1]', 2*true_commute)
@test isapprox(true_noncom[1], true_commute + Aij)
@test isapprox(Aij, -Aij', atol=1e-15)
@test isapprox(Jii, 1//2 .* W.dW .* W.dW, atol=1e-15)
@test isapprox(Aii, zeros(m), atol=1e-15)


# .. for KPW scheme
Aii = [KPWdiagiip.WikA[i, i] for i in 1:m]
@test isapprox(KPWdiagiip.WikJ + KPWdiagiip.WikJ', 2*true_commute)
@test !isapprox(KPWdiagiip.WikJ, true_commute + KPWdiagiip.WikA) # ! because of correction with \mu term
@test isapprox(KPWdiagiip.WikA, -KPWdiagiip.WikA', atol=1e-15)
@test isapprox(KPWdiagonly, 1//2 .* W.dW .* W.dW, atol=1e-15)
@test isapprox(Aii, zeros(m), atol=1e-15)



# moment conditions
# E1(I_{j1} I_{j2}) = Δ δ_{j1,j2}
# E2(I_{j1, j2} I_{j1, j2}) = 1/2 Δ^2
# E3(I_{j1} I_{j2} I_{j3, j4}) = {Δ^2 if j1=j2=j3=j4, 1/2 Δ^2 if j3!=j4 and j1=j3, j2=j4 or j1=j4, j2=j3, 0 otherwise}

function moments!(tmp, E1, E2, E3, W, Wik, Δ, samples, p=nothing)
  for _ in 1:samples
    calculate_step!(W,Δ,nothing,nothing)
    accept_step!(W,Δ,nothing,nothing)
    mul!(tmp,vec(W.dW),vec(W.dW)')

    StochasticDiffEq.get_iterated_I!(Δ, W.dW, Wik, p, 1)
    @. E1 = E1 + tmp
    @. E2 = E2 + Wik.WikJ * Wik.WikJ
    @. E3 = E3 + tmp * Wik.WikJ
  end
  @. E1 = E1/samples
  @. E2 = E2/samples
  @. E3 = E3/samples

  return nothing
end


samples = Int(1e5)
E1 = false .* vec(W.dW) .* vec(W.dW)'
E2 = zero(E1)
E3 = zero(E1)
tmp = zero(E1)
# generate new random dW
W = WienerProcess!(0.0,zeros(m),nothing)
Wikmom = StochasticDiffEq.WikJGeneral_iip(W.dW)
@time moments!(tmp, E1, E2, E3, W, Wikmom, dt, samples)

@test maximum(abs.(E1 - dt*I)) < 1e-3
@test maximum(abs.(E2- dt^2*(zero(E2).+1)/2)) < 4e-3
@test maximum(abs.(E3 - dt^2*(one(E2).+1)/2)) < 6e-3


"""
Problem 2.3.3 from
Kloeden, P. E., Platen, E., & Schurz, H. Numerical solution of SDE through computer
experiments. Springer Science & Business Media. (2012)
"""

Random.seed!(seed)
dt = 1.0
ps = [1,Int(1e2),Int(5e3),Int(1e4)]
m = 2
W2 = WienerProcess(0.0,zeros(m),nothing)
calculate_step!(W2,dt,nothing,nothing)
for i in 1:10
  accept_step!(W2,dt,nothing,nothing)
end

function path_convergence(dt, ps, W, Wik)
    sample_path = []
    for (i, p) in enumerate(ps)
      Random.seed!(seed)
      @show p
      StochasticDiffEq.get_iterated_I!(dt, W.dW, Wik, p, 1)
      #@show Wik.WikJ
      push!(sample_path,Wik.WikJ[1,2])
    end
    return sample_path
end

Wikpath = StochasticDiffEq.KPWJ_iip(W2.dW)
path = path_convergence(dt, ps, W2, Wikpath)

v = abs.(path .- path[end])
@test sort(v, rev = true) == v


"""
Exercise 2.3.1 from
Kloeden, P. E., Platen, E., & Schurz, H. Numerical solution of SDE through computer
experiments. Springer Science & Business Media. (2012)
"""


mutable struct StatsJ{AType}
  mean::Vector{AType}
  var::Vector{AType}
  meanold::Vector{AType}
  varold::Vector{AType}
  tmp::AType
  x::AType
end

function StatsJ(u)
  mean = Vector{typeof(u)}()
  var = Vector{typeof(u)}()
  meanold = Vector{typeof(u)}()
  varold = Vector{typeof(u)}()

  for k=1:2
    push!(mean,zero(u))
    push!(var,zero(u))
    push!(meanold,zero(u))
    push!(varold,zero(u))
  end

  tmp = zero(u)
  x = zero(u)
  StatsJ(mean,var,meanold,varold,tmp,x)
end

function Welford!(cache::StatsJ, W, Wik, Δ, samples, p=nothing)
  for k in 1:2
    fill!(cache.var[k], zero(eltype(cache.tmp)))
  end
  for i in 1:samples
    calculate_step!(W,Δ,nothing,nothing)
    accept_step!(W,Δ,nothing,nothing)
    mul!(cache.tmp,1//2*vec(W.dW),vec(W.dW)')

    StochasticDiffEq.get_iterated_I!(Δ, W.dW, Wik, p, 1)
    for k in 1:2
      if k==1
        copyto!(cache.x,cache.tmp)
      else
        copyto!(cache.x,Wik.WikJ)
      end
      if i == 1
        copyto!(cache.meanold[k], cache.x)
      else
        @. cache.mean[k] = cache.meanold[k] + (cache.x - cache.meanold[k]) / i
        @. cache.var[k] = cache.varold[k] + (cache.x - cache.meanold[k]) * (cache.x - cache.mean[k])

        copyto!(cache.meanold[k],cache.mean[k])
        copyto!(cache.varold[k],cache.var[k])
      end
    end
  end

  for k in 1:2
    @. cache.mean[k] = cache.mean[k]/samples
    @. cache.var[k] = cache.var[k]/(samples-1)
  end

  return nothing
end

cache = StatsJ(false .* vec(W2.dW) .* vec(W2.dW)')
@time Welford!(cache, W2, Wikpath, dt, Int(1e5), ps[2])

@test maximum(abs.(cache.mean[1])) < 1e-5
@test maximum(abs.(cache.mean[2])) < 1e-5
@test maximum(abs.(cache.var[1]-cache.var[2])) > 1e-1
