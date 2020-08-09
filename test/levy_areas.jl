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
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikdiag, 1.0)
Wikdiagoop = StochasticDiffEq.WikJDiagonal_oop()

@test StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikdiagoop, 1.0) == Wikdiag.WikJ
@test Wikdiag.WikJ == true_diag


KPWdiagiip = StochasticDiffEq.KPWJ_iip(W.dW)
Random.seed!(seed)
StochasticDiffEq.get_iterated_I!(dt, W.dW, KPWdiagiip, Int(1e3))
KPWdiagoop = StochasticDiffEq.KPWJ_oop()
Random.seed!(seed)
@test isapprox(StochasticDiffEq.get_iterated_I!(dt, W.dW, KPWdiagoop, Int(1e3)), KPWdiagiip.WikJ, atol=1e-14)

KPWdiagonly = [KPWdiagiip.WikJ[i, i] for i in 1:m]

@test KPWdiagonly == true_diag
@test KPWdiagonly == Wikdiag.WikJ


# Commutative noise tests

true_commute = 1//2 .* W.dW .* W.dW'

Wikcommute = StochasticDiffEq.WikJCommute_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikcommute, 1.0)
Wikcommuteoop = StochasticDiffEq.WikJCommute_oop()

@test StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikcommuteoop, 1.0) == Wikcommute.WikJ
@test Wikcommute.WikJ == true_commute


# General noise tests

Random.seed!(seed)
Wikgeneral = StochasticDiffEq.WikJGeneral_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikgeneral, 1.0)

Random.seed!(seed)
Wikgeneraloop  = StochasticDiffEq.WikJGeneral_oop(W.dW)
@test StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikgeneraloop, 1.0) == Wikgeneral.WikJ

Random.seed!(seed)
true_noncom = true_general_function(dt, W.dW, 1.0, m)

@test isapprox(Wikgeneral.WikJ, true_noncom[1], atol=1e-15)
@test isapprox(Wikgeneral.WikJ2, true_noncom[3], atol=1e-15)
