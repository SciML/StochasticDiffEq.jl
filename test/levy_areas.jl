using StochasticDiffEq, DiffEqNoiseProcess, Test, Random, LinearAlgebra

function true_general_function(dt,dW,C,m)

  # (1) Choose p
  sum_dW² = dot(dW,dW)
  M = div(m*(m-1),2)
  p = Int(floor(sqrt(M/(12*dt*C))*sqrt(m + 4*sum_dW²/dt)/π + 1)) #
  @show p
  # Alternative choice of p which on average results in larger values:
  p = Int(floor(sqrt(5*m*M/(12*dt*C))/π+1))
  @show p
  # initialze Gp from distribution N(0, I_M)
  Gp₁ = randn(M)


  # compute J_{ij}^p
  Jp = 1//2 .* dW .* dW'
  Ap = zero(Jp)
  aᵢ₀ = zero(dW)

  𝑎ₚ = pi^2/6 # for tail approx.

  for r=1:p
    𝑎ₚ -= 1/r^2

    sd = sqrt(dt/(2*pi*r))
    ζ = randn(m)*sd
    η = randn(m)*sd
    aᵢ₀ -= (2/sqrt(π*r))*ζ
    @show η, ζ
    for i=1:m
      for j=1:m
        Ap[i,j] += ζ[i]*η[j]-η[i]*ζ[j]
      end
    end
  end
  @show aᵢ₀, Ap
  for i=1:m
    for j=1:m
      Jp[i,j] += Ap[i,j] - 1//2*(aᵢ₀[j]*dW[i]-dW[j]*aᵢ₀[i])
    end
  end
  @show Jp



  return Jp
end

seed = 10
Random.seed!(seed)

m = 5
W = WienerProcess(0.0,zeros(m),nothing)

dt = 0.1
calculate_step!(W,dt,nothing,nothing)

for i in 1:10
  accept_step!(W,dt,nothing,nothing)
end

true_diag = 1//2 .* W.dW .* W.dW

Wikdiag = StochasticDiffEq.WikJDiagonal_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikdiag, 1.0)

@test Wikdiag.WikJ == true_diag
@show true_diag

true_commute = 1//2 .* W.dW .* W.dW'

Wikcommute = StochasticDiffEq.WikJCommute_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikcommute, 1.0)

@test Wikcommute.WikJ == true_commute
@show true_commute


Random.seed!(seed)
Wikgeneral = StochasticDiffEq.WikJGeneral_iip(W.dW)
StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikgeneral, 1.0)

Random.seed!(seed)
Wikgeneraloop  = StochasticDiffEq.WikJGeneral_oop(W.dW)
@test StochasticDiffEq.get_iterated_I!(dt, W.dW, Wikgeneraloop, 1.0) == Wikgeneral.WikJ

Random.seed!(seed)
true_noncom = true_general_function(dt, W.dW, 1.0, m)

@test Wikgeneral.WikJ == true_noncom
@show true_noncom



Random.seed!(seed)
zeta=randn(m)*sqrt(dt/(2*pi))
eta=randn(m)*sqrt(dt/(2*pi))



K = zeros(m,m)
for i=1:m
  K[i,i]=0.5*W.dW[i]^2
end
@show K
K[:,1]-= zeta*(2/sqrt(pi))
@show K
for i=2:m
  for j=1:i-1
    K[i,j]+=zeta[i]*eta[j]-eta[i]*zeta[j]
  end
end
@show K
for i=2:m
  for j=1:i-1
    @show i,j
    tmp1=K[i,j]-0.5*(K[j,1]*W.dW[i]-K[i,1]*W.dW[j])
    tmp2=0.5*W.dW[i]*W.dW[j]
    K[i,j]=tmp2+tmp1
    K[j,i]=tmp2-tmp1
  end
end
@show K
