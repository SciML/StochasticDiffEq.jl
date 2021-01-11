"""
All stochastic iterated integrals are written in Stratonovich sense as indicated
by the J.
"""

abstract type AbstractJ end

"""
Kloeden, P. E., Platen, E., & Wright, I. W., The approximation of multiple stochastic integrals.
Stochastic analysis and applications, 10(4), pp. 431-441 (1992).
DOI: 10.1080/07362999208809281
"""
abstract type AbstractKPWJ <: AbstractJ end

struct KPWJ_oop <: AbstractKPWJ end

mutable struct KPWJ_iip{WikAType, WikJType, randType} <: AbstractKPWJ
  WikA::WikAType
  WikJ::WikJType
  μ::randType
  ζ::randType
  η::randType
  kronprod1::WikAType
  kronprod2::WikAType
end

function KPWJ_iip(ΔW)
  WikJ = false .* vec(ΔW) .* vec(ΔW)'
  WikA = false .* vec(ΔW) .* vec(ΔW)'

  μ = zero(ΔW)
  ζ = zero(ΔW)
  η = zero(ΔW)

  kronprod1 = zero(WikA)
  kronprod2 = zero(WikA)

  KPWJ_iip(WikA,WikJ,μ,ζ,η,kronprod1,kronprod2)
end


function get_iterated_I(dt, dW, dZ, Wik::KPWJ_oop, p=nothing, C=1, γ=1//1)
  m  = length(dW)

  WikJ = vec(dW) .* vec(dW)'
  WikA = zero(WikJ)

  p==nothing && (p = Int(floor(C*dt^(1//1-2//1*γ)) + 1))

  # Eq.(20)
  ρp = (π^2)/6
  for k in 1:p
    ρp -= (1/k^2)
  end
  ρp = ρp/(2*π^2)

  if dZ==nothing
    μ = randn(eltype(dW),m)
  else
    μ = convert(eltype(dW),sqrt(1/dt))*dZ[1:m]
  end

  for k in 1:p
    if dZ==nothing
      ζ = randn(eltype(dW),m)
      η = randn(eltype(dW),m)
    else
      ζ = convert(eltype(dW),sqrt(1/dt))*(dZ[Int(m*(2*k-1)+1):Int(m*2*k)])
      η = convert(eltype(dW),sqrt(1/dt))*(dZ[Int(m*2*k+1):Int(m*(2*k+1))])
    end

    kronprod1 = vec(ζ) .* vec(η+convert(eltype(dW),sqrt(2/dt))*dW)'
    kronprod2 = vec(η+convert(eltype(dW),sqrt(2/dt))*dW) .* vec(ζ)'

    WikA += (kronprod1 - kronprod2)/k
  end

  WikA *=  dt/(2*π)

  kronprod1 = vec(μ) .* vec(dW)'
  kronprod2 = vec(dW).* vec(μ)'

  WikJ = WikJ/2 + WikA + sqrt(dt*ρp)*(kronprod1-kronprod2)  # Ito: (UniformScaling(dt))/2
  return WikJ
end


function get_iterated_I!(dt, dW, dZ, Wik::KPWJ_iip, p=nothing, C=1, γ=1//1)
  m  = length(dW)
  @unpack WikA,WikJ,μ,ζ,η,kronprod1,kronprod2 = Wik

  mul!(WikJ,vec(dW),vec(dW)')
  fill!(WikA, zero(eltype(WikJ)))

  # Below Eq.(26): truncation
  p==nothing && (p = Int(floor(C*dt^(1//1-2//1*γ)) + 1))
  #@show p, dZ
  # Eq.(20)
  ρp = (π^2)/6
  for k in 1:p
    ρp -= (1/k^2)
  end
  ρp = ρp/(2*π^2)

  if dZ==nothing
    randn!(μ)
  else
    μ .= convert(eltype(dW),sqrt(1/dt))*@view(dZ[1:m])
  end

  # Eq. (21)
  for k in 1:p
    if dZ==nothing
      randn!(ζ)
      randn!(η)
    else
      ζ .= convert(eltype(dW),sqrt(1/dt))*@view(dZ[Int(m*(2*k-1)+1):Int(m*2*k)])
      η .= convert(eltype(dW),sqrt(1/dt))*@view(dZ[Int(m*2*k+1):Int(m*(2*k+1))])
    end
    mul!(kronprod1,vec(ζ),vec(η+convert(eltype(dW),sqrt(2/dt))*dW)')
    mul!(kronprod2,vec(η+convert(eltype(dW),sqrt(2/dt))*dW),vec(ζ)')

    @.. WikA = WikA + (kronprod1 - kronprod2)/k
  end

  @.. WikA = WikA*dt/(2*π) # check if factor of 2 is correct
  mul!(kronprod1,vec(μ),vec(dW)')
  mul!(kronprod2,vec(dW),vec(μ)')
  #@show WikA, sqrt(dt*ρp)*(kronprod1-kronprod2) , WikA + sqrt(dt*ρp)*(kronprod1-kronprod2)
  #@show sqrt(dt*ρp)
  @.. WikJ = WikJ/2 + WikA + sqrt(dt*ρp)*(kronprod1-kronprod2) # Ito: -UniformScaling(dt))/2
  return nothing
end



# Wiktorssson Approximation

abstract type AbstractWikJ <: AbstractJ end
abstract type AbstractWikJDiagonal <: AbstractWikJ end
abstract type AbstractWikJCommute <: AbstractWikJ end
abstract type AbstractWikJGeneral <: AbstractWikJ end

struct WikJDiagonal_oop <: AbstractWikJDiagonal end

mutable struct WikJDiagonal_iip{WikJType} <: AbstractWikJDiagonal
  WikJ::WikJType
  WikJDiagonal_iip(ΔW) = new{typeof(ΔW)}(false .* ΔW .* ΔW)
end

struct WikJCommute_oop <: AbstractWikJCommute end

mutable struct WikJCommute_iip{WikJType} <: AbstractWikJCommute
  WikJ::WikJType
  function WikJCommute_iip(ΔW)
    WikJ = false .* vec(ΔW) .* vec(ΔW)'
    new{typeof(WikJ)}(WikJ)
  end
end

struct WikJGeneral_oop <: AbstractWikJGeneral
  m_seq::Matrix{Int}
  function WikJGeneral_oop(ΔW)
    m = length(ΔW)
    M = div(m*(m-1),2)
    m_seq = Matrix{Int}(undef, M, 2)
    k = 1
    for i in 1:length(ΔW)
      for j in i+1:length(ΔW)
        m_seq[k,1] = i
        m_seq[k,2] = j
        k += 1
      end
    end
    new(m_seq)
  end
end

mutable struct WikJGeneral_iip{rateNoiseElTypeNoUnits, WikJType} <: AbstractWikJGeneral
  WikJ::WikJType
  WikJ2::WikJType
  WikJ3::WikJType
  m_seq::Matrix{Int}
  vec_ζ::Vector{rateNoiseElTypeNoUnits}
  vec_η::Vector{rateNoiseElTypeNoUnits}
  Gp₁::Vector{rateNoiseElTypeNoUnits}
  Gp₂::Vector{rateNoiseElTypeNoUnits}
  Aᵢ::Vector{rateNoiseElTypeNoUnits}
end

function WikJGeneral_iip(ΔW)
  WikJ = false .* vec(ΔW) .* vec(ΔW)'
  WikJ2 = false .* vec(ΔW) .* vec(ΔW)'
  WikJ3 = false .* vec(ΔW) .* vec(ΔW)'
  m = length(ΔW)
  M = div(m*(m-1),2)
  m_seq = Array{Int}(undef, M, 2)
  k = 1
  for i in 1:length(ΔW)
    for j in i+1:length(ΔW)
      m_seq[k,1] = i
      m_seq[k,2] = j
      k += 1
    end
  end
  vec_ζ = zero(ΔW)
  vec_η = zero(ΔW)
  Gp₁ = false .* Array{eltype(ΔW)}(undef, M)
  Gp₂ = false .* Array{eltype(ΔW)}(undef, M)
  Aᵢ = zero(ΔW)
  WikJGeneral_iip{eltype(ΔW), typeof(WikJ)}(WikJ, WikJ2, WikJ3, m_seq, vec_ζ, vec_η, Gp₁, Gp₂, Aᵢ)
end

function get_WikJ(ΔW,prob,alg)
  if isinplace(prob)
    if typeof(ΔW) <: Number || is_diagonal_noise(prob)
      return WikJDiagonal_iip(ΔW)
    elseif alg.ii_approx isa IICommutative
      return WikJCommute_iip(ΔW)
    else
      return KPWJ_iip(ΔW) # WikJGeneral_iip(ΔW)
    end
  else
    if typeof(ΔW) <: Number || is_diagonal_noise(prob)
      return WikJDiagonal_oop()
    elseif alg.ii_approx isa IICommutative
      return WikJCommute_oop()
    else
      return KPWJ_oop() #WikJGeneral_oop(ΔW)
    end
  end
end

function get_iterated_I(dt, dW, dZ, Wik::WikJDiagonal_oop, p=nothing, C=1, γ=1//1)
  WikJ = 1//2 .* dW .* dW
  WikJ
end

function get_iterated_I!(dt, dW, dZ, Wik::WikJDiagonal_iip, p=nothing, C=1, γ=1//1)
  @unpack WikJ = Wik
  if typeof(dW) <: Number
    Wik.WikJ = 1//2 .* dW .^ 2
  else
    @.. WikJ = 1//2*dW^2
  end
  return nothing
end

function get_iterated_I(dt, dW, dZ, Wik::WikJCommute_oop, p=nothing, C=1, γ=1//1)
  WikJ = 1//2 .* vec(dW) .* vec(dW)'
  WikJ
end

function get_iterated_I!(dt, dW, dZ, Wik::WikJCommute_iip, p=nothing, C=1, γ=1//1)
  @unpack WikJ = Wik
  mul!(WikJ,vec(dW),vec(dW)')
  @.. WikJ *= 1//2
  return nothing
end

"""

    get_iterated_I(dt, dW, Wik::WikJGeneral_oop, C=1)

This function calculates WikJ, a mxm Array for a m dimensional general noise problem, which is a approximation
to the second order iterated integrals.

WikJDiagonal and WikJCommute use the properties of respective noises to simplify the calculations.
While the calculation for General Noise case is taken from section 4 of [SDELab: A Package for solving stochastic differential
equations in MATLAB](https://doi.org/10.1016/j.cam.2006.05.037) and SDELAB2(https://github.com/tonyshardlow/SDELAB2)
which is the Implementation of SDELab in Julia.
```math
    𝒜ᵖ = (Iₘ² - Pₘ)Kₘᵀ Δt/(2π) √(𝑎ₚ) √(Σ∞) Gp₁
```

```math
    √(Σ∞) = (Σ∞ + 2αIₘ)/(√2 * (1 + α))
```

let the combined operators be,
```math
    F = Kₘ(Iₘ² - Pₘ)(Iₘ ⨂ W(Δt)W(Δt)ᵀ)(Iₘ² - Pₘ)Kₘᵀ
```

```math
    Σ∞ = 2Iₘ + (2/Δt)F
```

See the paper for further details of specific operators.
Here we've only shown in which order these are implemented in this code.

From above we can see:

```math
    Δt/(2π) √(𝑎ₚ) √(Σ∞) Gp₁ = Δt/π √(𝑎ₚ) (√(Σ∞)/2 Gp₁)
```

let,
```math
    Oper2(Gp₁) = (√(Σ∞)/2 Gp₁) = (Iₘ/√2 + F/(√2 * (1 + α) * Δt))(Gp₁)
    Oper2(Gp₁) = Iₘ*Gp₁/√2 + F(Gp₁/(√2*(1+α)*Δt))
    Oper2(Gp₁) = Gp₁/√2 + F(Gp₁/(√2*(1+α)*Δt))
```

we have,

    Gp₂ = Gp₁/(sqrt(2)*(1+α)*dt)

```math
    𝒜ᵖ = (Iₘ² - Pₘ)Kₘᵀ Δt/π √(𝑎ₚ) Oper2(Gp₁/√2 + F(Gp₂))
    𝒜ᵖ = √(𝑎ₚ)*Δt/π * (Iₘ² - Pₘ)Kₘᵀ Oper2(Gp₁/√2 + F(Gp₂))
```
In the code we have

```math
    WikJ2 = (Iₘ² - Pₘ)Kₘᵀ(Oper2(Gp₁))
```

"""
function get_iterated_I(dt, dW, dZ, Wik::WikJGeneral_oop, p=nothing, C=1, γ=1//1)
  @unpack m_seq = Wik
  m = length(dW)
  M = div(m*(m-1),2)
  sum_dW² = dot(dW,dW)

  WikJ = dW*dW'
  eltWikJ = eltype(WikJ)
  WikJ2 =zeros(eltWikJ, size(WikJ))
  Gp₁ = randn(eltWikJ,M)
  α = convert(eltWikJ,sqrt(1 + sum_dW²/dt))
  Gp₂ = Gp₁/convert(eltWikJ, sqrt(2)*(1+α)*dt)
  #operator (Iₘ² - Pₘ)Kₘᵀ
  for i in 1:M
    WikJ2[m_seq[i,1], m_seq[i,2]] = Gp₂[i]
    WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp₂[i]
  end

  #operator (Iₘ X W*Wᵀ)
  WikJ2 = WikJ*WikJ2

  #operator Kₘ(Iₘ² - Pₘ)
  WikJ2 = WikJ2 - WikJ2'
  for i in 1:M
    Gp₂[i] = WikJ2[m_seq[i,1], m_seq[i,2]]
  end
  Gp₁ = Gp₁/convert(eltWikJ,sqrt(2))+Gp₂

  #operator (Iₘ² - Pₘ)Kₘᵀ
  for i in 1:M
    WikJ2[m_seq[i,1], m_seq[i,2]] = Gp₁[i]
    WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp₁[i]
  end

  WikJ *= 1//2
  𝑎ₚ = (π^2)/6
  p==nothing && (p = Int(floor((1/(π))*sqrt(M/(12*dt*C))*sqrt(m + 4*sum_dW²/dt) + 1)))
  Aᵢ = false .* vec(dW)   # Aᵢ is vector of aᵢ₀
  for r in 1:p
    𝑎ₚ -= (1/r^2)
    var = convert(eltWikJ, sqrt(dt/(2*π*r)))
    vec_ζ = randn(eltWikJ,m)*var
    vec_η = randn(eltWikJ,m)*var
    WikJ += (vec_ζ*vec_η' - vec_η*vec_ζ')
    Aᵢ +=  convert(eltWikJ, 2/(sqrt(π*r)))*vec_ζ
  end

  WikJ -= 1//2*(dW*Aᵢ' - Aᵢ*dW')
  WikJ -= convert(eltWikJ, sqrt(𝑎ₚ)*dt/π)*WikJ2 #-= to agree with operator version from Wiktorssson
  WikJ
end

"""

    get_iterated_I!(dt, dW, Wik::WikJGeneral_iip, C=1)

This function calculates WikJ, a mxm Array for a m dimensional general noise problem, which is a approximation
to the second order iterated integrals.

WikJDiagonal and WikJCommute use the properties of respective noises to simplify the calculations.
While the calculation for General Noise case is taken from section 4 of [SDELab: A Package for solving stochastic differential
equations in MATLAB](https://doi.org/10.1016/j.cam.2006.05.037) and SDELAB2(https://github.com/tonyshardlow/SDELAB2)
which is the Implementation of SDELab in Julia.
```math
    𝒜ᵖ = (Iₘ² - Pₘ)Kₘᵀ Δt/(2π) √(𝑎ₚ) √(Σ∞) Gp₁
```

```math
    √(Σ∞) = (Σ∞ + 2αIₘ)/(√2 * (1 + α))
```

let the combined operators be,
```math
    F = Kₘ(Iₘ² - Pₘ)(Iₘ ⨂ W(Δt)W(Δt)ᵀ)(Iₘ² - Pₘ)Kₘᵀ
```

```math
    Σ∞ = 2Iₘ + (2/Δt)F
```

See the paper for further details of specific operators.
Here we've only shown in which order these are implemented in this code.

From above we can see:

```math
    Δt/(2π) √(𝑎ₚ) √(Σ∞) Gp₁ = Δt/π √(𝑎ₚ) (√(Σ∞)/2 Gp₁)
```

let,
```math
    Oper2(Gp₁) = (√(Σ∞)/2 Gp₁) = (Iₘ/√2 + F/(√2 * (1 + α) * Δt))(Gp₁)
    Oper2(Gp₁) = Iₘ*Gp₁/√2 + F(Gp₁/(√2*(1+α)*Δt))
    Oper2(Gp₁) = Gp₁/√2 + F(Gp₁/(√2*(1+α)*Δt))
```

we have,

    Gp₂ = Gp₁/(sqrt(2)*(1+α)*dt)

```math
    𝒜ᵖ = (Iₘ² - Pₘ)Kₘᵀ Δt/π √(𝑎ₚ) Oper2(Gp₁/√2 + F(Gp₂))
    𝒜ᵖ = √(𝑎ₚ)*Δt/π * (Iₘ² - Pₘ)Kₘᵀ Oper2(Gp₁/√2 + F(Gp₂))
```
In the code we have

```math
    WikJ2 = (Iₘ² - Pₘ)Kₘᵀ(Oper2(Gp₁))
```

"""
function get_iterated_I!(dt, dW, dZ, Wik::WikJGeneral_iip, p=nothing, C=1, γ=1//1)
  @unpack WikJ, WikJ2, WikJ3, m_seq, vec_ζ, vec_η, Gp₁, Gp₂, Aᵢ = Wik

  m = length(dW)
  M = div(m*(m-1),2)

  sum_dW² = dot(dW,dW) #zero(eltype(dW))
  # mul!(sum_dW²,dW', dW)

  Gp₁ .= randn(M)
  α = sqrt(1 + sum_dW²/dt)
  @.. Gp₂ = Gp₁/(sqrt(2)*(1+α)*dt)

  #operator (Iₘ² - Pₘ)Kₘᵀ
  for i in 1:M
    WikJ2[m_seq[i,1], m_seq[i,2]] = Gp₂[i]
    WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp₂[i]
  end

  #operator (Iₘ X W*Wᵀ)
  mul!(WikJ,dW,dW')
  mul!(WikJ3,WikJ,WikJ2)

  #operator Kₘ(Iₘ² - Pₘ)
  @.. WikJ2 = WikJ3 - WikJ3'
  for i in 1:M
    Gp₂[i] = WikJ2[m_seq[i,1], m_seq[i,2]]
  end
  @.. Gp₁ = Gp₁/sqrt(2) + Gp₂

  #operator (Iₘ² - Pₘ)Kₘᵀ
  for i in 1:M
    WikJ2[m_seq[i,1], m_seq[i,2]] = Gp₁[i]
    WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp₁[i]
  end

  @.. WikJ *= 1//2
  𝑎ₚ = (π^2)/6
  p==nothing && (p = Int(floor((1/(π))*sqrt(M/(12*dt*C))*sqrt(m + 4*sum_dW²/dt) + 1)))
  Aᵢ .= false .* vec(dW)    # Aᵢ is vector of aᵢ₀
  for r in 1:p
    𝑎ₚ -= (1/r^2)
    var = sqrt(dt/(2*π*r))
    vec_ζ .= randn(m) .* var
    vec_η .= randn(m) .* var
    mul!(WikJ3, vec_ζ, vec_η')
    @.. WikJ += WikJ3 - WikJ3'
    @.. Aᵢ += (2/sqrt(π*r))*vec_ζ
  end
  mul!(WikJ3, dW, Aᵢ')
  @.. WikJ = WikJ - 1//2*(WikJ3 - WikJ3')
  @.. WikJ2 *= -(sqrt(𝑎ₚ)*dt/π) # minus sign needed to agree with operator version
  #@show WikJ, (sqrt(𝑎ₚ)*dt/π)*WikJ2
  @.. WikJ = WikJ + WikJ2
  return nothing
end
