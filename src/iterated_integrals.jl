abstract type AbstractWikJ end
abstract type AbstractWikJDiagonal <: AbstractWikJ end
abstract type AbstractWikJCommute <: AbstractWikJ end
abstract type AbstractWikJGeneral <: AbstractWikJ end

struct WikJDiagonal_oop <: AbstractWikJDiagonal end

mutable struct WikJDiagonal_iip{WikJType} <: AbstractWikJDiagonal
    WikJ::WikJType
end

struct WikJCommute_oop <: AbstractWikJCommute end

mutable struct WikJCommute_iip{WikJType} <: AbstractWikJCommute
    WikJ::WikJType
end

struct WikJGeneral_oop <: AbstractWikJGeneral
    m_seq::Array{Int}
end

mutable struct WikJGeneral_iip{rateNoiseElTypeNoUnits, WikJType} <: AbstractWikJGeneral
    WikJ::WikJType
    WikJ2::WikJType
    WikJ3::WikJType
    m_seq::Array{Int}
    vec_ζ::Vector{eltype(rateNoiseElTypeNoUnits)}
    vec_η::Vector{eltype(rateNoiseElTypeNoUnits)}
    Gp₁::Vector{eltype(rateNoiseElTypeNoUnits)}
    Gp₂::Vector{eltype(rateNoiseElTypeNoUnits)}
    Aᵢ::Vector{eltype(rateNoiseElTypeNoUnits)}
end

function WikJDiagonal_oop(ΔW)
    new()
end

function WikJDiagonal_iip(ΔW)
    WikJ = false .* ΔW .* ΔW
    new{typeof(WikJ)}(WikJ)
end

function WikJCommute_oop(ΔW)
    new()
end

function WikJCommute_iip(ΔW)
    WikJ = false .* ΔW .* ΔW'
    new{typeof(WikJ)}(WikJ)
end

function WikJGeneral_oop(ΔW)
    m = length(ΔW)
    M = m*(m-1)/2
    m_seq = Array{Int}(undef, M, 2)
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

function WikJGeneral_iip(ΔW)
    WikJ = false .* ΔW .* ΔW'
    WikJ2 = false .* ΔW .* ΔW'
    WikJ3 = false .* ΔW .* ΔW'
    m = length(ΔW)
    M = m*(m-1)/2
    m_seq = Array{Int}(undef, M, 2)
    k = 1
    for i in 1:length(ΔW)
      for j in i+1:length(ΔW)
        m_seq[k,1] = i
        m_seq[k,2] = j
        k += 1
      end
    end
    vec_ζ = false .* vec(ΔW)
    vec_η = false .* vec(ΔW)
    Gp₁ = false .* Array{eltype(ΔW)}(undef, M)
    Gp₂ = false .* Array{eltype(ΔW)}(undef, M)
    Aᵢ = false .* vec(ΔW)
    new{eltype(ΔW), typeof(WikJ)}(WikJ, WikJ2, WikJ3, m_seq, vec_ζ, vec_η, Gp₁, Gp₂, Aᵢ)
end

fill_WikJ(ΔW,::Val{1},::Val{false}) = WikJDiagonal_oop(ΔW)
fill_WikJ(ΔW,::Val{1},::Val{true}) = WikJDiagonal_iip(ΔW)
fill_WikJ(ΔW,::Val{2},::Val{false}) = WikJCommute_oop(ΔW)
fill_WikJ(ΔW,::Val{2},::Val{true}) = WikJCommute_iip(ΔW)
fill_WikJ(ΔW,::Val{3},::Val{false}) = WikJGeneral_oop(ΔW)
fill_WikJ(ΔW,::Val{3},::Val{true}) = WikJGeneral_iip(ΔW)

"""

    get_iterated_I!(dW, Wik::AbstractWikJ, C=1)

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

```math
    Oper2(Gp₁) = (√(Σ∞)/2 Gp₁) = (Iₘ/√2 + F/(√2 * (1 + α) * Δt))(Gp₁)
```
then,

```math
    (√(Σ∞)/2 Gp₁) = Iₘ*Gp₁/√2 + F(Gp₁/(√2*(1+α)*Δt))
    (√(Σ∞)/2 Gp₁) = Gp₁/√2 + F(Gp₁/(√2*(1+α)*Δt))
```

initially we have

    Gp₂ = Gp₁/(sqrt(2)*(1+α)*dt)

thus applying operator `F` on `Gp₂`. And hence we have

    Gp₂ = F(Gp₂)
    Gp₁ = (Gp₁/√2) + Gp₂

```math
    𝒜ᵖ = (Iₘ² - Pₘ)Kₘᵀ Δt/π √(𝑎ₚ) Oper2(Gp₁)
    𝒜ᵖ = √(𝑎ₚ)*Δt/π * (Iₘ² - Pₘ)Kₘᵀ(Oper2(Gp₁))
```
In the code we have

```math
    WikJ2 = (Iₘ² - Pₘ)Kₘᵀ(Oper2(Gp₁))
```

"""

function get_iterated_I!(dW, Wik::WikJDiagonal_oop, C=1)
    WikJ = 1//2 .* dW .* dW
    WikJ
end

function get_iterated_I!(dW, Wik::WikJDiagonal_iip, C=1)
    @unpack WikJ = Wik
    if typeof(dW) <: Number
        Wik.WikJ = 1//2 .* dW .^ 2
    else
        @.. WikJ = 1//2*dW^2
    end
    return nothing
end

function get_iterated_I!(dW, Wik::WikJCommute_oop, C=1)
    WikJ = 1//2 .* vec(dW) .* vec(dW)'
    WikJ
end

function get_iterated_I!(dW, Wik::WikJCommute_iip, C=1)
    @unpack WikJ = Wik
    mul!(WikJ,vec(dW),vec(dW)')
    @.. WikJ *= 1//2
    return nothing
end

function get_iterated_I!(dW, Wik::WikJGeneral_oop, C=1)
    @unpack m_seq = Wik
    m      = length(dW)
    M      = m*(m-1)/2
    sum_dW² = dW'*dW

    WikJ = dW*dW'
    Gp₁ = randn(M)
    α = sqrt(1 + sum_dW²/dt)
    Gp₂ = Gp₁/(sqrt(2)*(1+α)*dt)

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
    Gp₁ = Gp₁/sqrt(2) + Gp₂

    #operator (Iₘ² - Pₘ)Kₘᵀ
    for i in 1:M
        WikJ2[m_seq[i,1], m_seq[i,2]] = Gp[i]
        WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp[i]
    end

    WikJ *= 1//2
    𝑎ₚ = (π^2)/6
    p = Int(floor((1/(C*π))*sqrt(M/(24*dt))*sqrt(m + 4*sum_dW²/dt) + 1))
    Aᵢ = false .* vec(dW)   # Aᵢ is vector of aᵢ₀
    for r in 1:p
        𝑎ₚ -= (1/r^2)
        var = sqrt(dt/(2*π*r))
        vec_ζ = randn(m)*var
        vec_η = randn(m)*var
        WikJ += (vec_ζ*vec_η' - vec_η*vec_ζ')
        Aᵢ -= (2/sqrt(π*r))*vec_ζ
    end

    WikJ -= 1//2*(dW*Aᵢ' - Aᵢ*dW')
    WikJ += (sqrt(𝑎ₚ)*dt/π)*WikJ2
    WikJ
end

function get_iterated_I!(dW, Wik::WikJGeneral_iip, C=1)
    @unpack WikJ, WikJ2, WikJ3, m_seq, vec_ζ, vec_η, Gp₁, Gp₂, Aᵢ = Wik

    m      = length(dW)
    M      = m*(m-1)/2

    sum_dW² = zero(eltype(dW))
    mul!(sum_dW²,dW', dW)

    @.. Gp₁ = randn(M)
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
    p = Int(floor((1/(C*π))*sqrt(M/(24*dt))*sqrt(m + 4*sum_dW²/dt) + 1))
    @.. Aᵢ = false .* vec(dW)    # Aᵢ is vector of aᵢ₀
    for r in 1:p
        𝑎ₚ -= (1/r^2)
        var = sqrt(dt/(2*π*r))
        @.. vec_ζ = randn(m)*var
        @.. vec_η = randn(m)*var
        mul!(WikJ3, vec_ζ, vec_η')
        @.. WikJ += WikJ3 - WikJ3'
        @.. Aᵢ -= (2/sqrt(π*r))*vec_ζ
    end
    mul!(WikJ3, dW, Aᵢ')
    @.. WikJ -= 1//2*(WikJ3 - WikJ3')
    @.. WikJ += (sqrt(𝑎ₚ)*dt/π)*WikJ2
    return nothing
end
