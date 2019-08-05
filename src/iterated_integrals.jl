# This function calculates WikJ, a mxm Array for a m dimensional general noise problem,
# which is a approximation to the second order iterated integrals
# this is implementation of the section 4 of the paper doi:10.1016/j.cam.2006.05.037
abstract type AbstractWikJ end
abstract type AbstractWikJDiagonal <: AbstractWikJ end
abstract type AbstractWikJCommute <: AbstractWikJ end
abstract type AbstractWikJGeneral <: AbstractWikJ end

struct WikJDiagonal_oop{WikJType} <: AbstractWikJDiagonal
    WikJ::WikJType
end

mutable struct WikJDiagonal_iip{WikJType} <: AbstractWikJDiagonal
    WikJ::WikJType
end

struct WikJCommute_oop{WikJType} <: AbstractWikJCommute
    WikJ::WikJType
end

mutable struct WikJCommute_iip{WikJType} <: AbstractWikJCommute
    WikJ::WikJType
end

struct WikJGeneral_oop{rateNoiseElTypeNoUnits, WikJType} <: AbstractWikJGeneral
    WikJ::WikJType
    m_seq::Array{Int}
end

mutable struct WikJGeneral_iip{rateNoiseElTypeNoUnits, WikJType} <: AbstractWikJGeneral
    WikJ::WikJType
    WikJ2::WikJType
    WikJ3::WikJType
    m_seq::Array{Int}
    vec_ζ::Vector{eltype(rateNoiseElTypeNoUnits)}
    vec_η::Vector{eltype(rateNoiseElTypeNoUnits)}
    Gp1::Vector{eltype(rateNoiseElTypeNoUnits)}
    Gp2::Vector{eltype(rateNoiseElTypeNoUnits)}
end

function fill_WikJDiagonal_oop(ΔW)
    WikJ = false .* ΔW .* ΔW
    WikJDiagonal_oop{typeof(WikJ)}(WikJ)
end

function fill_WikJDiagonal_iip(ΔW)
    WikJ = false .* ΔW .* ΔW
    WikJDiagonal_iip{typeof(WikJ)}(WikJ)
end

function fill_WikJCommute_oop(ΔW)
    WikJ = false .* ΔW .* ΔW'
    WikJCommute_oop{typeof(WikJ)}(WikJ)
end

function fill_WikJCommute_iip(ΔW)
    WikJ = false .* ΔW .* ΔW'
    WikJCommute_iip{typeof(WikJ)}(WikJ)
end

function fill_WikJGeneral_oop(ΔW)
    WikJ = false .* ΔW .* ΔW'
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
    WikJCommute_oop{eltype(ΔW), typeof(WikJ)}(WikJ, m_seq)
end

function fill_WikJGeneral_iip(ΔW)
    WikJ = false .* ΔW .* ΔW'
    WikJ2 = similar(WikJ)
    WikJ3 = similar(WikJ)
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
    vec_η = similar(vec_ζ)
    Gp1 = zeros(M)
    Gp2 = similar(Gp1)
    WikJCommute_iip{eltype(ΔW), typeof(WikJ)}(WikJ, WikJ2, WikJ3, m_seq, vec_ζ, vec_η, Gp1, Gp2)
end

function get_iterated_I!(dW, Wik::WikJDiagonal_oop)
    @unpack WikJ = Wik
    WikJ = 1//2 .* dW .* dW
    WikJ
end

function get_iterated_I!(dW, Wik::WikJDiagonal_iip)
    @unpack WikJ = Wik
    if typeof(dW) <: Number
        Wik.WikJ = 1//2 .* dW .^ 2
    else
        @.. WikJ = 1//2*dW^2
    end
    return nothing
end

function get_iterated_I!(dW, Wik::WikJCommute_oop)
    @unpack WikJ = Wik
    WikJ = 1//2 .* vec(dW) .* vec(dW)'
    WikJ
end

function get_iterated_I!(dW, Wik::WikJCommute_iip)
    @unpack WikJ = Wik
    mul!(WikJ,vec(dW),vec(dW)')
    @.. WikJ *= 1//2
    return nothing
end

function get_iterated_I!(dW, Wik::WikJGeneral_oop)
    @unpack WikJ, m_seq = Wik

    m      = length(dW)
    M      = m*(m-1)/2
    sum_dW² = dW'*dW

    WikJ = dW*dW'
    Gp1 = randn(M)
    α = sqrt(1 + sum_dW²/dt)
    Gp2 = Gp1/(sqrt(2)*(1+α)*dt)

    #operator (Iₘ² - Pₘ)Kₘᵀ
    for i in 1:M
        WikJ2[m_seq[i,1], m_seq[i,2]] = Gp2[i]
        WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp2[i]
    end

    #operator (Iₘ X W*Wᵀ)
    WikJ2 = WikJ*WikJ2

    #operator Kₘ(Iₘ² - Pₘ)
    WikJ2 = WikJ2 - WikJ2'
    for i in 1:M
        Gp2[i] = WikJ2[m_seq[i,1], m_seq[i,2]]
    end
    Gp = Gp/sqrt(2) + Gp2

    #operator (Iₘ² - Pₘ)Kₘᵀ
    for i in 1:M
        WikJ2[m_seq[i,1], m_seq[i,2]] = Gp[i]
        WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp[i]
    end

    WikJ *= 1//2
    a2ₚ = (π^2)/6
    p = Int(floor((1/π)*sqrt(M/(24*dt))*sqrt(m + 4*sum_dW²/dt) + 1))
    for i in 1:p
        a2ₚ -= (1/i^2)
        var = sqrt(dt/(2*π*i))
        vec_ζ = randn(m)*var
        vec_η = randn(m)*var
        WikJ += (vec_ζ*vec_η' - vec_η*vec_ζ')
        Aₚ -= (2/sqrt(π*i))*vec_ζ
    end

    WikJ -= 1//2*(dW*Aₚ' - Aₚ*dW')
    WikJ += (sqrt(a2ₚ)*dt/π)*WikJ2

    WikJ
end

function get_iterated_I!(dW, Wik::WikJGeneral_iip)
    @unpack WikJ, WikJ2, WikJ3, m_seq, vec_ζ, vec_η, Gp1, Gp2 = Wik

    m      = length(dW)
    M      = m*(m-1)/2

    sum_dW² = zero(eltype(dW))
    mul!(sum_dW²,dW', dW)

    @.. Gp1 = randn(M)
    α = sqrt(1 + sum_dW²/dt)
    @.. Gp2 = Gp1/(sqrt(2)*(1+α)*dt)

    #operator (Iₘ² - Pₘ)Kₘᵀ
    for i in 1:M
        WikJ2[m_seq[i,1], m_seq[i,2]] = Gp2[i]
        WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp2[i]
    end

    #operator (Iₘ X W*Wᵀ)
    mul!(WikJ,dW,dW')
    mul!(WikJ3,WikJ,WikJ2)

    #operator Kₘ(Iₘ² - Pₘ)
    @.. WikJ2 = WikJ3 - WikJ3'
    for i in 1:M
        Gp2[i] = WikJ2[m_seq[i,1], m_seq[i,2]]
    end
    Gp1 = Gp1/sqrt(2) + Gp2

    #operator (Iₘ² - Pₘ)Kₘᵀ
    for i in 1:M
        WikJ2[m_seq[i,1], m_seq[i,2]] = Gp1[i]
        WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp1[i]
    end

    @.. WikJ *= 1//2
    a2ₚ = (π^2)/6
    p = Int(floor((1/π)*sqrt(M/(24*dt))*sqrt(m + 4*sum_dW²/dt) + 1))
    for i in 1:p
        a2ₚ -= (1/i^2)
        var = sqrt(dt/(2*π*i))
        @.. vec_ζ = randn(m)*var
        @.. vec_η = randn(m)*var
        mul!(WikJ3, vec_ζ, vec_η')
        @.. WikJ += WikJ3 - WikJ3'
        @.. Aₚ -= (2/sqrt(π*i))*vec_ζ
    end
    mul!(WikJ3, dW, Aₚ')
    @.. WikJ -= 1//2*(WikJ3 - WikJ3')
    @.. WikJ += (sqrt(a2ₚ)*dt/π)*WikJ2

    return nothing
end
