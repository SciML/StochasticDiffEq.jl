"""
All stochastic iterated integrals are written in Stratonovich sense as indicated
by the J.
"""

abstract type AbstractJ end

# Wiktorssson Approximation (for diagonal and commutative noise processes where the Levy area approximation is not required.)

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

# algs from LevyArea.jl # LevyArea.levyarea allocates random variables and then mutates these, see e.g. 
# https://github.com/stochastics-uni-luebeck/LevyArea.jl/blob/68c5cb08ab103b4dcd3178651f7a5dd9ce8c666d/src/milstein.jl#L25
function get_iterated_I(dt, dW, dZ, alg::LevyArea.AbstractIteratedIntegralAlgorithm, p=nothing, c=1, γ=1//1)
  if isnothing(p)
      ε = c*dt^(γ+1//2)
      p = terms_needed(length(dW), dt, ε, alg, MaxL2())
  end
  I = LevyArea.levyarea(dW/√dt, p, alg)
  I .= 1//2*dW.*dW' .+ dt.*I
end


# Default algorithms, keep KPWJ_oop() to have a non-mutating version
function get_WikJ(ΔW,dt,prob,alg)
  if alg.ii_approx isa IILevyArea
    if isinplace(prob)
      if typeof(ΔW) <: Number || is_diagonal_noise(prob)
        return WikJDiagonal_iip(ΔW)
      else
        # optimal_algorithm(dim, stepsize, eps=stepsize^(3/2), norm=MaxL2())
        return LevyArea.optimal_algorithm(length(ΔW), dt)
      end
    else
      if typeof(ΔW) <: Number || is_diagonal_noise(prob)
        return WikJDiagonal_oop()
      else
        return LevyArea.optimal_algorithm(length(ΔW), dt) 
      end
    end
  elseif alg.ii_approx isa IICommutative
    if isinplace(prob)
      if typeof(ΔW) <: Number || is_diagonal_noise(prob)
        return WikJDiagonal_iip(ΔW)
      else
        return WikJCommute_iip(ΔW)
      end
    else
      if typeof(ΔW) <: Number || is_diagonal_noise(prob)
        return WikJDiagonal_oop()
      else
        return WikJCommute_oop()
      end
    end
  else
    return alg.ii_approx
  end
end

# Specific Levy area alg for an SDE solver
# function StochasticDiffEq.get_WikJ(ΔW,prob,alg::SOLVER)
#  return MronRoe()
# end