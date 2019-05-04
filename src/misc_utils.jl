struct DiffEqNLSolveTag end

struct DiffCache{T<:AbstractArray, S<:AbstractArray}
    du::T
    dual_du::S
end

Base.@pure function DiffCache(T, size, ::Type{Val{chunk_size}}) where chunk_size
    DiffCache(fill(zero(T), size...), fill(zero(Dual{typeof(ForwardDiff.Tag(DiffEqNLSolveTag(),T)),T,chunk_size}), size...))
end

Base.@pure DiffCache(u::AbstractArray) = DiffCache(eltype(u),size(u),Val{ForwardDiff.pickchunksize(length(u))})
Base.@pure DiffCache(u::AbstractArray,nlsolve) = DiffCache(eltype(u),size(u),Val{get_chunksize(nlsolve)})
Base.@pure DiffCache(u::AbstractArray,T::Type{Val{CS}}) where {CS} = DiffCache(eltype(u),size(u),T)

get_du(dc::DiffCache, ::Type{T}) where {T<:Dual} = dc.dual_du
get_du(dc::DiffCache, T) = dc.du

# Default nlsolve behavior, should move to DiffEqDiffTools.jl

Base.@pure determine_chunksize(u,alg::DEAlgorithm) = determine_chunksize(u,get_chunksize(alg))
Base.@pure function determine_chunksize(u,CS)
  if CS != 0
    return CS
  else
    return ForwardDiff.pickchunksize(length(u))
  end
end

struct NLSOLVEJL_SETUP{CS,AD} end
Base.@pure NLSOLVEJL_SETUP(;chunk_size=0,autodiff=true) = NLSOLVEJL_SETUP{chunk_size,autodiff}()
(::NLSOLVEJL_SETUP)(f,u0; kwargs...) = (res=NLsolve.nlsolve(f,u0; kwargs...); res.zero)
function (p::NLSOLVEJL_SETUP{CS,AD})(::Type{Val{:init}},f,u0_prototype) where {CS,AD}
  AD ? autodiff = :forward : autodiff = :central
  OnceDifferentiable(f, u0_prototype, u0_prototype, autodiff,
                     ForwardDiff.Chunk(determine_chunksize(u0_prototype,CS)))
end

get_chunksize(x) = 0
get_chunksize(x::NLSOLVEJL_SETUP{CS,AD}) where {CS,AD} = CS

"""
    calculate_residuals(ũ, u₀, u₁, α, ρ, scalarnorm)

Return element-wise residuals
```math
\\frac{ũ}{α+\\max{scalarnorm(u₀),scalarnorm(u₁)}*ρ}.
```
"""
@inline @muladd function calculate_residuals(ũ::Number, u₀::Number, u₁::Number, α::Real,
                                             ρ::Real, scalarnorm, t)
    ũ / (α + max(scalarnorm(u₀,t), scalarnorm(u₁,t)) * ρ)
end

"""
    calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm)

Return element-wise residuals
```math
\\frac{δ E₁ + E₂}{α+\\max{scalarnorm(u₀),scalarnorm(u₁)}*ρ}.
```
"""
@inline @muladd function calculate_residuals(E₁::Number, E₂::Number, u₀::Number, u₁::Number,
                                             α::Real, ρ::Real, δ::Number, scalarnorm, t)
    (δ * E₁ + E₂) / (α + max(scalarnorm(u₀,t), scalarnorm(u₁,t)) * ρ)
end

"""
    calculate_residuals!(out, ũ, u₀, u₁, α, ρ, scalarnorm)

Same as [`calculate_residuals`](@ref) but save result in `out`.
"""
@inline function calculate_residuals!(out, ũ, u₀, u₁, α, ρ, scalarnorm, t)
  @.. out = calculate_residuals(ũ, u₀, u₁, α, ρ, scalarnorm, t)
  out
end

@inline function calculate_residuals!(out::Array{<:Number}, ũ::Array{<:Number},
                                      u₀::Array{<:Number}, u₁::Array{<:Number}, α::Real,
                                      ρ::Real, scalarnorm, t)
  @tight_loop_macros @inbounds for i in eachindex(out)
    out[i] = calculate_residuals(ũ[i], u₀[i], u₁[i], α, ρ, scalarnorm, t)
  end
  out
end

@inline function calculate_residuals(ũ, u₀, u₁, α, ρ, scalarnorm, t)
  @.. calculate_residuals(ũ, u₀, u₁, α, ρ, scalarnorm, ts)
end

@inline function calculate_residuals(ũ::Array{<:Number}, u₀::Array{<:Number},
                                     u₁::Array{<:Number}, α::Real, ρ::Real, scalarnorm, t)
  out = similar(ũ)
  calculate_residuals!(out, ũ, u₀, u₁, α, ρ, scalarnorm, t)
  out
end

"""
    calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm)

Same as [`calculate_residuals`](@ref) but save result in `out`.
"""
@inline function calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
  @.. out = calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
  out
end

@inline function calculate_residuals!(out::Array{<:Number}, E₁::Array{<:Number},
                                      E₂::Array{<:Number}, u₀::Array{<:Number},
                                      u₁::Array{<:Number}, α::Real, ρ::Real, δ::Number,
                                      scalarnorm, t)
  @tight_loop_macros @inbounds for i in eachindex(out)
      out[i] = calculate_residuals(E₁[i], E₂[i], u₀[i], u₁[i], α, ρ, δ, scalarnorm, t)
  end
  out
end

@inline function calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
  @.. calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
end

@inline function calculate_residuals(E₁::Array{<:Number}, E₂::Array{<:Number},
                                     u₀::Array{<:Number}, u₁::Array{<:Number}, α::Real,
                                     ρ::Real, δ::Number, scalarnorm, t)
    out = similar(u₀)
    calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
    out
end

macro cache(expr)
  name = expr.args[2].args[1].args[1]
  fields = expr.args[3].args[2:2:end]
  cache_vars = Expr[]
  rand_vars = Expr[]
  jac_vars = Pair{Symbol,Expr}[]
  ratenoise_vars = Expr[]
  for x in fields
    if x.args[2] == :uType || x.args[2] == :rateType ||
       x.args[2] == :kType || x.args[2] == :uNoUnitsType #|| x.args[2] == :possibleRateType
      push!(cache_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :JCType
      push!(cache_vars,:(c.$(x.args[1]).duals...))
    elseif x.args[2] == :GCType
      push!(cache_vars,:(c.$(x.args[1]).duals))
    elseif x.args[2] == :DiffCacheType
      push!(cache_vars,:(c.$(x.args[1]).du))
      push!(cache_vars,:(c.$(x.args[1]).dual_du))
    elseif x.args[2] == :JType || x.args[2] == :WType
      push!(jac_vars,x.args[1] => :(c.$(x.args[1])))
    elseif x.args[2] == :randType
      push!(rand_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :rateNoiseType || x.args[2] == :rateNoiseCollectionType
      # Should be a pair for handling non-diagonal
      push!(ratenoise_vars,:(c.$(x.args[1])))
    end
  end
  quote
    $expr
    $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
    $(esc(:jac_iter))($(esc(:c))::$name) = tuple($(jac_vars...))
    $(esc(:rand_cache))($(esc(:c))::$name) = tuple($(rand_vars...))
    $(esc(:ratenoise_cache))($(esc(:c))::$name) = tuple($(ratenoise_vars...))
  end
end

_reshape(v, siz) = reshape(v, siz)
_reshape(v::Number, siz) = v
_reshape(v::AbstractVector, siz) = v
_vec(v) = vec(v)
_vec(v::Number) = v
_vec(v::AbstractVector) = v
