function calc_J!(integrator, cache::StochasticDiffEqConstantCache, is_compos)
  @unpack t,dt,uprev,u,f,p = integrator
  if has_jac(f)
    J = f.jac(uprev, p, t)
  else
    error("Jacobian wrapper for constant caches not yet implemented") #TODO
  end
  return J
end
function calc_J!(integrator, cache::StochasticDiffEqConstantCache, is_compos)
  @unpack t,dt,uprev,u,f,p = integrator
  J = cache.J
  if has_jac(f)
    f.jac(J, uprev, p, t)
  else
    @unpack du1,uf,jac_config = cache
    uf.t = t
    uf.p = p
    jacobian!(J, uf, uprev, du1, integrator, jac_config)
    if is_compos
      integrator.eigen_est = norm(J, Inf)
    end
  end
end

"""
    WOperator(mass_matrix,gamma,J[;transform=false])

A linear operator that represents the W matrix of an ODEProblem, defined as

```math
W = MM - \\gamma J
```

or, if `transform=true`:

```math
W = \\frac{1}{\\gamma}MM - J
```

where `MM` is the mass matrix (a regular `AbstractMatrix` or a `UniformScaling`),
`γ` is a real number proportional to the time step, and `J` is the Jacobian
operator (must be a `AbstractDiffEqLinearOperator`). A `WOperator` can also be
constructed using a `*DEFunction` directly as

    WOperator(f,gamma[;transform=false])

`f` needs to have a jacobian and `jac_prototype`, but the prototype does not need
to be a diffeq operator --- it will automatically be converted to one.

`WOperator` supports lazy `*` and `mul!` operations, the latter utilizing an
internal cache (can be specified in the constructor; default to regular `Vector`).
It supports all of `AbstractDiffEqLinearOperator`'s interface.
"""
mutable struct WOperator{T,
  MType <: Union{UniformScaling,AbstractMatrix},
  GType <: Real,
  JType <: DiffEqBase.AbstractDiffEqLinearOperator{T}
  } <: DiffEqBase.AbstractDiffEqLinearOperator{T}
  mass_matrix::MType
  gamma::GType
  J::JType
  transform::Bool       # true => W = mm/gamma - J; false => W = mm - gamma*J
  _func_cache           # cache used in `mul!`
  _concrete_form        # non-lazy form (matrix/number) of the operator
  WOperator(mass_matrix, gamma, J; transform=false) = new{eltype(J),typeof(mass_matrix),
    typeof(gamma),typeof(J)}(mass_matrix,gamma,J,transform,nothing,nothing)
end
function WOperator(f::DiffEqBase.AbstractODEFunction, gamma; transform=false)
  @assert DiffEqBase.has_jac(f) "f needs to have an associated jacobian"
  if isa(f, Union{SplitFunction, DynamicalODEFunction})
    error("WOperator does not support $(typeof(f)) yet")
  end
  # Convert mass matrix, if needed
  mass_matrix = f.mass_matrix
  if !isa(mass_matrix, Union{AbstractMatrix,UniformScaling})
    mass_matrix = convert(AbstractMatrix, mass_matrix)
  end
  # Convert jacobian, if needed
  J = deepcopy(f.jac_prototype)
  if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
    J = DiffEqArrayOperator(J; update_func=f.jac)
  end
  return WOperator(mass_matrix, gamma, J; transform=transform)
end

set_gamma!(W::WOperator, gamma) = (W.gamma = gamma; W)
DiffEqBase.update_coefficients!(W::WOperator,u,p,t) = (update_coefficients!(W.J,u,p,t); W)
function Base.convert(::Type{AbstractMatrix}, W::WOperator)
  if W._concrete_form == nothing
    # Allocating
    if W.transform
      W._concrete_form = W.mass_matrix / W.gamma - convert(AbstractMatrix,W.J)
    else
      W._concrete_form = W.mass_matrix - W.gamma * convert(AbstractMatrix,W.J)
    end
  else
    # Non-allocating
    if W.transform
      rmul!(copyto!(W._concrete_form, W.mass_matrix), 1/W.gamma)
      axpy!(-1, convert(AbstractMatrix,W.J), W._concrete_form)
    else
      copyto!(W._concrete_form, W.mass_matrix)
      axpy!(-W.gamma, convert(AbstractMatrix,W.J), W._concrete_form)
    end
  end
  W._concrete_form
end
function Base.convert(::Type{Number}, W::WOperator)
  if W.transform
    W._concrete_form = W.mass_matrix / W.gamma - convert(Number,W.J)
  else
    W._concrete_form = W.mass_matrix - W.gamma * convert(Number,W.J)
  end
  W._concrete_form
end
Base.size(W::WOperator, args...) = size(W.J, args...)
function Base.getindex(W::WOperator, i::Int)
  if W.transform
    W.mass_matrix[i] / W.gamma - W.J[i]
  else
    W.mass_matrix[i] - W.gamma * W.J[i]
  end
end
function Base.getindex(W::WOperator, I::Vararg{Int,N}) where {N}
  if W.transform
    W.mass_matrix[I...] / W.gamma - W.J[I...]
  else
    W.mass_matrix[I...] - W.gamma * W.J[I...]
  end
end
function Base.:*(W::WOperator, x::Union{AbstractVecOrMat,Number})
  if W.transform
    (W.mass_matrix*x) / W.gamma - W.J*x
  else
    W.mass_matrix*x - W.gamma * (W.J*x)
  end
end
function Base.:\(W::WOperator, x::Union{AbstractVecOrMat,Number})
  if size(W) == () # scalar operator
    convert(Number,W) \ x
  else
    convert(AbstractMatrix,W) \ x
  end
end
function LinearAlgebra.mul!(Y::AbstractVecOrMat, W::WOperator, B::AbstractVecOrMat)
  if W._func_cache == nothing
    # Allocate cache only if needed
    W._func_cache = Vector{eltype(W)}(undef, size(Y, 1))
  end
  if W.transform
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      a = W.mass_matrix.λ / W.gamma
      @. Y = a * B
    else
      mul!(Y, W.mass_matrix, B)
      lmul!(1/W.gamma, Y)
    end
    # Compute J * B and subtract
    mul!(W._func_cache, W.J, B)
    Y .-= W._func_cache
  else
    # Compute mass_matrix * B
    if isa(W.mass_matrix, UniformScaling)
      @. Y = W.mass_matrix.λ * B
    else
      mul!(Y, W.mass_matrix, B)
    end
    # Compute J * B
    mul!(W._func_cache, W.J, B)
    # Subtract result
    axpy!(-W.gamma, W._func_cache, Y)
  end
end

function calc_W!(integrator, cache::StochasticDiffEqMutableCache, γdt, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p, alg = integrator
    @unpack J,W,jac_config = cache
    is_compos = is_composite(alg)
    alg = unwrap_alg(integrator, true)
    mass_matrix = integrator.f.mass_matrix

    new_W = true
    if has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && f(Val{:invW},W,uprev,p,γdt,t) # W == inverse W
      is_compos && calc_J!(integrator, cache, true)
    else
      # skip calculation of J if step is repeated
      if repeat_step || (!integrator.last_stepfail && cache.newton_iters == 1 && cache.ηold < alg.new_jac_conv_bound)
        new_jac = false
      else # Compute a new Jacobian
        new_jac = true
        calc_J!(integrator, cache, is_compos)
      end
      # skip calculation of W if step is repeated
      if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t)))
        for j in 1:length(u), i in 1:length(u)
            @inbounds W[i,j] = mass_matrix[i,j]-γdt*J[i,j]
        end
      else
        new_W = false
      end
    end
    return new_W
  end
end

function calc_W!(integrator, cache::StochasticDiffEqConstantCache, γdt, repeat_step)
  uprev = integrator.uprev
  uf = cache.uf
  is_compos = is_composite(integrator.alg)
  if typeof(uprev) <: AbstractArray
    J = ForwardDiff.jacobian(uf,uprev)
    is_compos && ( integrator.eigen_est = norm(J, Inf) )
    W = I - γdt*J
  else
    J = ForwardDiff.derivative(uf,uprev)
    W = 1 - γdt*J
    is_compos && ( integrator.eigen_est = J )
  end
  J, W
end
