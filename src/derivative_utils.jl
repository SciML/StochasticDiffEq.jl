function calc_tderivative!(integrator, cache, dtd1, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack du2,fsalfirst,dT,tf,linsolve_tmp = cache

    # Time derivative
    if !repeat_step # skip calculation if step is repeated
      if DiffEqBase.has_tgrad(f)
        f.tgrad(dT, uprev, p, t)
      else
        tf.uprev = uprev
        tf.p = p
        derivative!(dT, tf, t, du2, integrator, cache.grad_config)
      end
    end

    f(fsalfirst, uprev, p, t)
    @.. linsolve_tmp = fsalfirst + dtd1*dT
  end
end

function calc_tderivative(integrator, cache)
  @unpack t,dt,uprev,u,f,p = integrator

  # Time derivative
  if DiffEqBase.has_tgrad(f)
    dT = f.tgrad(uprev, p, t)
  else
    tf = cache.tf
    tf.u = uprev
    tf.p = p
    dT = derivative(tf, t, integrator)
  end
  dT
end

"""
    calc_J!(integrator,cache,is_compos)

Interface for calculating the jacobian.

For constant caches, a new jacobian object is returned whereas for mutable
caches `cache.J` is updated. In both cases, if `integrator.f` has a custom
jacobian update function, then it will be called for the update. Otherwise,
either ForwardDiff or finite difference will be used depending on the
`jac_config` of the cache.
"""
function calc_J(integrator, cache::StochasticDiffEqConstantCache, is_compos)
  @unpack t,dt,uprev,u,f,p = integrator
  if DiffEqBase.has_jac(f)
    J = f.jac(uprev, p, t)
  else
    J = jacobian(cache.uf,uprev,integrator)
  end
  is_compos && (integrator.eigen_est = opnorm(J, Inf))
  return J
end

function calc_J!(integrator, cache::StochasticDiffEqMutableCache, is_compos)
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
    WOperator(mass_matrix,gamma,J)

A linear operator that represents the W matrix of an ODEProblem, defined as

```math
W = MM - \\gamma J
```

where `MM` is the mass matrix (a regular `AbstractMatrix` or a `UniformScaling`),
`γ` is a real number proportional to the time step, and `J` is the Jacobian
operator (must be a `AbstractDiffEqLinearOperator`). A `WOperator` can also be
constructed using a `*DEFunction` directly as

    WOperator(f,gamma)

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
  _func_cache           # cache used in `mul!`
  _concrete_form        # non-lazy form (matrix/number) of the operator
  WOperator(mass_matrix, gamma, J) = new{eltype(J),typeof(mass_matrix),
    typeof(gamma),typeof(J)}(mass_matrix,gamma,J,nothing,nothing)
end
function WOperator(f::DiffEqBase.AbstractSDEFunction, gamma)
  @assert DiffEqBase.has_jac(f) "f needs to have an associated jacobian"
  if isa(f, SplitSDEFunction)
    error("WOperator does not support SplitSDEFunction yet")
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
  return WOperator(mass_matrix, gamma, J)
end

set_gamma!(W::WOperator, gamma) = (W.gamma = gamma; W)
DiffEqBase.update_coefficients!(W::WOperator,u,p,t) = (update_coefficients!(W.J,u,p,t); W)
function Base.convert(::Type{AbstractMatrix}, W::WOperator)
  if W._concrete_form === nothing
    # Allocating
    W._concrete_form = W.mass_matrix - W.gamma * convert(AbstractMatrix,W.J)
  else
    # Non-allocating
    copyto!(W._concrete_form, W.mass_matrix)
    axpy!(-W.gamma, convert(AbstractMatrix,W.J), W._concrete_form)
  end
  W._concrete_form
end
Base.convert(::Type{Number}, W::WOperator) = W.mass_matrix - W.gamma * convert(Number,W.J)
Base.size(W::WOperator, args...) = size(W.J, args...)
Base.getindex(W::WOperator, i::Int) = W.mass_matrix[i] - W.gamma * W.J[i]
Base.getindex(W::WOperator, I::Vararg{Int,N}) where {N} =
  W.mass_matrix[I...] - W.gamma * W.J[I...]
Base.:*(W::WOperator, x::Union{AbstractVecOrMat,Number}) = W.mass_matrix*x - W.gamma * (W.J*x)
function Base.:\(W::WOperator, x::Union{AbstractVecOrMat,Number})
  if size(W) == () # scalar operator
    convert(Number,W) \ x
  else
    convert(AbstractMatrix,W) \ x
  end
end
function LinearAlgebra.mul!(Y::AbstractVecOrMat, W::WOperator, B::AbstractVecOrMat)
  if W._func_cache === nothing
    # Allocate cache only if needed
    W._func_cache = Vector{eltype(W)}(undef, size(Y, 1))
  end
  # Compute mass_matrix * B
  if isa(W.mass_matrix, UniformScaling)
    @.. Y = W.mass_matrix.λ * B
  else
    mul!(Y, W.mass_matrix, B)
  end
  # Compute J * B
  mul!(W._func_cache, W.J, B)
  # Subtract result
  axpy!(-W.gamma, W._func_cache, Y)
end

function calc_W!(integrator, cache::StochasticDiffEqMutableCache, γdt, repeat_step)
  @inbounds begin
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack J,W = cache
    is_compos = is_composite(integrator.alg)
    alg = unwrap_alg(integrator, true)
    mass_matrix = integrator.f.mass_matrix
    nlcache = cache.nlsolve.cache

    new_W = true
    if has_invW(f)
      # skip calculation of inv(W) if step is repeated
      !repeat_step && f.invW(W,uprev,p,γdt,t) # W == inverse W
      is_compos && calc_J!(integrator, cache, true)
    elseif has_jac(f) && f.jac_prototype !== nothing
      # skip calculation of J if step is repeated
      if repeat_step || (!integrator.last_stepfail && nlcache.nl_iters == 1 && nlcache.ηold < alg.new_jac_conv_bound)
        new_jac = false
      else # Compute a new Jacobian
        new_jac = true
        DiffEqBase.update_coefficients!(W,uprev,p,t)
      end
      # skip calculation of W if step is repeated
      if !repeat_step && (integrator.iter < 1 || new_jac || abs(dt - (t-integrator.tprev)) > 100eps(typeof(integrator.t)))
        set_gamma!(W, γdt)
      else
        new_W = false
      end
    else
      # skip calculation of J if step is repeated
      if repeat_step || (!integrator.last_stepfail && nlcache.nl_iters == 1 && nlcache.ηold < alg.new_jac_conv_bound)
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
    nlcache.new_W = new_W
    return nothing
  end
end

function calc_W!(integrator, cache::StochasticDiffEqConstantCache, γdt, repeat_step)
  @unpack t,uprev,p,f = integrator
  uf = cache.uf
  isarray = typeof(uprev) <: AbstractArray
  is_compos = is_composite(integrator.alg)
  mass_matrix = integrator.f.mass_matrix
  if has_jac(f)
    J = f.jac(uprev, p, t)
    if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
      J = DiffEqArrayOperator(J)
    end
    W = WOperator(mass_matrix, γdt, J)
  else
    J = calc_J(integrator, cache, is_compos)
    W = mass_matrix - γdt*J
  end
  is_compos && (integrator.eigen_est = isarray ? opnorm(J, Inf) : J)
  J, W
end
