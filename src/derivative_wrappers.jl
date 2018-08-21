function derivative(f, x::Union{Number,AbstractArray{<:Number}},
                    integrator)
    local d
    alg = unwrap_alg(integrator, true)
    if get_current_alg_autodiff(integrator.alg, integrator.cache)
      d = ForwardDiff.derivative(f, x)
    else
      d = DiffEqDiffTools.finite_difference_gradient(f, x, alg.diff_type, eltype(x), Val{false})
    end
    d
end

function derivative!(df::AbstractArray{<:Number}, f, x::Union{Number,AbstractArray{<:Number}}, fx::AbstractArray{<:Number}, integrator::DEIntegrator)
    if alg_autodiff(integrator.alg)
        ForwardDiff.derivative!(df, f, fx, x)
    else
        RealOrComplex = eltype(integrator.u) <: Complex ? Val{:Complex} : Val{:Real}
        DiffEqDiffTools.finite_difference!(df, f, x, integrator.alg.diff_type, RealOrComplex, fx)
    end
    nothing
end

function jacobian(f, x,
                  integrator::DiffEqBase.DEIntegrator)
    local J
    alg = unwrap_alg(integrator, true)
    isarray = typeof(x) <: AbstractArray
    if get_current_alg_autodiff(integrator.alg, integrator.cache)
      if isarray
        J = ForwardDiff.jacobian(f,x)
      else
        J = ForwardDiff.derivative(f,x)
      end
    else
      if isarray
        J = DiffEqDiffTools.finite_difference_jacobian(f, x, alg.diff_type, eltype(x), Val{false})
      else
        J = DiffEqDiffTools.finite_difference_derivative(f, x, alg.diff_type, eltype(x))
      end
    end
    J
end

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DEIntegrator, jac_config)
    if alg_autodiff(integrator.alg)
      ForwardDiff.jacobian!(J, f, fx, x, jac_config)
    else
      DiffEqDiffTools.finite_difference_jacobian!(J, f, x, jac_config)
    end
    nothing
end

function build_jac_config(alg,f,uf,du1,uprev,u,tmp,du2)
  if !has_jac(f)
    if alg_autodiff(alg)
      jac_config = ForwardDiff.JacobianConfig(uf,du1,uprev,ForwardDiff.Chunk{determine_chunksize(u,alg)}())
    else
      if alg.diff_type != Val{:complex}
        jac_config = DiffEqDiffTools.JacobianCache(tmp,du1,du2,alg.diff_type)
      else
        jac_config = DiffEqDiffTools.JacobianCache(Complex{eltype(tmp)}.(tmp),Complex{eltype(du1)}.(du1),nothing,alg.diff_type)
      end
    end
  else
    jac_config = nothing
  end
  jac_config
end
