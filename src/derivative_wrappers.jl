function derivative(f, x::Union{Number,AbstractArray{<:Number}},
                    integrator)
    local d
    alg = unwrap_alg(integrator, true)
    if get_current_alg_autodiff(integrator.alg, integrator.cache)
      d = ForwardDiff.derivative(f, x)
    else
      d = DiffEqDiffTools.finite_difference_gradient(f, x, alg.diff_type, eltype(x), Val(false))
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

jacobian_autodiff(f,x,_,_)=ForwardDiff.derivative(f,x)
jacobian_autodiff(f,x::AbstractArray,_,hascolorvec::Val{false})=ForwardDiff.jacobian(f,x)
function jacobian_autodiff(f,x::AbstractArray,integrator,hascolorvec::Val{true})
  colorvec=integrator.f.colorvec
  jac=integrator.f.jac_prototype
  J=jac isa SparseMatrixCSC ? similar(jac) : fill(0.,size(jac))
  forwarddiff_color_jacobian!(J,f,x,colorvec=colorvec,sparsity=jac)
  J
end
jacobian_finitediff(f,x,difftype,_,_)=DiffEqDiffTools.finite_difference_derivative(f, x, difftype, eltype(x))
jacobian_finitediff(f,x::AbstractArray,difftype,_,hascolorvec::Val{false})=DiffEqDiffTools.finite_difference_jacobian(f, x, difftype, eltype(x), Val{false})
jacobian_finitediff(f,x::AbstractArray,difftype,integrator,hascolorvec::Val{true})=
  DiffEqDiffTools.finite_difference_jacobian(f, x, difftype, eltype(x), Val(false),colorvec=integrator.f.colorvec,sparsity=integrator.f.jac_prototype)

function jacobian(f, x,
                  integrator::DiffEqBase.DEIntegrator)
    local J
    alg = unwrap_alg(integrator, true)
    if get_current_alg_autodiff(integrator.alg, integrator.cache)
      J = jacobian_autodiff(f,x,integrator,Val(DiffEqBase.has_colorvec(integrator.f)))
    else
      J = jacobian_finitediff(f,x,alg.diff_type,integrator,Val(DiffEqBase.has_colorvec(integrator.f)))
    end
    J
end

jacobian_autodiff!(J,f,fx,x,jac_config::ForwardColorJacCache)=forwarddiff_color_jacobian!(J,f,x,jac_config)
jacobian_autodiff!(J,f,fx,x,jac_config::ForwardDiff.JacobianConfig)=ForwardDiff.jacobian!(J, f, fx, x, jac_config)

function jacobian!(J::AbstractMatrix{<:Number}, f, x::AbstractArray{<:Number}, fx::AbstractArray{<:Number}, integrator::DEIntegrator, jac_config)
    if alg_autodiff(integrator.alg)
      jacobian_autodiff!(J, f, fx, x, jac_config)
    else
      DiffEqDiffTools.finite_difference_jacobian!(J, f, x, jac_config)
    end
    nothing
end

jac_cache_autodiff(alg,f,uf,du1,uprev,u,hascolorvec::Val{false})=ForwardDiff.JacobianConfig(uf,du1,uprev,ForwardDiff.Chunk{determine_chunksize(u,alg)}())
jac_cache_autodiff(alg,f,uf,du1,uprev,u,hascolorvec::Val{true})=ForwardColorJacCache(uf,uprev,colorvec=f.colorvec,sparsity=f.jac_prototype)

function DiffEqBase.build_jac_config(alg::StochasticDiffEqAlgorithm,f,uf,du1,uprev,u,tmp,du2)
  if !has_jac(f)
    if alg_autodiff(alg)
      jac_config = jac_cache_autodiff(alg,f,uf,du1,uprev,u,Val(DiffEqBase.has_colorvec(f)))
    else
      colorvec= f.colorvec isa Nothing ? Base.OneTo(length(u)) : f.colorvec
      if alg.diff_type != Val{:complex}
        jac_config = DiffEqDiffTools.JacobianCache(tmp,du1,du2,alg.diff_type,colorvec=colorvec,sparsity=f.jac_prototype)
      else
        jac_config = DiffEqDiffTools.JacobianCache(Complex{eltype(tmp)}.(tmp),Complex{eltype(du1)}.(du1),nothing,alg.diff_type,colorvec=colorvec,sparsity=f.jac_prototype)
      end
    end
  else
    jac_config = nothing
  end
  jac_config
end
