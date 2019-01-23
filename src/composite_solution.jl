struct RODECompositeSolution{T,N,uType,uType2,EType,tType,randType,P,A,IType} <: DiffEqBase.AbstractRODESolution{T,N}
  u::uType
  u_analytic::uType2
  errors::EType
  t::tType
  W::randType
  prob::P
  alg::A
  interp::IType
  alg_choice::Vector{Int}
  dense::Bool
  tslocation::Int
  retcode::Symbol
  seed::UInt64
end
(sol::RODECompositeSolution)(t,deriv::Type=Val{0};idxs=nothing,continuity=:left) = sol.interp(t,idxs,deriv,sol.prob.p,continuity)
(sol::RODECompositeSolution)(v,t,deriv::Type=Val{0};idxs=nothing,continuity=:left) = sol.interp(v,t,idxs,deriv,sol.prob.p,continuity)

function DiffEqBase.build_solution(
        prob::DiffEqBase.AbstractRODEProblem,
        alg::Union{StochasticDiffEqCompositeAlgorithm,
                   StochasticDiffEqRODECompositeAlgorithm},t,u;
        W=[],timeseries_errors=length(u)>2,
        dense=false,dense_errors=dense,
        calculate_error = true,
        interp = LinearInterpolation(t,u),
        alg_choice=[1], retcode = :Default,seed = UInt64(0),
        kwargs...)

  T = eltype(eltype(u))
  if typeof(prob.u0) <: Tuple
    N = length((size(ArrayPartition(prob.u0))..., length(u)))
  else
    N = length((size(prob.u0)..., length(u)))
  end

  if typeof(prob.f) <: Tuple
    f = prob.f[1]
  else
    f = prob.f
  end

  if DiffEqBase.has_analytic(f)
    u_analytic = Vector{typeof(prob.u0)}()
    errors = Dict{Symbol,real(eltype(prob.u0))}()
    sol = RODECompositeSolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(W),
                       typeof(prob),typeof(alg),typeof(interp)}(u,u_analytic,errors,t,W,prob,alg,interp,alg_choice,dense,0,retcode,seed)
    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return RODECompositeSolution{T,N,typeof(u),typeof(nothing),typeof(nothing),typeof(t),typeof(W),
                       typeof(prob),typeof(alg),typeof(interp)}(u,nothing,nothing,t,W,prob,alg,interp,alg_choice,dense,0,retcode,seed)
  end
end

function DiffEqBase.solution_new_retcode(sol::RODECompositeSolution{T,N,uType,uType2,EType,tType,randType,P,A,IType},retcode) where {T,N,uType,uType2,EType,tType,randType,P,A,IType}
  RODECompositeSolution{T,N,uType,uType2,EType,tType,randType,P,A,IType}(
                        sol.u,sol.u_analytic,sol.errors,sol.t,sol.W,sol.prob,
                        sol.alg,sol.interp,sol.alg_choice,sol.dense,sol.tslocation,
                        retcode,sol.seed)
end

function DiffEqBase.solution_new_tslocation(sol::RODECompositeSolution{T,N,uType,uType2,EType,tType,randType,P,A,IType},tslocation) where {T,N,uType,uType2,EType,tType,randType,P,A,IType}
  RODECompositeSolution{T,N,uType,uType2,EType,tType,randType,P,A,IType}(
                        sol.u,sol.u_analytic,sol.errors,sol.t,sol.W,sol.prob,
                        sol.alg,sol.interp,sol.alg_choice,sol.dense,tslocation,
                        sol.retcode,sol.seed)
end
