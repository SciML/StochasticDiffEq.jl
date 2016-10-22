"""
`SDESolution`

Holds the data for the solution to a SDE problem.

### Fields

* `u::Array{Float64}`: The solution (at the final timepoint)
* `trueknown::Bool`: Boolean flag for if the true solution is given.
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `timeseries`::AbstractArrayOrVoid`: u over time. Only saved if `save_timeseries=true`
  is specified in the solver.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `Ws`: All of the W's in the solution. Only saved if `save_timeseries=true` is specified
  in the solver.
* `timeseries_analytic`: If `save_timeseries=true`, saves the solution at each save point.
* `prob::DEProblem`: Holds the problem object used to define the problem.
* `save_timeseries::Bool`: True if solver saved the extra timepoints.
* `appxtrue::Bool`: Boolean flag for if u_analytic was an approximation.

"""
type SDESolution <: AbstractODESolution
  u#::AbstractArrayOrNumber
  trueknown::Bool
  u_analytic#::AbstractArrayOrNumber
  errors#::Dict{}
  timeseries#::AbstractArrayOrVoid
  t#::AbstractArrayOrVoid
  Δt#::AbstractArrayOrVoid
  Ws#::AbstractArrayOrVoid
  timeseries_analytic#::AbstractArrayOrVoid
  appxtrue::Bool
  save_timeseries::Bool
  maxstacksize::Int
  W
  function SDESolution(u::Union{AbstractArray,Number};timeseries=[],timeseries_analytic=[],t=[],Δt=[],Ws=[],maxstacksize=0,W=0.0)
    save_timeseries = timeseries == nothing
    trueknown = false
    return(new(u,trueknown,nothing,Dict(),timeseries,t,Δt,Ws,timeseries_analytic,false,save_timeseries,maxstacksize,W))
  end
  function SDESolution(u,u_analytic;timeseries=[],timeseries_analytic=[],t=[],Δt=nothing,Ws=[],maxstacksize=0,W=0.0)
    save_timeseries = timeseries != []
    trueknown = true
    errors = Dict(:final=>mean(abs.(u-u_analytic)))
    if save_timeseries
      errors = Dict(:final=>mean(abs.(u-u_analytic)),:l∞=>maximum(vecvecapply((x)->abs.(x),timeseries-timeseries_analytic)),:l2=>sqrt(mean(vecvecapply((x)->x.^2,timeseries-timeseries_analytic))))
    end
    return(new(u,trueknown,u_analytic,errors,timeseries,t,Δt,Ws,timeseries_analytic,false,save_timeseries,maxstacksize,W))
  end
  #Required to convert pmap results
  SDESolution(a::Any) = new(a.u,a.trueknown,a.u_analytic,a.errors,a.timeseries,a.t,a.Δt,a.Ws,a.timeseries_analytic,a.appxtrue,a.save_timeseries,a.maxstacksize,a.W)
end
