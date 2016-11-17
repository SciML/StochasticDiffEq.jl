"""
`monteCarloSim(dt::Number,prob::SDEProblem,alg::AbstractSDEAlgorithm)`

Performs a parallel Monte-Carlo simulation to solve the SDE problem with dt numMonte times.
Returns a vector of solution objects.

### Keyword Arguments
* `T` - Final time. Default is 1.
* `numMonte` - Number of Monte-Carlo simulations to run. Default is 10000
* `save_timeseries` - Denotes whether save_timeseries should be turned on in each run. Default is false.
"""
function monteCarloSim(prob::AbstractSDEProblem,alg;numMonte=10000,save_timeseries=false,kwargs...)
  elapsedTime = @elapsed solutions = pmap((i)->solve(prob,alg;save_timeseries=save_timeseries,kwargs...),1:numMonte)
  if typeof(prob) <: SDEProblem
    solutions = convert(Array{SDESolution},solutions)
  elseif typeof(prob) <: SDETestProblem
    solutions = convert(Array{SDETestSolution},solutions)
  end
  if typeof(prob) <: SDETestProblem
    N = size(solutions,1)
    errors = Dict() #Should add type information
    error_means  = Dict()
    error_medians= Dict()
    for k in keys(solutions[1].errors)
      errors[k] = reshape(Float64[sol.errors[k] for sol in solutions],size(solutions)...)
      error_means[k] = mean(errors[k])
      error_medians[k]=median(errors[k])
    end
  end
  return(MonteCarloSimulation(solutions,errors,error_means,error_medians,N,elapsedTime))
end

type MonteCarloSimulation
  solutions#::Array{T}
  errors
  error_means
  error_medians
  N
  elapsedTime
end

Base.length(sim::MonteCarloSimulation) = sim.N
Base.endof( sim::MonteCarloSimulation) = length(sim)
Base.getindex(sim::MonteCarloSimulation,i::Int) = sim.solutions[i]
Base.getindex(sim::MonteCarloSimulation,i::Int,I::Int...) = sim.solutions[i][I]

function print(io::IO, sim::MonteCarloSimulation)
  println(io,"$(typeof(sim)) of length $(length(sim)).")
  println(io,"\n-----------Errors-----------")
  for (k,v) in sim.errors
    println(io,"$k: $v")
  end
end

function show(io::IO,sim::MonteCarloSimulation)
  println(io,"$(typeof(sim)) of length $(length(sim)).")
end
