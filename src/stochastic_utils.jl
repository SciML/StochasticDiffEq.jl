"""
`monte_carlo_simulation(dt::Number,prob::AbstractSDEProblem,alg::AbstractSDEAlgorithm)`

Performs a parallel Monte-Carlo simulation to solve the SDE problem with dt numMonte times.
Returns a vector of solution objects.

### Keyword Arguments
* `T` - Final time. Default is 1.
* `numMonte` - Number of Monte-Carlo simulations to run. Default is 10000
* `save_timeseries` - Denotes whether save_timeseries should be turned on in each run. Default is false.
"""
function monte_carlo_simulation(prob::SDETestProblem,alg;numMonte=10000,save_timeseries=false,kwargs...)
  elapsedTime = @elapsed solutions = pmap((i)->solve(prob,alg;save_timeseries=save_timeseries,kwargs...),1:numMonte)
  solutions = convert(Array{SDETestSolution},solutions)
  N = size(solutions,1)
  errors = Dict() #Should add type information
  error_means  = Dict()
  error_medians= Dict()
  for k in keys(solutions[1].errors)
    errors[k] = reshape(Float64[sol.errors[k] for sol in solutions],size(solutions)...)
    error_means[k] = mean(errors[k])
    error_medians[k]=median(errors[k])
  end
  return(MonteCarloTestSimulation(solutions,errors,error_means,error_medians,elapsedTime))
end

function monte_carlo_simulation(prob::SDEProblem,alg;numMonte=10000,save_timeseries=false,kwargs...)
  elapsedTime = @elapsed solutions = pmap((i)->solve(prob,alg;save_timeseries=save_timeseries,kwargs...),1:numMonte)
  solutions = convert(Array{SDESolution},solutions)
  return(MonteCarloSimulation(solutions,elapsedTime))
end

type MonteCarloTestSimulation
  solutions#::Array{T}
  errors
  error_means
  error_medians
  elapsedTime
end

type MonteCarloSimulation
  solutions#::Array{T}
  elapsedTime
end


Base.length(sim::AbstractMonteCarloSimulation) = length(sim.solutions)
Base.endof( sim::AbstractMonteCarloSimulation) = length(sim)
Base.getindex(sim::AbstractMonteCarloSimulation,i::Int) = sim.solutions[i]
Base.getindex(sim::AbstractMonteCarloSimulation,i::Int,I::Int...) = sim.solutions[i][I]
