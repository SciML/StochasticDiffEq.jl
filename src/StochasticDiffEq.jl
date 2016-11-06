module StochasticDiffEq

  using DiffEqBase, Parameters, ChunkedArrays, RecursiveArrayTools
  import DiffEqBase: solve

  macro def(name, definition)
      quote
          macro $name()
              esc($(Expr(:quote, definition)))
          end
      end
  end


  include("sde/sde_noise_process.jl")
  include("problems.jl")
  include("solutions.jl")

  include("stochastic_utils.jl")
  include("sde/sde_solve.jl")
  include("sde/sde_integrators.jl")
  include("sde/sde_tableaus.jl")

  export SDEProblem, SDESolution

   #Stochastic Utils
   export monteCarloSim, construct_correlated_noisefunc, WHITE_NOISE, NoiseProcess

   #General Functions
   export solve

   #Misc Tools
   export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
