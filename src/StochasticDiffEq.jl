module StochasticDiffEq

  using DiffEqBase, Parameters, ChunkedArrays, RecursiveArrayTools, Juno
  import DiffEqBase: solve

  macro def(name, definition)
      quote
          macro $name()
              esc($(Expr(:quote, definition)))
          end
      end
  end

  include("algorithms.jl")
  include("sde/sde_noise_process.jl")

  include("stochastic_utils.jl")
  include("sde/sde_solve.jl")
  include("sde/sde_integrators.jl")
  include("sde/sde_tableaus.jl")

   export StochasticDiffEqAlgorithm, EM, RKMil, SRA, SRI, SRIW1Optimized,
          SRA1Optimized, SRAVectorized, SRIVectorized

   #Stochastic Utils
   export monteCarloSim, construct_correlated_noisefunc, WHITE_NOISE, NoiseProcess

   #General Functions
   export solve

   #Misc Tools
   export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
