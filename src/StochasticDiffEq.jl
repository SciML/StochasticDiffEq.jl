module StochasticDiffEq

  using DiffEqBase, ParameterizedFunctions, Parameters, ChunkedArrays
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
  include("premade_problems.jl")

  include("stochastic_utils.jl")
  include("sde/sde_solve.jl")
  include("sde/sde_integrators.jl")
  include("sde/sde_tableaus.jl")

  export SDEProblem, SDESolution

  #SDE Example Problems
  export prob_sde_wave, prob_sde_linear, prob_sde_cubic, prob_sde_2Dlinear, prob_sde_lorenz,
         prob_sde_2Dlinear, prob_sde_additive, prob_sde_additivesystem, oval2ModelExample

   #Stochastic Utils
   export monteCarloSim, construct_correlated_noisefunc, WHITE_NOISE, NoiseProcess

   #General Functions
   export solve

   #Misc Tools
   export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
