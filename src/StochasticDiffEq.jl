__precompile__()

module StochasticDiffEq

  using DiffEqBase, Parameters, ChunkedArrays, RecursiveArrayTools, Juno,
        DataStructures, ResettableStacks
  import DiffEqBase: solve


  abstract AbstractMonteCarloSimulation

  macro def(name, definition)
      quote
          macro $name()
              esc($(Expr(:quote, definition)))
          end
      end
  end

  include("algorithms.jl")
  include("stochastic_utils.jl")
  include("solve.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/low_order.jl")
  include("integrators/sri.jl")
  include("integrators/sra.jl")
  include("tableaus.jl")

   export StochasticDiffEqAlgorithm, EM, RKMil, SRA, SRI, SRIW1,
          SRA1

   #Stochastic Utils
   export monte_carlo_simulation

   #General Functions
   export solve

   #Misc Tools
   export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
