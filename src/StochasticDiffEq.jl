__precompile__()

module StochasticDiffEq

  using DiffEqBase, Parameters, ChunkedArrays, RecursiveArrayTools, Juno,
        DataStructures, ResettableStacks
  import DiffEqBase: solve

  include("algorithms.jl")
  include("solve.jl")
  include("initdt.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/low_order.jl")
  include("integrators/sri.jl")
  include("integrators/sra.jl")
  include("tableaus.jl")

   export StochasticDiffEqAlgorithm, EM, RKMil, SRA, SRI, SRIW1,
          SRA1

   #General Functions
   export solve

   #Misc Tools
   export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
