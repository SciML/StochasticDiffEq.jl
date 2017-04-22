__precompile__()

module StochasticDiffEq

  import Base: linspace

  import Base: start, next, done, eltype

  using DiffEqBase, Parameters, RecursiveArrayTools, Juno,
        DataStructures, Iterators, Roots, DiffEqNoiseProcess

  using Compat

  import DiffEqBase: solve, solve!, init, step!, build_solution, initialize!

  import DiffEqBase: resize!,deleteat!,addat!,full_cache,user_cache, u_cache,du_cache,
                     resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
                     terminate!,get_du, get_dt,get_proposed_dt,
                     u_modified!,savevalues!,add_tstop!,add_saveat!,set_reltol!,
                     set_abstol!

  include("algorithms.jl")
  include("options_type.jl")
  include("interp_func.jl")
  include("constants.jl")
  include("caches.jl")
  include("integrators/type.jl")
  include("dense.jl")
  include("callbacks.jl")
  include("alg_utils.jl")
  include("integrators/integrator_utils.jl")
  include("integrators/integrator_interface.jl")
  include("iterator_interface.jl")
  include("solve.jl")
  include("initdt.jl")
  include("integrators/low_order.jl")
  include("integrators/sri.jl")
  include("integrators/sra.jl")
  include("integrators/split.jl")
  include("integrators/composite_integrator.jl")
  include("tableaus.jl")

   export StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
          StochasticCompositeAlgorithm

  export EM, RKMil, SRA, SRI, SRIW1, SRA1

  export EulerHeun

  export SplitEM

  export StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
         StochasticDiffEqRODECompositeAlgorithm

  export RandomEM

  #General Functions
  export solve, init, solve!, step!

  #Misc Tools
  export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
