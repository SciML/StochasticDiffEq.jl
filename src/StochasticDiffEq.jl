__precompile__()

module StochasticDiffEq

  import Base: linspace

  import Base: start, next, done, eltype

  import RandomNumbers: Xorshifts

  using Reexport
  @reexport using DiffEqBase

  using Parameters, RecursiveArrayTools, Juno,
        DataStructures, Roots, DiffEqNoiseProcess,
        NLsolve, ForwardDiff, StaticArrays, MuladdMacro

  import DiffEqBase: ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

  import RecursiveArrayTools: chain

  using Compat

  import ForwardDiff.Dual

  import DiffEqBase: solve, solve!, init, step!, build_solution, initialize!

  import DiffEqBase: resize!,deleteat!,addat!,full_cache,user_cache, u_cache,du_cache,
                     resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
                     terminate!,get_du, get_dt,get_proposed_dt,
                     u_modified!,savevalues!,add_tstop!,add_saveat!,set_reltol!,
                     set_abstol!

  macro tight_loop_macros(ex)
   :($(esc(ex)))
  end

  include("misc_utils.jl")
  include("algorithms.jl")
  include("options_type.jl")
  include("interp_func.jl")
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
  include("integrators/iif.jl")
  include("integrators/sri.jl")
  include("integrators/sra.jl")
  include("integrators/split.jl")
  include("integrators/composite_integrator.jl")
  include("tableaus.jl")

   export StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
          StochasticCompositeAlgorithm

  export EM, RKMil, SRA, SRI, SRIW1, SRA1

  export EulerHeun

  export SplitEM, IIF1M, IIF2M, IIF1Mil

  export StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
         StochasticDiffEqRODECompositeAlgorithm

  export RandomEM

  #General Functions
  export solve, init, solve!, step!

  #Misc Tools
  export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1

end # module
