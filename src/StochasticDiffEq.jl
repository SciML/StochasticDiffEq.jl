__precompile__()

module StochasticDiffEq

  import RandomNumbers: Xorshifts

  using Reexport
  @reexport using DiffEqBase

  using Parameters, RecursiveArrayTools, DataStructures
  using DiffEqNoiseProcess, Random
  using NLsolve, ForwardDiff, StaticArrays, MuladdMacro, DiffEqDiffTools

  import DiffEqBase: ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN,
         ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

  using DiffEqOperators: DiffEqArrayOperator

  import RecursiveArrayTools: chain

  using Logging

  using LinearAlgebra, Random

  import ForwardDiff.Dual

  import DiffEqBase: step!, initialize!, DEAlgorithm,
                     AbstractSDEAlgorithm, AbstractRODEAlgorithm, DEIntegrator, AbstractDiffEqInterpolation,
                     DECache, AbstractSDEIntegrator, AbstractRODEIntegrator, AbstractContinuousCallback,
                     Tableau

  # Integrator Interface
  import DiffEqBase: resize!,deleteat!,addat!,full_cache,user_cache,u_cache,du_cache,
                     rand_cache,ratenoise_cache,
                     resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
                     terminate!,get_du, get_dt,get_proposed_dt,set_proposed_dt!,
                     u_modified!,savevalues!,add_tstop!,add_saveat!,set_reltol!,
                     set_abstol!, postamble!, last_step_failed, has_invW, has_jac

  using DiffEqBase: check_error!, is_diagonal_noise, @..

  using DiffEqBase: nlsolvefail, isnewton, set_new_W!, get_W, @iipnlsolve, @oopnlsolve, _vec, _reshape

  using DiffEqBase: NLSolver

  using DiffEqBase: FastConvergence, Convergence, SlowConvergence, VerySlowConvergence, Divergence

  import DiffEqBase: calculate_residuals, calculate_residuals!, nlsolve_f, unwrap_cache, islinear


  const CompiledFloats = Union{Float32,Float64}

  include("misc_utils.jl")
  include("algorithms.jl")
  include("options_type.jl")
  include("derivative_wrappers.jl")
  include("interp_func.jl")
  include("caches/cache_types.jl")
  include("caches/basic_method_caches.jl")
  include("caches/lamba_caches.jl")
  include("caches/iif_caches.jl")
  include("caches/sdirk_caches.jl")
  include("caches/implicit_split_step_caches.jl")
  include("caches/sra_caches.jl")
  include("caches/rossler_caches.jl")
  include("caches/kencarp_caches.jl")
  include("caches/predcorr_caches.jl")
  include("integrators/type.jl")
  include("dense.jl")
  include("alg_utils.jl")
  include("integrators/integrator_utils.jl")
  include("cache_utils.jl")
  include("integrators/integrator_interface.jl")
  include("iterator_interface.jl")
  include("composite_solution.jl")
  include("solve.jl")
  include("initdt.jl")
  include("perform_step/low_order.jl")
  include("perform_step/lamba.jl")
  include("perform_step/iif.jl")
  include("perform_step/sri.jl")
  include("perform_step/sra.jl")
  include("perform_step/sdirk.jl")
  include("perform_step/implicit_split_step.jl")
  include("perform_step/kencarp.jl")
  include("perform_step/predcorr.jl")
  include("perform_step/split.jl")
  include("perform_step/composite.jl")
  include("tableaus.jl")
  include("derivative_utils.jl")
  include("composite_algs.jl")

   export StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
          StochasticCompositeAlgorithm

  export EM, LambaEM, PCEuler, RKMil, SRA, SRI, SRIW1,
         SRA1, SRA2, SRA3,
         SOSRA, SOSRA2, RKMilCommute,
         SRIW2, SOSRI, SOSRI2, SKenCarp,
         AutoSOSRI2, AutoSOSRA2

  export EulerHeun, LambaEulerHeun

  export SplitEM, IIF1M, IIF2M, IIF1Mil

  export ImplicitEM, ImplicitEulerHeun, ISSEM, ISSEulerHeun,
         ImplicitRKMil

  export StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
         StochasticDiffEqRODECompositeAlgorithm

  export RandomEM

  #General Functions
  export solve, init, solve!, step!

  #Misc Tools
  export checkSRIOrder, checkSRAOrder,  constructSRIW1, constructSRA1
end # module
