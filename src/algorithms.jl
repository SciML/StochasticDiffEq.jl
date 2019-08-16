abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

abstract type StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller} <: StochasticDiffEqAdaptiveAlgorithm end
abstract type StochasticDiffEqNewtonAlgorithm{CS,AD,Controller} <: StochasticDiffEqAlgorithm end

abstract type IteratedIntegralApprox end
struct IICommutative <:  IteratedIntegralApprox end
struct IIWiktorsson <:  IteratedIntegralApprox end

################################################################################

# Basics

struct EM{split} <: StochasticDiffEqAlgorithm end
EM(split=true) = EM{split}()

struct SplitEM <: StochasticDiffEqAlgorithm end
struct EulerHeun <: StochasticDiffEqAlgorithm end

struct LambaEM{split} <: StochasticDiffEqAdaptiveAlgorithm end
LambaEM(split=true) = LambaEM{split}()

struct LambaEulerHeun <: StochasticDiffEqAdaptiveAlgorithm end

struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(;interpretation=:Ito) = RKMil{interpretation}()

struct RKMilCommute{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMilCommute(;interpretation=:Ito) = RKMilCommute{interpretation}()

struct RKMil_General{T<:IteratedIntegralApprox} <: StochasticDiffEqAdaptiveAlgorithm
  interpretation::Symbol
  ii_approx::T
  c::Int
end
function RKMil_General(;interpretation=:Ito, ii_approx=IIWiktorsson(), c = 1)
  RKMil_General(interpretation, ii_approx, c)
end

struct WangLi3SMil_A <: StochasticDiffEqAlgorithm end
struct WangLi3SMil_B <: StochasticDiffEqAlgorithm end
struct WangLi3SMil_C <: StochasticDiffEqAlgorithm end
struct WangLi3SMil_D <: StochasticDiffEqAlgorithm end
struct WangLi3SMil_E <: StochasticDiffEqAlgorithm end
struct WangLi3SMil_F <: StochasticDiffEqAlgorithm end

#SROCK methods
struct SROCK1{interpretation} <: StochasticDiffEqAlgorithm end
SROCK1(;interpretation=:Ito) = SROCK1{interpretation}()

# Weak Order 2
struct SROCK2 <: StochasticDiffEqAlgorithm end
struct KomBurSROCK2 <: StochasticDiffEqAlgorithm end
struct SROCKC2 <: StochasticDiffEqAlgorithm end

# ROCK stabilization for EM
struct SROCKEM <: StochasticDiffEqAlgorithm
  strong_order_1::Bool
end
SROCKEM(;strong_order_1=true) = SROCKEM(strong_order_1)

struct SKSROCK <: StochasticDiffEqAlgorithm
  post_processing::Bool
end
SKSROCK(;post_processing=false) = SKSROCK(post_processing)

struct TangXiaoSROCK2 <: StochasticDiffEqAlgorithm
  version_num::Int
end
TangXiaoSROCK2(;version_num=5) = TangXiaoSROCK2(version_num)
###############################################################################

# Predictor Corrector
struct PCEuler{T<:Real, F} <: StochasticDiffEqAlgorithm
  theta::T
  eta::T
  ggprime::F
end


"""
    PCEuler(ggprime; theta=1/2, eta=1/2)

Predictor Corrector Euler

# Arguments
- `ggprime::Function`:
  For scalar problems, `ggprime` ``= b\\partial_x(b)``
  For multi-dimensional problems
  `bbprime_k` ``= \\sum_{j=1...M, i=1...D} b^(j)_i \\partial_i b^(j)_k``
  where ``b^(j)`` correspond to the noise vector due to the j'th noise channel.
  If problem is in place - a in place ggprime should be supplied - and
  vice versa for not in place speicification of problem.
- `theta::Real`:
  Degree of implicitness in the drift term. Set to 0.5 by default.
- `eta::Real`:
  Degree of implicitness in the diffusion term. Set to 0.5 by default.

Reference: Stochastics and Dynamics, Vol. 8, No. 3 (2008) 561–581
Note that the original paper has a typo in the definition of ggprime...
"""
PCEuler(ggprime; theta=1/2, eta=1/2) = PCEuler(theta,eta,ggprime)

################################################################################

# Rossler

struct SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType
end
SRA(;tableau=constructSRA1()) = SRA(tableau)

struct SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType
  error_terms::Int
end
SRI(;tableau=constructSRIW1(),error_terms=4) = SRI(tableau,error_terms)

struct SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end
struct SRIW2 <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRI <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRI2 <: StochasticDiffEqAdaptiveAlgorithm end
struct SRA1 <: StochasticDiffEqAdaptiveAlgorithm end
struct SRA2 <: StochasticDiffEqAdaptiveAlgorithm end
struct SRA3 <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRA <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRA2 <: StochasticDiffEqAdaptiveAlgorithm end

################################################################################

# IIF

struct IIF1M{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
IIF1M(;nlsolve=NLSOLVEJL_SETUP()) = IIF1M{typeof(nlsolve)}(nlsolve)

struct IIF2M{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
IIF2M(;nlsolve=NLSOLVEJL_SETUP()) = IIF2M{typeof(nlsolve)}(nlsolve)

struct IIF1Mil{F} <: StochasticDiffEqAlgorithm
  nlsolve::F
end
IIF1Mil(;nlsolve=NLSOLVEJL_SETUP()) = IIF1Mil{typeof(nlsolve)}(nlsolve)

################################################################################

# SDIRK

struct ImplicitEM{CS,AD,F,F2,S,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  nlsolve::F2
  diff_type::S
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ImplicitEM(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          theta = 1/2,symplectic=false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEM{chunk_size,autodiff,
                          typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                          typeof(new_jac_conv_bound),controller}(
                          linsolve,nlsolve,diff_type,
                          symplectic ? 1/2 : theta,
                          extrapolant,new_jac_conv_bound,symplectic)

struct ImplicitEulerHeun{CS,AD,F,S,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::S
  nlsolve::N
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ImplicitEulerHeun(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          theta = 1/2,symplectic = false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEulerHeun{chunk_size,autodiff,
                          typeof(linsolve),typeof(diff_type),
                          typeof(nlsolve),
                          typeof(new_jac_conv_bound),controller}(
                          linsolve,diff_type,nlsolve,
                          symplectic ? 1/2 : theta,
                          extrapolant,
                          new_jac_conv_bound,symplectic)

struct ImplicitRKMil{CS,AD,F,S,N,T2,Controller,interpretation} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::S
  nlsolve::N
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ImplicitRKMil(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          theta = 1/2,symplectic = false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive,interpretation=:Ito) =
                          ImplicitRKMil{chunk_size,autodiff,
                          typeof(linsolve),typeof(diff_type),
                          typeof(nlsolve),typeof(new_jac_conv_bound),
                          controller,interpretation}(
                          linsolve,diff_type,nlsolve,
                          symplectic ? 1/2 : theta,
                          extrapolant,
                          new_jac_conv_bound,symplectic)

struct ISSEM{CS,AD,F,S,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::S
  nlsolve::N
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ISSEM(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                       linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                       extrapolant=:constant,
                       theta = 1,symplectic=false,
                       new_jac_conv_bound = 1e-3,
                       controller = :Predictive) =
                       ISSEM{chunk_size,autodiff,
                       typeof(linsolve),typeof(diff_type),
                       typeof(nlsolve),
                       typeof(new_jac_conv_bound),controller}(
                       linsolve,diff_type,nlsolve,
                       symplectic ? 1/2 : theta,
                       extrapolant,
                       new_jac_conv_bound,symplectic)

struct ISSEulerHeun{CS,AD,F,S,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
 linsolve::F
 diff_type::S
 nlsolve::N
 theta::T2
 extrapolant::Symbol
 new_jac_conv_bound::T2
 symplectic::Bool
end
ISSEulerHeun(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                      extrapolant=:constant,
                      theta = 1,symplectic=false,
                      new_jac_conv_bound = 1e-3,
                      controller = :Predictive) =
                      ISSEulerHeun{chunk_size,autodiff,
                      typeof(linsolve),typeof(diff_type),
                      typeof(nlsolve),typeof(new_jac_conv_bound),controller}(
                      linsolve,diff_type,nlsolve,
                      symplectic ? 1/2 : theta,
                      extrapolant,
                      new_jac_conv_bound,symplectic)

struct SKenCarp{CS,AD,F,FDT,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller}
  linsolve::F
  diff_type::FDT
  nlsolve::N
  smooth_est::Bool
  extrapolant::Symbol
  new_jac_conv_bound::T2
  ode_error_est::Bool
end

SKenCarp(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   linsolve=DEFAULT_LINSOLVE,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:min_correct,
                   new_jac_conv_bound = 1e-3,controller = :Predictive,
                   ode_error_est = true) =
 SKenCarp{chunk_size,autodiff,typeof(linsolve),typeof(diff_type),
        typeof(nlsolve),typeof(new_jac_conv_bound),controller}(
        linsolve,diff_type,nlsolve,smooth_est,extrapolant,new_jac_conv_bound,
        ode_error_est)

################################################################################

# Etc.

struct StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

struct RandomEM <: StochasticDiffEqRODEAlgorithm end

const SplitSDEAlgorithms = Union{IIF1M,IIF2M,IIF1Mil,SKenCarp,SplitEM}
