abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

abstract type StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,Controller} <: StochasticDiffEqAdaptiveAlgorithm end
abstract type StochasticDiffEqNewtonAlgorithm{CS,AD,Controller} <: StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonAdaptiveAlgorithm{CS,AD,Controller} <: StochasticDiffEqJumpAdaptiveAlgorithm end

abstract type StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{CS,AD,Controller} <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm end

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

struct SimplifiedEM <: StochasticDiffEqAlgorithm end

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)
"""
struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(;interpretation=:Ito) = RKMil{interpretation}()

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)
"""
struct RKMilCommute{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMilCommute(;interpretation=:Ito) = RKMilCommute{interpretation}()

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)
"""
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
struct SROCK1{interpretation,E} <: StochasticDiffEqAlgorithm
  eigen_est::E
end
SROCK1(;interpretation=:Ito,eigen_est=nothing) = SROCK1{interpretation,typeof(eigen_est)}(eigen_est)

# Weak Order 2
for Alg in [:SROCK2, :KomBurSROCK2, :SROCKC2]
  @eval begin
    struct $Alg{E} <: StochasticDiffEqAlgorithm
      eigen_est::E
    end
    $Alg(;eigen_est=nothing) = $Alg(eigen_est)
  end
end

# ROCK stabilization for EM
struct SROCKEM{E} <: StochasticDiffEqAlgorithm
  strong_order_1::Bool
  eigen_est::E
end
SROCKEM(;strong_order_1=true,eigen_est=nothing) = SROCKEM(strong_order_1,eigen_est)

struct SKSROCK{E} <: StochasticDiffEqAlgorithm
  post_processing::Bool
  eigen_est::E
end
SKSROCK(;post_processing=false,eigen_est=nothing) = SKSROCK(post_processing,eigen_est)

struct TangXiaoSROCK2{E} <: StochasticDiffEqAlgorithm
  version_num::Int
  eigen_est::E
end
TangXiaoSROCK2(;version_num=5,eigen_est=nothing) = TangXiaoSROCK2(version_num,eigen_est)
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

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType
end
SRA(;tableau=constructSRA1()) = SRA(tableau)

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRI{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType
  error_terms::Int
end
SRI(;tableau=constructSRIW1(),error_terms=4) = SRI(tableau,error_terms)

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRIW2 <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRI <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRI2 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRA2 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X
"""
struct SRA3 <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRA <: StochasticDiffEqAdaptiveAlgorithm end
struct SOSRA2 <: StochasticDiffEqAdaptiveAlgorithm end

################################################################################

# Rossler second order for weak approx.

struct DRI1 <: StochasticDiffEqAdaptiveAlgorithm end
struct RI1 <: StochasticDiffEqAdaptiveAlgorithm end
struct RI3 <: StochasticDiffEqAdaptiveAlgorithm end
struct RI5 <: StochasticDiffEqAdaptiveAlgorithm end
struct RI6 <: StochasticDiffEqAdaptiveAlgorithm end
struct RDI1WM <: StochasticDiffEqAlgorithm end
struct RDI2WM <: StochasticDiffEqAdaptiveAlgorithm end
struct RDI3WM <: StochasticDiffEqAdaptiveAlgorithm end
struct RDI4WM <: StochasticDiffEqAdaptiveAlgorithm end

# Stratonovich sense

struct RS1 <: StochasticDiffEqAlgorithm end
struct RS2 <: StochasticDiffEqAlgorithm end

struct PL1WM <: StochasticDiffEqAlgorithm end
struct PL1WMA <: StochasticDiffEqAlgorithm end

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
                          theta = 1,symplectic=false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEM{chunk_size,autodiff,
                          typeof(linsolve),typeof(nlsolve),typeof(diff_type),
                          typeof(new_jac_conv_bound),controller}(
                          linsolve,nlsolve,diff_type,
                          symplectic ? 1/2 : theta,
                          extrapolant,new_jac_conv_bound,symplectic)

STrapezoid(;kwargs...) = ImplicitEM(;theta=1/2,kwargs...)
SImplicitMidpoint(;kwargs...) = ImplicitEM(;theta=1/2,symplectic=true,kwargs...)

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
                          theta = 1,symplectic = false,
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
                          theta = 1,symplectic = false,
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

# Jumps

struct TauLeaping <: StochasticDiffEqJumpAdaptiveAlgorithm end
struct CaoTauLeaping <: StochasticDiffEqJumpAdaptiveAlgorithm end

################################################################################

# Etc.

struct StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

struct RandomEM <: StochasticDiffEqRODEAlgorithm end

const SplitSDEAlgorithms = Union{IIF1M,IIF2M,IIF1Mil,SKenCarp,SplitEM}
