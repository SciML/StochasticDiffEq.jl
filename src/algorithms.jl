abstract type StochasticDiffEqAlgorithm <: AbstractSDEAlgorithm end
abstract type StochasticDiffEqAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqCompositeAlgorithm <: StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqRODEAlgorithm <: AbstractRODEAlgorithm end
abstract type StochasticDiffEqRODEAdaptiveAlgorithm <: StochasticDiffEqRODEAlgorithm end
abstract type StochasticDiffEqRODECompositeAlgorithm <: StochasticDiffEqRODEAlgorithm end

abstract type StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller} <: StochasticDiffEqAdaptiveAlgorithm end
abstract type StochasticDiffEqNewtonAlgorithm{CS,AD,FDT,ST,CJ,Controller} <: StochasticDiffEqAlgorithm end

abstract type StochasticDiffEqJumpAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller} <: StochasticDiffEqJumpAdaptiveAlgorithm end

abstract type StochasticDiffEqJumpDiffusionAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpDiffusionAdaptiveAlgorithm <: StochasticDiffEqAlgorithm end
abstract type StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller} <: StochasticDiffEqJumpDiffusionAdaptiveAlgorithm end

abstract type IteratedIntegralApprox end
struct IICommutative <:  IteratedIntegralApprox end
struct IILevyArea <:  IteratedIntegralApprox end

################################################################################

# Basics
"""
EM: Nonstiff Method
The Euler-Maruyama method. Strong Order 0.5 in the Ito sense.
Has an optional argument split=true for controlling step splitting.
When splitting is enabled, the stability with large diffusion eigenvalues is improved.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
Fixed time step only.
"""
struct EM{split} <: StochasticDiffEqAlgorithm end
EM(split=true) = EM{split}()

struct SplitEM <: StochasticDiffEqAlgorithm end
"""
EulerHeun: Nonstiff Method
The Euler-Heun method.
Strong Order 0.5 in the Stratonovich sense.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
Fixed time step only.
"""
struct EulerHeun <: StochasticDiffEqAlgorithm end
"""
LambaEM: Nonstiff Method
A modified Euler-Maruyama method with adaptive time stepping with an error estimator based on Lamba and Rackauckas.
Has an optional argument split=true for controlling step splitting.
When splitting is enabled, the stability with large diffusion eigenvalues is improved.
Strong Order 0.5 in the Ito sense. Can handle all forms of noise, including non-diagonal, scalar, and colored noise
"""
struct LambaEM{split} <: StochasticDiffEqAdaptiveAlgorithm end
LambaEM(split=true) = LambaEM{split}()
"""
LambaEulerHeun: Nonstiff Method
A modified Euler-Heun method with adaptive time stepping with an error estimator based on Lamba due to Rackauckas.
Strong order 0.5 in the Stratonovich sense.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
"""
struct LambaEulerHeun <: StochasticDiffEqAdaptiveAlgorithm end

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)

SimplifiedEM: High Weak Order Method
A simplified Euler-Maruyama method with weak order 1.0 and fixed step size.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise.
"""
struct SimplifiedEM <: StochasticDiffEqAlgorithm end

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)

RKMil: Nonstiff Method
An explicit Runge-Kutta discretization of the strong order 1.0 Milstein method.
Defaults to solving the Ito problem, but RKMil(interpretation=SciMLBase.AlgorithmInterpretation.Stratonovich) makes it solve the Stratonovich problem.
Only handles scalar and diagonal noise.
"""
struct RKMil{interpretation} <: StochasticDiffEqAdaptiveAlgorithm end
RKMil(;interpretation=SciMLBase.AlgorithmInterpretation.Ito) = RKMil{interpretation}()

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)

RKMilCommute: Nonstiff Method
An explicit Runge-Kutta discretization of the strong order 1.0 Milstein method for commutative noise problems.
Defaults to solving the Ito problem, but RKMilCommute(interpretation=SciMLBase.AlgorithmInterpretation.Stratonovich) makes it solve the Stratonovich problem.
Uses a 1.5/2.0 error estimate for adaptive time stepping.
Default: ii_approx=IICommutative() does not approximate the Levy area.
"""
struct RKMilCommute{T} <: StochasticDiffEqAdaptiveAlgorithm
  interpretation::SciMLBase.AlgorithmInterpretation.T
  ii_approx::T
end
RKMilCommute(;interpretation=SciMLBase.AlgorithmInterpretation.Ito, ii_approx=IICommutative()) = RKMilCommute(interpretation,ii_approx)

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)

RKMilGeneral: Nonstiff Method
RKMilGeneral(;interpretation=SciMLBase.AlgorithmInterpretation.Ito, ii_approx=IILevyArea()
An explicit Runge-Kutta discretization of the strong order 1.0 Milstein method for general non-commutative noise problems.
Allows for a choice of interpretation between SciMLBase.AlgorithmInterpretation.Ito and SciMLBase.AlgorithmInterpretation.Stratonovich.
Allows for a choice of iterated integral approximation.
Default: ii_approx=IILevyArea() uses LevyArea.jl to choose optimal algorithm. See
Kastner, F. and Rößler, A., arXiv: 2201.08424
Kastner, F. and Rößler, A., LevyArea.jl, 10.5281/ZENODO.5883748, https://github.com/stochastics-uni-luebeck/LevyArea.jl
"""
struct RKMilGeneral{T, TruncationType} <: StochasticDiffEqAdaptiveAlgorithm
  interpretation::SciMLBase.AlgorithmInterpretation.T
  ii_approx::T
  c::Int
  p::TruncationType
end

function RKMilGeneral(;interpretation=SciMLBase.AlgorithmInterpretation.Ito,ii_approx=IILevyArea(), c=1, p=nothing, dt=nothing)
  γ = 1//1
  p==true && (p = Int(floor(c*dt^(1//1-2//1*γ)) + 1))
  RKMilGeneral{typeof(ii_approx), typeof(p)}(interpretation, ii_approx, c, p)
end

"""
WangLi3SMil_A: Nonstiff Method
Fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
"""
struct WangLi3SMil_A <: StochasticDiffEqAlgorithm end
"""
WangLi3SMil_B: Nonstiff Method
Fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
"""
struct WangLi3SMil_B <: StochasticDiffEqAlgorithm end
"""
WangLi3SMil_C: Nonstiff Method
Fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
"""
struct WangLi3SMil_C <: StochasticDiffEqAlgorithm end
"""
WangLi3SMil_D: Nonstiff Method
Fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
"""
struct WangLi3SMil_D <: StochasticDiffEqAlgorithm end
"""
WangLi3SMil_E: Nonstiff Method
Fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
"""
struct WangLi3SMil_E <: StochasticDiffEqAlgorithm end
"""
WangLi3SMil_F: Nonstiff Method
Fixed step-size explicit 3-stage Milstein methods for Ito problem with strong and weak order 1.0
"""
struct WangLi3SMil_F <: StochasticDiffEqAlgorithm end

#SROCK methods
"""
SROCK1: S-ROCK Method
Is a fixed step size stabilized explicit method for stiff problems.
Defaults to solving th Ito problem but SROCK1(interpretation=SciMLBase.AlgorithmInterpretation.Stratonovich) can make it solve the Stratonovich problem.
Strong order of convergence is 0.5 and weak order 1, but is optimised to get order 1 in case os scalar/diagonal noise.
"""
struct SROCK1{interpretation,E} <: StochasticDiffEqAlgorithm
  eigen_est::E
end
SROCK1(;interpretation=SciMLBase.AlgorithmInterpretation.Ito,eigen_est=nothing) = SROCK1{interpretation,typeof(eigen_est)}(eigen_est)

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
"""
SROCKEM: S-ROCK Method
Is fixed step Euler-Mayurama with first order ROCK stabilization thus can handle stiff problems.
Only for Ito problems. Defaults to strong and weak order 1.0, but can solve with weak order 0.5 as SROCKEM(strong_order_1=false).
This method can handle 1-dimensional, diagonal and multi-dimensional noise.
"""
struct SROCKEM{E} <: StochasticDiffEqAlgorithm
  strong_order_1::Bool
  eigen_est::E
end
SROCKEM(;strong_order_1=true,eigen_est=nothing) = SROCKEM(strong_order_1,eigen_est)
"""
SKSROCK: S-ROCK Method
Is fixed step stabilized explicit method for stiff Ito problems.
Strong order 0.5 and weak order 1.
This method has a better stability domain then SROCK1.
Also it allows special post-processing techniques in case of ergodic dynamical systems, in the context of ergodic Brownian dynamics, to achieve order 2 accuracy.
SKSROCK(;post_processing=true) will make use of post processing.
By default it doesn't use post processing.
Post processing is optional and under development.
The rest of the method is completely functional and can handle 1-dimensional, diagonal and multi-dimensional noise.
"""
struct SKSROCK{E} <: StochasticDiffEqAlgorithm
  post_processing::Bool
  eigen_est::E
end
SKSROCK(;post_processing=false,eigen_est=nothing) = SKSROCK(post_processing,eigen_est)
"""
TangXiaoSROCK2: S-ROCK Method
Is a fixed step size stabilized expicit method for stiff problems.
Only for Ito problems. Weak order of 2 and strog order of 1.
Has 5 versions with different stability domains which can be used as TangXiaoSROCK2(version_num=i) where i is 1-5. Under Development.
"""
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

SRA: Nonstiff Method
Adaptive strong order 1.5 methods for additive Ito and Stratonovich SDEs.
Default tableau is for SRA1. Can handle diagonal, non-diagonal and scalar additive noise.
"""
struct SRA{TabType} <: StochasticDiffEqAdaptiveAlgorithm
  tableau::TabType
end
SRA(;tableau=constructSRA1()) = SRA(tableau)

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X

SRI: Nonstiff Method
Adaptive strong order 1.5 methods for diagonal/scalar Ito SDEs.
Default tableau is for SRIW1.
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

SRIW1: Nonstiff Method
Adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs.
"""
struct SRIW1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X

SRIW2: Nonstiff Method
Adaptive strong order 1.5 and weak order 3.0 for diagonal/scalar Ito SDEs.
"""
struct SRIW2 <: StochasticDiffEqAdaptiveAlgorithm end
"""
SOSRI: Nonstiff Method
Stability-optimized adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs.
Stable at high tolerances and robust to stiffness.
"""
struct SOSRI <: StochasticDiffEqAdaptiveAlgorithm end
"""
SOSRI2: Nonstiff Method
Stability-optimized adaptive strong order 1.5 and weak order 2.0 for diagonal/scalar Ito SDEs.
Stable at high tolerances and robust to stiffness.
"""
struct SOSRI2 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X

SRA1: Nonstiff Method
Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2.
Can handle diagonal, non-diagonal, and scalar additive noise.
"""
struct SRA1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X

SRA2: Nonstiff Method
Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2.
Can handle diagonal, non-diagonal, and scalar additive noise.
"""
struct SRA2 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Runge–Kutta Methods for the Strong Approximation of Solutions of
Stochastic Differential Equations, SIAM J. Numer. Anal., 48 (3), pp. 922–952.
DOI:10.1137/09076636X

SRA3: Nonstiff Method
Adaptive strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 3.
Can handle non-diagonal and scalar additive noise.
"""
struct SRA3 <: StochasticDiffEqAdaptiveAlgorithm end
"""
SOSRA: Nonstiff Method
A stability-optimized adaptive SRA. Strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2.
Can handle diagonal, non-diagonal, and scalar additive noise.
Stable at high tolerances and robust to stiffness.
"""
struct SOSRA <: StochasticDiffEqAdaptiveAlgorithm end
"""
SOSRA2: Nonstiff Method
A stability-optimized adaptive SRA. Strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2.
Can handle diagonal, non-diagonal, and scalar additive noise.
Stable at high tolerances and robust to stiffness.
"""
struct SOSRA2 <: StochasticDiffEqAdaptiveAlgorithm end

################################################################################

# Rossler second order for weak approx.

"""
Debrabant, K. and Rößler A., Families of efficient second order Runge–Kutta methods
for the weak approximation of Itô stochastic differential equations,
Applied Numerical Mathematics 59, pp. 582–594 (2009)
DOI:10.1016/j.apnum.2008.03.012

DRI1: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs with minimized error constants (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct DRI1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Debrabant, K. and Rößler A., Families of efficient second order Runge–Kutta methods
for the weak approximation of Itô stochastic differential equations,
Applied Numerical Mathematics 59, pp. 582–594 (2009)
DOI:10.1016/j.apnum.2008.03.012

DRI1NM: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs with minimized error constants (deterministic order 3).
Can handle non-mixing diagonal (i.e., du[k] = f(u[k])) and scalar additive noise.
"""
struct DRI1NM <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations,
SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009)
DOI:10.1137/060673308

RI1: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RI1 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations,
SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009)
DOI:10.1137/060673308

RI3: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RI3 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations,
SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009)
DOI:10.1137/060673308

RI5: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RI5 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Rößler A., Second Order Runge–Kutta Methods for Itô Stochastic Differential Equations,
SIAM J. Numer. Anal., 47, pp. 1713-1738 (2009)
DOI:10.1137/060673308

RI6: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RI6 <: StochasticDiffEqAdaptiveAlgorithm end

"""
Debrabant, K. and Rößler A., Classification of Stochastic Runge–Kutta Methods for
the Weak Approximation of Stochastic Differential Equations,
Mathematics and Computers in Simulation 77, pp. 408-420 (2008)
DOI:10.1016/j.matcom.2007.04.016

RDI1WM: High Weak Order Method

Fixed step weak order 1.0 for Ito SDEs (deterministic order 2).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RDI1WM <: StochasticDiffEqAlgorithm end

"""
Debrabant, K. and Rößler A., Classification of Stochastic Runge–Kutta Methods for
the Weak Approximation of Stochastic Differential Equations,
Mathematics and Computers in Simulation 77, pp. 408-420 (2008)
DOI:10.1016/j.matcom.2007.04.016

RDI2WM: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RDI2WM <: StochasticDiffEqAdaptiveAlgorithm end

"""
Debrabant, K. and Rößler A., Classification of Stochastic Runge–Kutta Methods for
the Weak Approximation of Stochastic Differential Equations,
Mathematics and Computers in Simulation 77, pp. 408-420 (2008)
DOI:10.1016/j.matcom.2007.04.016

RDI3WM: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.†
"""
struct RDI3WM <: StochasticDiffEqAdaptiveAlgorithm end

"""
Debrabant, K. and Rößler A., Classification of Stochastic Runge–Kutta Methods for
the Weak Approximation of Stochastic Differential Equations,
Mathematics and Computers in Simulation 77, pp. 408-420 (2008)
DOI:10.1016/j.matcom.2007.04.016

RDI4WM: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RDI4WM <: StochasticDiffEqAdaptiveAlgorithm end


"""
Tang, X., & Xiao, A., Efficient weak second-order stochastic Runge–Kutta methods
for Itô stochastic differential equations,
BIT Numerical Mathematics, 57, 241-260 (2017)
DOI: 10.1007/s10543-016-0618-9

W2Ito1: High Weak Order Method
Adaptive step weak order 2.0 for Ito SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct W2Ito1 <: StochasticDiffEqAdaptiveAlgorithm end

# Stratonovich sense

"""
Rößler A., Second order Runge–Kutta methods for Stratonovich stochastic differential
equations, BIT Numerical Mathematics 47, pp. 657-680 (2007)
DOI:10.1007/s10543-007-0130-3

RS1: High Weak Order Method
Fixed step weak order 2.0 for Stratonovich SDEs (deterministic order 2).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RS1 <: StochasticDiffEqAlgorithm end

"""
Rößler A., Second order Runge–Kutta methods for Stratonovich stochastic differential
equations, BIT Numerical Mathematics 47, pp. 657-680 (2007)
DOI:10.1007/s10543-007-0130-3

RS2: High Weak Order Method
Fixed step weak order 2.0 for Stratonovich SDEs (deterministic order 3).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct RS2 <: StochasticDiffEqAlgorithm end

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)

PL1WM: High Weak Order Method
Fixed step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct PL1WM <: StochasticDiffEqAlgorithm end

"""
Kloeden, P.E., Platen, E., Numerical Solution of Stochastic Differential Equations.
Springer. Berlin Heidelberg (2011)

PL1WMA: High Weak Order Method
Fixed step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle additive noise.
"""
struct PL1WMA <: StochasticDiffEqAlgorithm end

"""
Komori, Y., Weak second-order stochastic Runge–Kutta methods for non-commutative
stochastic differential equations, Journal of Computational and Applied
Mathematics 206, pp. 158 – 173 (2007)
DOI:10.1016/j.cam.2006.06.006

NON: High Weak Order Method
Fixed step weak order 2.0 for Stratonovich SDEs (deterministic order 4).
Can handle diagonal, non-diagonal, non-commuting, and scalar additive noise.
"""
struct NON <: StochasticDiffEqAlgorithm end

"""
Komori, Y., Weak order stochastic Runge–Kutta methods for commutative stochastic
differential equations, Journal of Computational and Applied Mathematics 203,
pp. 57 – 79 (2007)
DOI:10.1016/j.cam.2006.03.010
"""
struct COM <: StochasticDiffEqAlgorithm end

"""
Komori, Y., & Burrage, K. (2011). Supplement: Efficient weak second order stochastic
Runge–Kutta methods for non-commutative Stratonovich stochastic differential equations.
Journal of computational and applied mathematics, 235(17), pp. 5326 - 5329 (2011)
DOI:10.1016/j.cam.2011.04.021
"""
struct NON2 <: StochasticDiffEqAlgorithm end


"""
Tocino, A. and Vigo-Aguiar, J., Weak Second Order Conditions for Stochastic Runge-
Kutta Methods, SIAM Journal on Scientific Computing 24, pp. 507 - 523 (2002)
DOI:10.1137/S1064827501387814

SIEA:High Weak Order Method
Fixed step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal and scalar additive noise.†
Stochastic generalization of the improved Euler method.
"""
struct SIEA <: StochasticDiffEqAlgorithm end

"""
Tocino, A. and Vigo-Aguiar, J., Weak Second Order Conditions for Stochastic Runge-
Kutta Methods, SIAM Journal on Scientific Computing 24, pp. 507 - 523 (2002)
DOI:10.1137/S1064827501387814

SMEA: High Weak Order Method
Fixed step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal and scalar additive noise.†
Stochastic generalization of the modified Euler method.
"""
struct SMEA <: StochasticDiffEqAlgorithm end

"""
Tocino, A. and Vigo-Aguiar, J., Weak Second Order Conditions for Stochastic Runge-
Kutta Methods, SIAM Journal on Scientific Computing 24, pp. 507 - 523 (2002)
DOI:10.1137/S1064827501387814

SIEB: High Weak Order Method
Fixed step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal and scalar additive noise.†
Stochastic generalization of the improved Euler method.
"""
struct SIEB <: StochasticDiffEqAlgorithm end

"""
Tocino, A. and Vigo-Aguiar, J., Weak Second Order Conditions for Stochastic Runge-
Kutta Methods, SIAM Journal on Scientific Computing 24, pp. 507 - 523 (2002)
DOI:10.1137/S1064827501387814

SMEB: High Order Weak Method
Fixed step weak order 2.0 for Ito SDEs (deterministic order 2).
Can handle diagonal and scalar additive noise.†
Stochastic generalization of the modified Euler method.
"""
struct SMEB <: StochasticDiffEqAlgorithm end


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
"""
ImplicitEM: Stiff Method
An order 0.5 Ito drift-implicit method.
This is a theta method which defaults to theta=1 or the Trapezoid method on the drift term.
This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.
"""
struct ImplicitEM{CS,AD,F,F2,P,FDT,ST,CJ,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller}
  linsolve::F
  nlsolve::F2
  precs::P
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ImplicitEM(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          standardtag = Val{true}(),concrete_jac = nothing,
                          precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                          linsolve=nothing,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          theta = 1,symplectic=false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEM{chunk_size,autodiff,
                          typeof(linsolve),typeof(nlsolve),typeof(precs),diff_type,
                          SciMLBase._unwrap_val(standardtag),
                          SciMLBase._unwrap_val(concrete_jac),
                          typeof(new_jac_conv_bound),controller}(
                          linsolve,nlsolve,precs,
                          symplectic ? 1/2 : theta,
                          extrapolant,new_jac_conv_bound,symplectic)

STrapezoid(;kwargs...) = ImplicitEM(;theta=1/2,kwargs...)
SImplicitMidpoint(;kwargs...) = ImplicitEM(;theta=1/2,symplectic=true,kwargs...)
"""
ImplicitEulerHeun: Stiff Method
An order 0.5 Stratonovich drift-implicit method.
This is a theta method which defaults to theta=1/2 or the Trapezoid method on the drift term.
This method defaults to symplectic=false, but when true and theta=1 this is the implicit Midpoint method on the drift term and is symplectic in distribution.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.
"""
struct ImplicitEulerHeun{CS,AD,F,P,FDT,ST,CJ,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller}
  linsolve::F
  nlsolve::N
  precs::P
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ImplicitEulerHeun(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          standardtag = Val{true}(),concrete_jac = nothing,
                          precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                          linsolve=nothing,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          theta = 1,symplectic = false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive) =
                          ImplicitEulerHeun{chunk_size,autodiff,
                          typeof(linsolve),typeof(precs),diff_type,
                          SciMLBase._unwrap_val(standardtag),
                          SciMLBase._unwrap_val(concrete_jac),
                          typeof(nlsolve),
                          typeof(new_jac_conv_bound),controller}(
                          linsolve,nlsolve,precs,
                          symplectic ? 1/2 : theta,
                          extrapolant,
                          new_jac_conv_bound,symplectic)

"""
ImplicitRKMil: Stiff Method
An order 1.0 drift-implicit method.
This is a theta method which defaults to theta=1 or the Trapezoid method on the drift term.
Defaults to solving the Ito problem, but ImplicitRKMil(interpretation=SciMLBase.AlgorithmInterpretation.Stratonovich) makes it solve the Stratonovich problem.
This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution.
Handles diagonal and scalar noise. Uses a 1.5/2.0 heuristic for adaptive time stepping.
"""
struct ImplicitRKMil{CS,AD,F,P,FDT,ST,CJ,N,T2,Controller,interpretation} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller}
  linsolve::F
  nlsolve::N
  precs::P
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ImplicitRKMil(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                          standardtag = Val{true}(),concrete_jac = nothing,
                          precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                          linsolve=nothing,nlsolve=NLNewton(),
                          extrapolant=:constant,
                          theta = 1,symplectic = false,
                          new_jac_conv_bound = 1e-3,
                          controller = :Predictive,interpretation=SciMLBase.AlgorithmInterpretation.Ito) =
                          ImplicitRKMil{chunk_size,autodiff,
                          typeof(linsolve),typeof(precs),diff_type,
                          SciMLBase._unwrap_val(standardtag),
                          SciMLBase._unwrap_val(concrete_jac),
                          typeof(nlsolve),typeof(new_jac_conv_bound),
                          controller,interpretation}(
                          linsolve,nlsolve,precs,
                          symplectic ? 1/2 : theta,
                          extrapolant,
                          new_jac_conv_bound,symplectic)
"""
ISSEM: Stiff Method
An order 0.5 split-step Ito implicit method.
It is fully implicit, meaning it can handle stiffness in the noise term.
This is a theta method which defaults to theta=1 or the Trapezoid method on the drift term.
This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution.
Can handle all forms of noise, including non-diagonal, scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.
"""
struct ISSEM{CS,AD,F,P,FDT,ST,CJ,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller}
  linsolve::F
  nlsolve::N
  precs::P
  theta::T2
  extrapolant::Symbol
  new_jac_conv_bound::T2
  symplectic::Bool
end
ISSEM(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                       standardtag = Val{true}(),concrete_jac = nothing,
                       precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                       linsolve=nothing,nlsolve=NLNewton(),
                       extrapolant=:constant,
                       theta = 1,symplectic=false,
                       new_jac_conv_bound = 1e-3,
                       controller = :Predictive) =
                       ISSEM{chunk_size,autodiff,
                       typeof(linsolve),typeof(precs),diff_type,
                       SciMLBase._unwrap_val(standardtag),
                       SciMLBase._unwrap_val(concrete_jac),
                       typeof(nlsolve),
                       typeof(new_jac_conv_bound),controller}(
                       linsolve,nlsolve,precs,
                       symplectic ? 1/2 : theta,
                       extrapolant,
                       new_jac_conv_bound,symplectic)
"""
ISSEulerHeun: Stiff Method
An order 0.5 split-step Stratonovich implicit method.
It is fully implicit, meaning it can handle stiffness in the noise term.
This is a theta method which defaults to theta=1 or the Trapezoid method on the drift term.
This method defaults to symplectic=false, but when true and theta=1/2 this is the implicit Midpoint method on the drift term and is symplectic in distribution.
Can handle all forms of noise, including non-diagonal,Q scalar, and colored noise. Uses a 1.0/1.5 heuristic for adaptive time stepping.
"""
struct ISSEulerHeun{CS,AD,F,P,FDT,ST,CJ,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller}
 linsolve::F
 nlsolve::N
 precs::P
 theta::T2
 extrapolant::Symbol
 new_jac_conv_bound::T2
 symplectic::Bool
end
ISSEulerHeun(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                      standardtag = Val{true}(),concrete_jac = nothing,
                      precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                      linsolve=nothing,nlsolve=NLNewton(),
                      extrapolant=:constant,
                      theta = 1,symplectic=false,
                      new_jac_conv_bound = 1e-3,
                      controller = :Predictive) =
                      ISSEulerHeun{chunk_size,autodiff,
                      typeof(linsolve),typeof(precs),diff_type,
                      SciMLBase._unwrap_val(standardtag),
                      SciMLBase._unwrap_val(concrete_jac),
                      typeof(nlsolve),typeof(new_jac_conv_bound),controller}(
                      linsolve,nlsolve,precs,
                      symplectic ? 1/2 : theta,
                      extrapolant,
                      new_jac_conv_bound,symplectic)
"""
SKenCarp: Stiff Method
Adaptive L-stable drift-implicit strong order 1.5 for additive Ito and Stratonovich SDEs with weak order 2.
Can handle diagonal, non-diagonal and scalar additive noise.
"""
struct SKenCarp{CS,AD,F,P,FDT,ST,CJ,N,T2,Controller} <: StochasticDiffEqNewtonAdaptiveAlgorithm{CS,AD,FDT,ST,CJ,Controller}
  linsolve::F
  nlsolve::N
  precs::P
  smooth_est::Bool
  extrapolant::Symbol
  new_jac_conv_bound::T2
  ode_error_est::Bool
end

SKenCarp(;chunk_size=0,autodiff=true,diff_type=Val{:central},
                   standardtag = Val{true}(),concrete_jac = nothing,
                   precs = OrdinaryDiffEqCore.DEFAULT_PRECS,
                   linsolve=nothing,nlsolve=NLNewton(),
                   smooth_est=true,extrapolant=:min_correct,
                   new_jac_conv_bound = 1e-3,controller = :Predictive,
                   ode_error_est = true) =
 SKenCarp{chunk_size,autodiff,typeof(linsolve),typeof(precs),diff_type,
        SciMLBase._unwrap_val(standardtag),SciMLBase._unwrap_val(concrete_jac),
        typeof(nlsolve),typeof(new_jac_conv_bound),controller}(
        linsolve,nlsolve,precs,smooth_est,extrapolant,new_jac_conv_bound,
        ode_error_est)


################################################################################

# Jumps

function TauLeaping_docstring(
        description::String,
        name::String;
        references::String = "",
        extra_keyword_description::String = "",
        extra_keyword_default::String = "")
    keyword_default = """
        adaptive = true,
        """ * "\n" * extra_keyword_default

    keyword_default_description = """
    - `adaptive`: Boolean to enable/disable adaptive step sizing. When `true`, the step size `τ` is adjusted dynamically based on error estimates or bounds. Defaults to `true`.
    """ * "\n" * extra_keyword_description

    docstring = """
    $description

    ### Algorithm Type
    Stochastic Jump Method

    ### References
    $references

    ### Keyword Arguments
    $keyword_default_description

    ### Default Values
    $keyword_default
    """
    return docstring
end

@doc TauLeaping_docstring(
    "An explicit tau-leaping method for stochastic jump processes with optional post-leap step size adaptivity. " *
    "This algorithm approximates the stochastic simulation algorithm (SSA) by advancing the system state over " *
    "a fixed time step `τ` using Poisson-distributed jump counts based on initial propensities. When `adaptive=true`, " *
    "it adjusts `τ` dynamically based on post-leap error estimates derived from propensity changes.",
    "TauLeaping",
    references = """@article{gillespie2001approximate,
    title={Approximate accelerated stochastic simulation of chemically reacting systems},
    author={Gillespie, Daniel T},
    journal={The Journal of Chemical Physics},
    volume={115},
    number={4},
    pages={1716--1733},
    year={2001},
    publisher={AIP Publishing}}""",
    extra_keyword_description = """
    - `dtmax`: Maximum allowed step size.
    - `dtmin`: Minimum allowed step size.
    - `qmax`: Maximum step size increase factor.
    - `qmin`: Minimum step size reduction factor.
    - `gamma`: Safety factor for step size adjustment.
    """,
    extra_keyword_default = """
    dtmax = 10.0,
    dtmin = 1e-6
    """)
struct TauLeaping <: StochasticDiffEqJumpAdaptiveAlgorithm
    adaptive::Bool
end

function TauLeaping(; adaptive=true)
    TauLeaping(adaptive)
end

@doc TauLeaping_docstring(
    "An adaptive tau-leaping method for stochastic jump processes that selects the step size `τ` prior to each leap " *
    "based on bounds on the expected change in state variables. Introduced by Cao et al., this method ensures stability " *
    "and accuracy by constraining the relative change in propensities, controlled by the `epsilon` parameter. " *
    "When `adaptive=false`, a fixed step size is used.",
    "CaoTauLeaping",
    references = """@article{cao2006efficient,
    title={Efficient step size selection for the tau-leaping simulation method},
    author={Cao, Yang and Gillespie, Daniel T and Petzold, Linda R},
    journal={The Journal of Chemical Physics},
    volume={124},
    number={4},
    pages={044109},
    year={2006},
    publisher={AIP Publishing}}""",
    extra_keyword_description = """
    - `epsilon`: Tolerance parameter controlling the relative change in state variables for step size selection.
    - `dtmax`: Maximum allowed step size.
    - `dtmin`: Minimum allowed step size.
    """,
    extra_keyword_default = """
    epsilon = 0.03,
    dtmax = 10.0,
    dtmin = 1e-6
    """)
struct CaoTauLeaping <: StochasticDiffEqJumpAdaptiveAlgorithm
    adaptive::Bool
    epsilon::Float64
end

function CaoTauLeaping(; adaptive=true, epsilon=0.03)
    CaoTauLeaping(adaptive, epsilon)
end

################################################################################

# Etc.

struct StochasticCompositeAlgorithm{T,F} <: StochasticDiffEqCompositeAlgorithm
  algs::T
  choice_function::F
end

struct RandomEM <: StochasticDiffEqRODEAlgorithm end

struct RandomHeun <: StochasticDiffEqRODEAlgorithm end

struct RandomTamedEM <: StochasticDiffEqRODEAlgorithm end

const SplitSDEAlgorithms = Union{IIF1M,IIF2M,IIF1Mil,SKenCarp,SplitEM}

@doc raw"""
Leimkuhler B., Matthews C., Robust and efficient configurational molecular sampling via
Langevin dynamics, J. Chem. Phys. 138, 174102 (2013)
DOI:10.1063/1.4802990

```math
du = vdt \\
dv = f(v,u) dt - \gamma v dt + g(u) \sqrt{2\gamma} dW
```
"""
struct BAOAB{T} <: StochasticDiffEqAlgorithm
  gamma::T
  scale_noise::Bool
end
BAOAB(;gamma=1.0, scale_noise=true) = BAOAB(gamma, scale_noise)
