"""
`SDEProblem`

Wraps the data which defines an SDE problem

```math
u = f(u,t)dt + Σgᵢ(u,t)dWⁱ
```

with initial condition ``u0``.

### Constructors

`SDEProblem(f,g,u0;analytic=nothing)` : Defines the SDE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the SDE.
* `g`: The noise function in the SDE.
* `u0`: The initial condition.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system
* `sizeu`: The size of the initial condition (and thus `u`)
* `noise`: The noise process applied to the noise upon generation.

"""
type SDEProblem <: AbstractSDEProblem
  f::Function
  g::Function
  u0#::AbstractArray
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  sizeu#::Tuple
  isinplace::Bool
  noise::NoiseProcess
  function SDEProblem(f,g,u0;analytic=nothing,noise=WHITE_NOISE)
    isinplace = numparameters(f)>=3
    if analytic==nothing
      knownanalytic = false
      analytic=(t,u,W)->0
    else
      knownanalytic = true
    end
    if typeof(u0) <: Number
      sizeu = (1,)
      numvars = 1
    else
      sizeu = size(u0)
      numvars = size(u0)[end]
    end
    new(f,g,u0,analytic,knownanalytic,numvars,sizeu,isinplace,noise)
  end
end
