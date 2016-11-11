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
