"""
RosslerSRI

Holds the Butcher tableaus for a Rosser SRI method.
"""
immutable RosslerSRI{T,T2} <: Tableau
  c₀::Vector{T2}
  c₁::Vector{T2}
  A₀::Matrix{T}
  A₁::Matrix{T}
  B₀::Matrix{T}
  B₁::Matrix{T}
  α::Vector{T}
  β₁::Vector{T}
  β₂::Vector{T}
  β₃::Vector{T}
  β₄::Vector{T}
  order::Rational{Int}
end

"""
RosslerSRA

Holds the Butcher tableaus for a Rosser SRA method.
"""
immutable RosslerSRA{T,T2} <: Tableau
  c₀::Vector{T2}
  c₁::Vector{T2}
  A₀::Matrix{T}
  B₀::Matrix{T}
  α::Vector{T}
  β₁::Vector{T}
  β₂::Vector{T}
  order::Rational{Int}
end

"""
constructSRIW1()

Constructs the tableau type for the SRIW1 method.
"""
function constructSRIW1(T=Float64,T2=Float64)
  c₀ = [0;3//4;0;0]
  c₁ = [0;1//4;1;1//4]
  A₀ = [0 0 0 0
      3//4 0 0 0
      0 0 0 0
      0 0 0 0]
  A₁ = [0 0 0 0
      1//4 0 0 0
      1 0 0 0
      0 0 1//4 0]
  B₀ = [0 0 0 0
      3//2 0 0 0
      0 0 0 0
      0 0 0 0]
  B₁ = [0 0 0 0
      1//2 0 0 0
      -1 0 0 0
      -5 3 1//2 0]

  α = [1//3;2//3;0;0]

  β₁ = [-1;4//3;2//3;0]
  β₂ = -[1;-4//3;1//3;0]
  β₃ = [2;-4//3;-2//3;0]
  β₄ = [-2;5//3;-2//3;1]
  RosslerSRI(map(T2,c₀),map(T2,c₁),
             map(T,A₀),map(T,A₁),
             map(T,B₀),map(T,B₁),
             map(T,α),map(T,β₁),map(T,β₂),
             map(T,β₃),map(T,β₄),3//2)
end

"""
constructSRA1()

Constructs the taleau type for the SRA1 method.
"""
function constructSRA1(T=Float64,T2=Float64)
  α  = [1//3;2//3]
  β₁ = [1;0]
  β₂ = [-1.0;1]
  A₀ = [0 0
       3//4 0]
  B₀ = [0 0
       3//2 0]
  c₀ = [0;3//4]
  c₁ = [1;0]
  RosslerSRA(map(T2,c₀),map(T2,c₁),
             map(T,A₀),map(T,B₀),
             map(T,α),map(T,β₁),map(T,β₂),4//2)
end

"""
checkSRIOrder(RosslerSRI)

Determines whether the order conditions are met via the tableaus of the SRI method.
"""
function checkSRIOrder(RosslerSRI;tol=1e-6)
  @unpack c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄ = RosslerSRI
  e = ones(size(α))
  conditions = Vector{Bool}(25)
  conditions[1] = abs(dot(α,e)-1)<tol
  conditions[2] = abs(dot(β₁,e)-1)<tol
  conditions[3] = abs(dot(β₂,e)-0)<tol
  conditions[4] = abs(dot(β₃,e)-0)<tol
  conditions[5] = abs(dot(β₄,e)-0)<tol
  conditions[6] = abs(sum(β₁'*B₁)-0)<tol
  conditions[7] = abs(sum(β₂'*B₁)-1)<tol
  conditions[8] = abs(sum(β₃'*B₁)-0)<tol
  conditions[9] = abs(sum(β₄'*B₁)-0)<tol
  conditions[10]= abs(sum(α'*A₀)- .5) < tol
  conditions[11]= abs(sum(α'*B₀)- 1) < tol
  conditions[12]= abs(sum(α'*(B₀*e).^2)- 3/2) < tol
  conditions[13]= abs(sum(β₁'*A₁)-1) < tol
  conditions[14]= abs(sum(β₂'*A₁)-0) < tol
  conditions[15]= abs(sum(β₃'*A₁)+1) < tol
  conditions[16]= abs(sum(β₄'*A₁)-0) < tol
  conditions[17]= abs(sum(β₁'*(B₁*e).^2)-1) < tol
  conditions[18]= abs(sum(β₂'*(B₁*e).^2)-0) < tol
  conditions[19]= abs(sum(β₃'*(B₁*e).^2)+1) < tol
  conditions[20]= abs(sum(β₄'*(B₁*e).^2)-2) < tol
  conditions[22]= abs(sum(β₂'*B₁*(B₁*e))-0) < tol
  conditions[21]= abs(sum(β₁'*B₁*(B₁*e))-0) < tol
  conditions[23]= abs(sum(β₃'*B₁*(B₁*e))-0) < tol
  conditions[24]= abs(sum(β₄'*B₁*(B₁*e))-1) < tol
  conditions[25]= abs.(.5*β₁'*(A₁*(B₀*e)) + (1/3)*β₃'*(A₁*(B₀*e)))[1] .< tol
  return(conditions)
end

"""
checkSRAOrder(RosslerSRI)

Determines whether the order conditions are met via the tableaus of the SRA method.
"""
function checkSRAOrder(SRA;tol=1e-6)
  @unpack c₀,c₁,A₀,B₀,α,β₁,β₂ = SRA
  e = ones(size(α))
  conditions = Vector{Bool}(8)
  conditions[1] = abs(dot(α,e)-1)<tol
  conditions[2] = abs(dot(β₁,e)-1)<tol
  conditions[3] = abs(dot(β₂,e)-0)<tol
  conditions[4] = abs(dot(α,B₀*e)-1).<tol
  conditions[5] = abs(dot(α,A₀*e)-1/2).<tol
  conditions[6] = abs(dot(α,(B₀*e).^2) - 3/2).<tol
  conditions[7] = abs(dot(β₁,c₁)-1)<tol
  conditions[8] = abs(dot(β₂,c₁)+1)<tol
  return(conditions)
end
