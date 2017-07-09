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
constructSRIOpt1()

Opti6-12-11-10-01-47
"""
function constructSRIOpt1(T=Float64,T2=Float64)
  A₀ = [0.0 0.0 0.0 0.0; -0.04199224421316468 0.0 0.0 0.0; 2.842612915017106 -2.0527723684000727 0.0 0.0; 4.338237071435815 -2.8895936137439793 2.3017575594644466 0.0]
  A₁ = [0.0 0.0 0.0 0.0; 0.26204282091330466 0.0 0.0 0.0; 0.20903646383505375 -0.1502377115150361 0.0 0.0; 0.05836595312746999 0.6149440396332373 0.08535117634046772 0.0]
  B₀ = [0.0 0.0 0.0 0.0; -0.21641093549612528 0.0 0.0 0.0; 1.5336352863679572 0.26066223492647056 0.0 0.0; -1.0536037558179159 1.7015284721089472 -0.20725685784180017 0.0]
  B₁ = [0.0 0.0 0.0 0.0; -0.5119011827621657 0.0 0.0 0.0; 2.67767339866713 -4.9395031322250995 0.0 0.0; 0.15580956238299215 3.2361551006624674 -1.4223118283355949 0.0]
  α  = [1.140099274172029,-0.6401334255743456,0.4736296532772559,0.026404498125060714]
  β₁ = [-1.8453464565104432,2.688764531100726,-0.2523866501071323,0.40896857551684956]
  β₂ = [0.4969658141589478,-0.5771202869753592,-0.12919702470322217,0.2093514975196336]
  β₃ = [2.8453464565104425,-2.688764531100725,0.2523866501071322,-0.40896857551684945]
  β₄ = [0.11522663875443433,-0.57877086147738,0.2857851028163886,0.17775911990655704]
  e = ones(size(α))
  c₀ = A₀*e
  c₁ = A₁*e
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
  β₂ = [-1;1]
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
constructSRA2()

Constructs the taleau type for the SRA2 method.
"""
function constructSRA2(T=Float64,T2=Float64)
  α  = [1//3;2//3]
  β₁ = [0;1]
  β₂ = [3//2;-3//2]
  A₀ = [0 0
       3//4 0]
  B₀ = [0 0
       3//2 0]
  c₀ = [0;3//4]
  c₁ = [1//3;1]
  RosslerSRA(map(T2,c₀),map(T2,c₁),
             map(T,A₀),map(T,B₀),
             map(T,α),map(T,β₁),map(T,β₂),4//2)
end

"""
constructSRA3()

Constructs the taleau type for the SRA3 method.
"""
function constructSRA3(T=Float64,T2=Float64)
  α  = [1//6;1//6;2//3]
  β₁ = [1;0;0]
  β₂ = [-1;1;0]
  A₀ = [0 0 0
        1 0 0
        1//4 1//4 0]
  B₀ = [0 0 0
        0 0 0
        1 1//2 0]
  c₀ = [0;1;1//2]
  c₁ = [1;0;0]
  RosslerSRA(map(T2,c₀),map(T2,c₁),
             map(T,A₀),map(T,B₀),
             map(T,α),map(T,β₁),map(T,β₂),4//2)
end

"""
constructSOSRA()

Constructs the taleau type for the SOSRA method.
"""
function constructSOSRA(T=Float64,T2=Float64)
  a1 = 0.2889874966892885;
  a2 = 0.6859880440839937;
  a3 = 0.025024459226717772;
  c01 = 0;
  c02 = 0.6923962376159507;
  c03 = 1;
  c11 = 0;
  c12 = 0.041248171110700504;
  c13 = 1;
  b11 = -16.792534242221663;
  b12 = 17.514995785380226;
  b13 = 0.27753845684143835;
  b21 = 0.4237535769069274;
  b22 = 0.6010381474428539;
  b23 = -1.0247917243497813;
  A021 = 0.6923962376159507;
  A031 = -3.1609142252828395;
  A032 = 4.1609142252828395;
  B021 = 1.3371632704399763;
  B031 = 1.442371048468624;
  B032 = 1.8632741501139225;
  α  = [a1;a2;a3]
  β₁ = [b11;b12;b13]
  β₂ = [b21;b22;b23]
  A₀ = [0 0 0
        A021 0 0
        A031 A032 0]
  B₀ = [0 0 0
        B021 0 0
        B031 B032 0]
  c₀ = [c01;c02;c03]
  c₁ = [c11;c12;c13]
  RosslerSRA(map(T2,c₀),map(T2,c₁),
             map(T,A₀),map(T,B₀),
             map(T,α),map(T,β₁),map(T,β₂),4//2)
end

"""
constructSOSRA2()

Constructs the taleau type for the SOSRA method.
"""
function constructSOSRA2(T=Float64,T2=Float64)
  a1 = 0.4999999999999998;
  a2 = -0.9683897375354181;
  a3 = 1.4683897375354185;
  c01 = 0;
  c02 = 1;
  c03 = 1;
  c11 = 0;
  c12 = 1.0;
  c13 = 1;
  b11 = 0.0;
  b12 = 0.92438032145683;
  b13 = 0.07561967854316998;
  b21 = 1.0;
  b22 = -0.8169981105823436;
  b23 = -0.18300188941765633;
  A021 = 1;
  A031 = 0.9511849235504364;
  A032 = 0.04881507644956362;
  B021 = 0.7686101171003622;
  B031 = 0.43886792994934987;
  B032 = 0.7490415909204886;
  α  = [a1;a2;a3]
  β₁ = [b11;b12;b13]
  β₂ = [b21;b22;b23]
  A₀ = [0 0 0
        A021 0 0
        A031 A032 0]
  B₀ = [0 0 0
        B021 0 0
        B031 B032 0]
  c₀ = [c01;c02;c03]
  c₁ = [c11;c12;c13]
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
