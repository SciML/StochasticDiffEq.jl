"""
RosslerSRI

Holds the Butcher tableaus for a Roessler SRI method.
"""
struct RosslerSRI{T, T2} <: Tableau
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
struct RosslerSRA{T, T2} <: Tableau
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
function constructSRIW1(T = Float64, T2 = Float64)
    c₀ = [0; 3 // 4; 0; 0]
    c₁ = [0; 1 // 4; 1; 1 // 4]
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

    α = [1 // 3; 2 // 3; 0; 0]

    β₁ = [-1; 4 // 3; 2 // 3; 0]
    β₂ = -[1; -4 // 3; 1 // 3; 0]
    β₃ = [2; -4 // 3; -2 // 3; 0]
    β₄ = [-2; 5 // 3; -2 // 3; 1]
    RosslerSRI(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, A₁),
        map(T, B₀), map(T, B₁),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄), 3 // 2)
end

"""
constructSRIW2()

Constructs the tableau type for the SRIW1 method.
"""
function constructSRIW2(T = Float64, T2 = Float64)
    c₀ = [0; 1; 1 // 2; 0]
    c₁ = [0; 1 // 4; 1; 1 // 4]
    A₀ = [0 0 0 0
          1 0 0 0
          1//4 1//4 0 0
          0 0 0 0]
    A₁ = [0 0 0 0
          1//4 0 0 0
          1 0 0 0
          0 0 1//4 0]
    B₀ = [0 0 0 0
          0 0 0 0
          1 1//2 0 0
          0 0 0 0]
    B₁ = [0 0 0 0
          -1//2 0 0 0
          1 0 0 0
          2 -1 1//2 0]

    α = [1 // 6; 1 // 6; 2 // 3; 0]

    β₁ = [-1; 4 // 3; 2 // 3; 0]
    β₂ = [1; -4 // 3; 1 // 3; 0]
    β₃ = [2; -4 // 3; -2 // 3; 0]
    β₄ = [-2; 5 // 3; -2 // 3; 1]
    RosslerSRI(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, A₁),
        map(T, B₀), map(T, B₁),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄), 3 // 2)
end

"""
constructSRIOpt1()

Opti6-12-11-10-01-47
"""
function constructSRIOpt1(T = Float64, T2 = Float64)
    A₀ = [0.0 0.0 0.0 0.0; -0.04199224421316468 0.0 0.0 0.0;
          2.842612915017106 -2.0527723684000727 0.0 0.0;
          4.338237071435815 -2.8895936137439793 2.3017575594644466 0.0]
    A₁ = [0.0 0.0 0.0 0.0; 0.26204282091330466 0.0 0.0 0.0;
          0.20903646383505375 -0.1502377115150361 0.0 0.0;
          0.05836595312746999 0.6149440396332373 0.08535117634046772 0.0]
    B₀ = [0.0 0.0 0.0 0.0; -0.21641093549612528 0.0 0.0 0.0;
          1.5336352863679572 0.26066223492647056 0.0 0.0;
          -1.0536037558179159 1.7015284721089472 -0.20725685784180017 0.0]
    B₁ = [0.0 0.0 0.0 0.0; -0.5119011827621657 0.0 0.0 0.0;
          2.67767339866713 -4.9395031322250995 0.0 0.0;
          0.15580956238299215 3.2361551006624674 -1.4223118283355949 0.0]
    α = [1.140099274172029, -0.6401334255743456, 0.4736296532772559, 0.026404498125060714]
    β₁ = [-1.8453464565104432, 2.688764531100726, -0.2523866501071323, 0.40896857551684956]
    β₂ = [0.4969658141589478, -0.5771202869753592, -0.12919702470322217, 0.2093514975196336]
    β₃ = [2.8453464565104425, -2.688764531100725, 0.2523866501071322, -0.40896857551684945]
    β₄ = [0.11522663875443433, -0.57877086147738, 0.2857851028163886, 0.17775911990655704]
    e = ones(size(α))
    c₀ = A₀ * e
    c₁ = A₁ * e
    RosslerSRI(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, A₁),
        map(T, B₀), map(T, B₁),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄), 3 // 2)
end

"""
constructSRIOpt2()

Opti6-12-11-10-01-47
"""
function constructSRIOpt2(T = Float64, T2 = Float64)
    c0 = [0.0, 0.13804532298278663, 0.9999999999999992, 0.9999999999999994]
    c1 = [0.0, 0.45605532163856893, 0.999999999999996, 0.9999999999999962]
    A0 = [0.0 0.0 0.0 0.0
          0.13804532298278663 0.0 0.0 0.0
          0.5818361298250374 0.4181638701749618 0.0 0.0
          0.4670018408674211 0.8046204792187386 -0.27162232008616016 0.0]
    A1 = [0.0 0.0 0.0 0.0
          0.45605532163856893 0.0 0.0 0.0
          0.7555807846451692 0.24441921535482677 0.0 0.0
          0.6981181143266059 0.3453277086024727 -0.04344582292908241 0.0]
    B0 = [0.0 0.0 0.0 0.0
          0.08852381537667678 0.0 0.0 0.0
          1.0317752458971061 0.4563552922077882 0.0 0.0
          1.73078280444124 -0.46089678470929774 -0.9637509618944188 0.0]
    B1 = [0.0 0.0 0.0 0.0
          0.6753186815412179 0.0 0.0 0.0
          -0.07452812525785148 -0.49783736486149366 0.0 0.0
          -0.5591906709928903 0.022696571806569924 -0.8984927888368557 0.0]
    α = [-0.15036858140642623, 0.7545275856696072, 0.686995463807979, -0.2911544680711602]
    β1 = [-0.45315689727309133, 0.8330937231303951, 0.3792843195533544, 0.24077885458934192]
    β2 = [
        -0.4994383733810986, 0.9181786186154077, -0.25613778661003145, -0.16260245862427797]
    β3 = [
        1.4531568972730915, -0.8330937231303933, -0.3792843195533583, -0.24077885458934023]
    β4 = [-0.4976090683622265, 0.9148155835648892, -1.4102107084476505, 0.9930041932449877]
    RosslerSRI(map(T2, c0), map(T2, c1),
        map(T, A0), map(T, A1),
        map(T, B0), map(T, B1),
        map(T, α), map(T, β1), map(T, β2),
        map(T, β3), map(T, β4), 3 // 2)
end

"""
constructSRA1()

Constructs the taleau type for the SRA1 method.
"""
function constructSRA1(T = Float64, T2 = Float64)
    α = [1 // 3; 2 // 3]
    β₁ = [1; 0]
    β₂ = [-1; 1]
    A₀ = [0 0
          3//4 0]
    B₀ = [0 0
          3//2 0]
    c₀ = [0; 3 // 4]
    c₁ = [1; 0]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

"""
constructSRA2()

Constructs the taleau type for the SRA2 method.
"""
function constructSRA2(T = Float64, T2 = Float64)
    α = [1 // 3; 2 // 3]
    β₁ = [0; 1]
    β₂ = [3 // 2; -3 // 2]
    A₀ = [0 0
          3//4 0]
    B₀ = [0 0
          3//2 0]
    c₀ = [0; 3 // 4]
    c₁ = [1 // 3; 1]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

"""
constructSRA3()

Constructs the taleau type for the SRA3 method.
"""
function constructSRA3(T = Float64, T2 = Float64)
    α = [1 // 6; 1 // 6; 2 // 3]
    β₁ = [1; 0; 0]
    β₂ = [-1; 1; 0]
    A₀ = [0 0 0
          1 0 0
          1//4 1//4 0]
    B₀ = [0 0 0
          0 0 0
          1 1//2 0]
    c₀ = [0; 1; 1 // 2]
    c₁ = [1; 0; 0]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

"""
constructSOSRA()

Constructs the taleau type for the SOSRA method.
"""
function constructSOSRA(T = Float64, T2 = Float64)
    a1 = 0.2889874966892885
    a2 = 0.6859880440839937
    a3 = 0.025024459226717772
    c01 = 0
    c02 = 0.6923962376159507
    c03 = 1
    c11 = 0
    c12 = 0.041248171110700504
    c13 = 1
    b11 = -16.792534242221663
    b12 = 17.514995785380226
    b13 = 0.27753845684143835
    b21 = 0.4237535769069274
    b22 = 0.6010381474428539
    b23 = -1.0247917243497813
    A021 = 0.6923962376159507
    A031 = -3.1609142252828395
    A032 = 4.1609142252828395
    B021 = 1.3371632704399763
    B031 = 1.442371048468624
    B032 = 1.8632741501139225
    α = [a1; a2; a3]
    β₁ = [b11; b12; b13]
    β₂ = [b21; b22; b23]
    A₀ = [0 0 0
          A021 0 0
          A031 A032 0]
    B₀ = [0 0 0
          B021 0 0
          B031 B032 0]
    c₀ = [c01; c02; c03]
    c₁ = [c11; c12; c13]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

"""
constructSOSRA2()

Constructs the taleau type for the SOSRA method.
"""
function constructSOSRA2(T = Float64, T2 = Float64)
    a1 = 0.4999999999999998
    a2 = -0.9683897375354181
    a3 = 1.4683897375354185
    c01 = 0
    c02 = 1
    c03 = 1
    c11 = 0
    c12 = 1.0
    c13 = 1
    b11 = 0.0
    b12 = 0.92438032145683
    b13 = 0.07561967854316998
    b21 = 1.0
    b22 = -0.8169981105823436
    b23 = -0.18300188941765633
    A021 = 1
    A031 = 0.9511849235504364
    A032 = 0.04881507644956362
    B021 = 0.7686101171003622
    B031 = 0.43886792994934987
    B032 = 0.7490415909204886
    α = [a1; a2; a3]
    β₁ = [b11; b12; b13]
    β₂ = [b21; b22; b23]
    A₀ = [0 0 0
          A021 0 0
          A031 A032 0]
    B₀ = [0 0 0
          B021 0 0
          B031 B032 0]
    c₀ = [c01; c02; c03]
    c₁ = [c11; c12; c13]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

"""
checkSRIOrder(RosslerSRI)

Determines whether the order conditions are met via the tableaus of the SRI method.
"""
function checkSRIOrder(RosslerSRI; tol = 1e-6)
    @unpack c₀, c₁, A₀, A₁, B₀, B₁, α, β₁, β₂, β₃, β₄ = RosslerSRI
    e = ones(size(α))
    conditions = Vector{Bool}(undef, 25)
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₁, e) - 1) < tol
    conditions[3] = abs(dot(β₂, e) - 0) < tol
    conditions[4] = abs(dot(β₃, e) - 0) < tol
    conditions[5] = abs(dot(β₄, e) - 0) < tol
    conditions[6] = abs(sum(β₁' * B₁) - 0) < tol
    conditions[7] = abs(sum(β₂' * B₁) - 1) < tol
    conditions[8] = abs(sum(β₃' * B₁) - 0) < tol
    conditions[9] = abs(sum(β₄' * B₁) - 0) < tol
    conditions[10] = abs(sum(α' * A₀) - 0.5) < tol
    conditions[11] = abs(sum(α' * B₀) - 1) < tol
    conditions[12] = abs(sum(α' * (B₀ * e) .^ 2) - 3 / 2) < tol
    conditions[13] = abs(sum(β₁' * A₁) - 1) < tol
    conditions[14] = abs(sum(β₂' * A₁) - 0) < tol
    conditions[15] = abs(sum(β₃' * A₁) + 1) < tol
    conditions[16] = abs(sum(β₄' * A₁) - 0) < tol
    conditions[17] = abs(sum(β₁' * (B₁ * e) .^ 2) - 1) < tol
    conditions[18] = abs(sum(β₂' * (B₁ * e) .^ 2) - 0) < tol
    conditions[19] = abs(sum(β₃' * (B₁ * e) .^ 2) + 1) < tol
    conditions[20] = abs(sum(β₄' * (B₁ * e) .^ 2) - 2) < tol
    conditions[22] = abs(sum(β₂' * B₁ * (B₁ * e)) - 0) < tol
    conditions[21] = abs(sum(β₁' * B₁ * (B₁ * e)) - 0) < tol
    conditions[23] = abs(sum(β₃' * B₁ * (B₁ * e)) - 0) < tol
    conditions[24] = abs(sum(β₄' * B₁ * (B₁ * e)) - 1) < tol
    conditions[25] = abs.(0.5 * β₁' * (A₁ * (B₀ * e)) + (1 / 3) * β₃' * (A₁ * (B₀ * e)))[1] .<
                     tol
    return (conditions)
end

"""
checkSRAOrder(RosslerSRI)

Determines whether the order conditions are met via the tableaus of the SRA method.
"""
function checkSRAOrder(SRA; tol = 1e-6)
    @unpack c₀, c₁, A₀, B₀, α, β₁, β₂ = SRA
    e = ones(size(α))
    conditions = Vector{Bool}(undef, 8)
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₁, e) - 1) < tol
    conditions[3] = abs(dot(β₂, e) - 0) < tol
    conditions[4] = abs(dot(α, B₀ * e) - 1) .< tol
    conditions[5] = abs(dot(α, A₀ * e) - 1 / 2) .< tol
    conditions[6] = abs(dot(α, (B₀ * e) .^ 2) - 3 / 2) .< tol
    conditions[7] = abs(dot(β₁, c₁) - 1) < tol
    conditions[8] = abs(dot(β₂, c₁) + 1) < tol
    return (conditions)
end

"""
constructSOSRA2()

Constructs the taleau type for the SOSRA method.
"""
function constructSKenCarp(T = Float64, T2 = Float64)
    γ = convert(T, 0.435866521508459)
    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)
    # bhat1 = convert(T,2756255671327//12835298489170)
    # bhat2 = -convert(T,10771552573575//22201958757719)
    # bhat3 = convert(T,9247589265047//10645013368117)
    # bhat4 = convert(T,2193209047091//5459859503100)
    btilde1 = convert(T, 0.027099261876665316) # bhat1-a41
    btilde2 = convert(T, 0.11013520969201586) # bhat2-a42
    btilde3 = convert(T, -0.10306492520138458) # bhat3-a43
    btilde4 = convert(T, -0.0341695463672966) # bhat4-γ
    c3 = convert(T2, 0.6)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    nb021 = convert(T, -12.246764387585056)
    nb043 = convert(T, -14.432096958608753)
    α = [a41; a42; a43; γ]
    β₁ = [0, 0, 0, 1]
    β₂ = [1, 0, 0, -1]
    A₀ = [0 0 0 0
          γ γ 0 0
          a31 a32 γ 0
          a41 a42 a43 γ]
    B₀ = [0 0 0 0
          nb021 0 0 0
          0 0 0 0
          0 0 nb043 0]
    c₀ = [0, c2, c3, 1]
    c₁ = [0, 0, 0, 1]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

"""
constructSOSRA2()

Constructs the taleau type for the SOSRA method.
"""
function constructExplicitSKenCarp(T = Float64, T2 = Float64)
    γ = convert(T, 0.435866521508459)
    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)
    # bhat1 = convert(T,2756255671327//12835298489170)
    # bhat2 = -convert(T,10771552573575//22201958757719)
    # bhat3 = convert(T,9247589265047//10645013368117)
    # bhat4 = convert(T,2193209047091//5459859503100)
    btilde1 = convert(T, 0.027099261876665316) # bhat1-a41
    btilde2 = convert(T, 0.11013520969201586) # bhat2-a42
    btilde3 = convert(T, -0.10306492520138458) # bhat3-a43
    btilde4 = convert(T, -0.0341695463672966) # bhat4-γ
    c3 = convert(T2, 0.6)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    ea21 = convert(T, 1767732205903 // 2027836641118)
    ea31 = convert(T, 5535828885825 // 10492691773637)
    ea32 = convert(T, 788022342437 // 10882634858940)
    ea41 = convert(T, 6485989280629 // 16251701735622)
    ea42 = -convert(T, 4246266847089 // 9704473918619)
    ea43 = convert(T, 10755448449292 // 10357097424841)
    eb1 = convert(T, 1471266399579 // 7840856788654)
    eb2 = convert(T, -4482444167858 // 7529755066697)
    eb3 = convert(T, 11266239266428 // 11593286722821)
    eb4 = convert(T, 1767732205903 // 4055673282236)

    nb021 = convert(T, -12.246764387585056)
    nb043 = convert(T, -14.432096958608753)

    α = [eb1; eb2; eb3; eb4]
    β₁ = [0, 0, 0, 1]
    β₂ = [1, 0, 0, -1]
    A₀ = [0 0 0 0
          ea21 0 0 0
          ea31 ea32 0 0
          ea41 ea42 ea43 0]
    B₀ = [0 0 0 0
          nb021 0 0 0
          0 0 0 0
          0 0 nb043 0]
    c₀ = [0, c2, c3, 1]
    c₁ = [0, 0, 0, 1]
    RosslerSRA(map(T2, c₀), map(T2, c₁),
        map(T, A₀), map(T, B₀),
        map(T, α), map(T, β₁), map(T, β₂), 4 // 2)
end

struct SKenCarpTableau{T, T2}
    γ::T
    a31::T
    a32::T
    a41::T
    a42::T
    a43::T
    btilde1::T
    btilde2::T
    btilde3::T
    btilde4::T
    c3::T2
    α31::T
    α32::T
    α41::T
    α42::T
    ea21::T
    ea31::T
    ea32::T
    ea41::T
    ea42::T
    ea43::T
    eb1::T
    eb2::T
    eb3::T
    eb4::T
    ebtilde1::T
    ebtilde2::T
    ebtilde3::T
    ebtilde4::T
    nb021::T
    nb043::T
end

#=
# KenCarp3
# Predict z4 from Hermite z2 and z1
# Not z3 because c3 < c2 !

θ = c3/c2
dt = c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
θ = c4/c2
((1 + (-4θ + 3θ^2)) + (6θ*(1-θ)/c2)*γ)
((-2θ + 3θ^2) + (6θ*(1-θ)/c2)*γ)
=#
function SKenCarpTableau(
        ::Type{T}, ::Type{T2}) where {T <: CompiledFloats, T2 <: CompiledFloats}
    γ = convert(T, 0.435866521508459)
    a31 = convert(T, 0.2576482460664272)
    a32 = -convert(T, 0.09351476757488625)
    a41 = convert(T, 0.18764102434672383)
    a42 = -convert(T, 0.595297473576955)
    a43 = convert(T, 0.9717899277217721)
    # bhat1 = convert(T,2756255671327//12835298489170)
    # bhat2 = -convert(T,10771552573575//22201958757719)
    # bhat3 = convert(T,9247589265047//10645013368117)
    # bhat4 = convert(T,2193209047091//5459859503100)
    btilde1 = convert(T, 0.027099261876665316) # bhat1-a41
    btilde2 = convert(T, 0.11013520969201586) # bhat2-a42
    btilde3 = convert(T, -0.10306492520138458) # bhat3-a43
    btilde4 = convert(T, -0.0341695463672966) # bhat4-γ
    c3 = convert(T2, 0.6)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    ea21 = convert(T, 0.871733043016918)
    ea31 = convert(T, 0.5275890119763004)
    ea32 = convert(T, 0.0724109880236996)
    ea41 = convert(T, 0.3990960076760701)
    ea42 = -convert(T, 0.4375576546135194)
    ea43 = convert(T, 1.0384616469374492)
    eb1 = convert(T, 0.18764102434672383)
    eb2 = convert(T, -0.595297473576955)
    eb3 = convert(T, 0.9717899277217721)
    eb4 = convert(T, 0.435866521508459)
    ebtilde1 = convert(T, 0.027099261876665316)
    ebtilde2 = convert(T, 0.11013520969201586)
    ebtilde3 = -convert(T, 0.10306492520138458)
    ebtilde4 = -convert(T, 0.0341695463672966)

    nb021 = convert(T, -12.246764387585056)
    nb043 = convert(T, -14.432096958608753)
    SKenCarpTableau(
        γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31,
        α32, α41, α42, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4,
        ebtilde1, ebtilde2, ebtilde3, ebtilde4, nb021, nb043)
end

function SKenCarpTableau(T, T2)
    γ = convert(T, 1767732205903 // 4055673282236)
    a31 = convert(T, 2746238789719 // 10658868560708)
    a32 = -convert(T, 640167445237 // 6845629431997)
    a41 = convert(T, 1471266399579 // 7840856788654)
    a42 = -convert(T, 4482444167858 // 7529755066697)
    a43 = convert(T, 11266239266428 // 11593286722821)
    # bhat1 = convert(T,2756255671327//12835298489170)
    # bhat2 = -convert(T,10771552573575//22201958757719)
    # bhat3 = convert(T,9247589265047//10645013368117)
    # bhat4 = convert(T,2193209047091//5459859503100)
    btilde1 = convert(
        T, BigInt(681815649026867975666107) // BigInt(25159934323302256049469295)) # bhat1-a41
    btilde2 = convert(
        T, BigInt(18411887981491912264464127) // BigInt(167175311446532472108584143)) # bhat2-a42
    btilde3 = convert(
        T, BigInt(-12719313754959329011138489) // BigInt(123410692144842870217698057)) # bhat3-a43
    btilde4 = convert(
        T, BigInt(-47289384293135913063989) // BigInt(1383962894467812063558225)) # bhat4-γ
    c3 = convert(T2, 3 // 5)
    c2 = 2γ
    θ = c3 / c2
    α31 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α32 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)
    θ = 1 / c2
    α41 = ((1 + (-4θ + 3θ^2)) + (6θ * (1 - θ) / c2) * γ)
    α42 = ((-2θ + 3θ^2) + (6θ * (1 - θ) / c2) * γ)

    # Explicit Tableau
    ea21 = convert(T, 1767732205903 // 2027836641118)
    ea31 = convert(T, 5535828885825 // 10492691773637)
    ea32 = convert(T, 788022342437 // 10882634858940)
    ea41 = convert(T, 6485989280629 // 16251701735622)
    ea42 = -convert(T, 4246266847089 // 9704473918619)
    ea43 = convert(T, 10755448449292 // 10357097424841)
    eb1 = convert(T, 1471266399579 // 7840856788654)
    eb2 = convert(T, -4482444167858 // 7529755066697)
    eb3 = convert(T, 11266239266428 // 11593286722821)
    eb4 = convert(T, 1767732205903 // 4055673282236)
    ebtilde1 = convert(
        T, BigInt(681815649026867975666107) // BigInt(25159934323302256049469295))
    ebtilde2 = convert(
        T, BigInt(18411887981491912264464127) // BigInt(167175311446532472108584143))
    ebtilde3 = -convert(
        T, BigInt(12719313754959329011138489) // BigInt(123410692144842870217698057))
    ebtilde4 = -convert(
        T, BigInt(47289384293135913063989) // BigInt(1383962894467812063558225))

    # Noise Tableau

    nb021 = convert(T,
        parse(BigFloat,
            "-12.246764387585055918338744103409192607986567514699471403397969732723452087723101"))
    nb043 = convert(T,
        parse(BigFloat,
            "-14.432096958608752822047165680776748797565142459789556194474191884258734697161106"))
    SKenCarpTableau(
        γ, a31, a32, a41, a42, a43, btilde1, btilde2, btilde3, btilde4, c3, α31,
        α32, α41, α42, ea21, ea31, ea32, ea41, ea42, ea43, eb1, eb2, eb3, eb4,
        ebtilde1, ebtilde2, ebtilde3, ebtilde4, nb021, nb043)
end

# Flip them all!

# ebtilde1 = big(1471266399579)//7840856788654 - big(2756255671327)//12835298489170
# ebtilde2 = -big(4482444167858)//7529755066697 + big(10771552573575)//22201958757719
# ebtilde3 = big(11266239266428)//11593286722821 - big(9247589265047)//10645013368117
# ebtilde4 = big(1767732205903)//4055673282236 - big(2193209047091)//5459859503100

"""

RoesslerRI

Holds the Butcher tableaus for a Roessler RI method. (high weak order)

"""

struct RoesslerRI{T, T2} <: Tableau
    c₀::Vector{T2}
    c₁::Vector{T2}
    c₂::Vector{T2}
    A₀::Matrix{T}
    A₁::Matrix{T}
    A₂::Matrix{T}
    B₀::Matrix{T}
    B₁::Matrix{T}
    B₂::Matrix{T}
    α::Vector{T}
    β₁::Vector{T}
    β₂::Vector{T}
    β₃::Vector{T}
    β₄::Vector{T}
    order::Rational{Int}
    quantile::T
end

function constructDRI1(T = Float64, T2 = Float64)
    c₀ = [0; 1 // 2; 1]
    c₁ = [0; 342 // 491; 342 // 491]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1//2 0 0
          -1 2 0]
    A₁ = [0 0 0
          342//491 0 0
          342//491 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          (6 - sqrt(6))/10 0 0
          (3 + 2 * sqrt(6))/5 0 0]
    B₁ = [0 0 0
          3*sqrt(38 // 491) 0 0
          -3*sqrt(38 / 491) 0 0]
    B₂ = [0 0 0
          -214 // 513*sqrt(1105 // 991) -491 // 513*sqrt(221 // 4955) -491 // 513*sqrt(221 // 4955)
          214 // 513*sqrt(1105 // 991) 491 // 513*sqrt(221 // 4955) 491 // 513*sqrt(221 // 4955)]
    α = [1 // 6; 2 // 3; 1 // 6]

    β₁ = [193 // 684; 491 // 1368; 491 // 1368]
    β₂ = [0; 1 // 6 * sqrt(491 // 38); -1 // 6 * sqrt(491 // 38)]
    β₃ = [-4955 // 7072; 4955 // 14144; 4955 // 14144]
    β₄ = [0; -1 // 8 * sqrt(4955 // 221); 1 // 8 * sqrt(4955 // 221)]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRI1(T = Float64, T2 = Float64)
    c₀ = [0; 2 // 3; 2 // 3]
    c₁ = [0; 1; 1]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          2//3 0 0
          -1//3 1 0]
    A₁ = [0 0 0
          1 0 0
          1 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          1 0 0
          0 0 0]
    B₁ = [0 0 0
          1 0 0
          -1 0 0]
    B₂ = [0 0 0
          1 0 0
          -1 0 0]
    α = [1 // 4; 1 // 2; 1 // 4]

    β₁ = [1 // 2; 1 // 4; 1 // 4]
    β₂ = [0; 1 // 2; -1 // 2]
    β₃ = [-1 // 2; 1 // 4; 1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRI3(T = Float64, T2 = Float64)
    c₀ = [0; 1; 1 // 2]
    c₁ = [0; 1; 1]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1 0 0
          1//4 1//4 0]
    A₁ = [0 0 0
          1 0 0
          1 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          (3 - 2 * sqrt(6))/5 0 0
          (6 + sqrt(6))/10 0 0]
    B₁ = [0 0 0
          1 0 0
          -1 0 0]
    B₂ = [0 0 0
          1 0 0
          -1 0 0]
    α = [1 // 6; 1 // 6; 2 // 3]

    β₁ = [1 // 2; 1 // 4; 1 // 4]
    β₂ = [0; 1 // 2; -1 // 2]
    β₃ = [-1 // 2; 1 // 4; 1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRI5(T = Float64, T2 = Float64)
    c₀ = [0; 1; 5 // 12]
    c₁ = [0; 1 // 4; 1 // 4]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1 0 0
          25//144 35//144 0]
    A₁ = [0 0 0
          1//4 0 0
          1//4 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          1//3 0 0
          -5//6 0 0]
    B₁ = [0 0 0
          1//2 0 0
          -1//2 0 0]
    B₂ = [0 0 0
          1 0 0
          -1 0 0]
    α = [1 // 10; 3 // 14; 24 // 35]

    β₁ = [1; -1; -1]
    β₂ = [0; 1; -1]
    β₃ = [1 // 2; -1 // 4; -1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRI6(T = Float64, T2 = Float64)
    c₀ = [0; 1; 0]
    c₁ = [0; 1; 1]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1 0 0
          0 0 0]
    A₁ = [0 0 0
          1 0 0
          1 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          1 0 0
          0 0 0]
    B₁ = [0 0 0
          1 0 0
          -1 0 0]
    B₂ = [0 0 0
          1 0 0
          -1 0 0]
    α = [1 // 2; 1 // 2; 0]

    β₁ = [1 // 2; 1 // 4; 1 // 4]
    β₂ = [0; 1 // 2; -1 // 2]
    β₃ = [-1 // 2; 1 // 4; 1 // 4]
    β₄ = [0; 1 // 2; -1 // 2]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRDI1WM(T = Float64, T2 = Float64)
    c₀ = [0; 2 // 3]
    c₁ = [0; 0]
    c₂ = [0; 0]
    A₀ = [0 0
          2//3 0]
    A₁ = [0 0
          0 0]
    A₂ = [0 0
          0 0]
    B₀ = [0 0
          2//3 0]
    B₁ = [0 0
          0 0]
    B₂ = [0 0
          0 0]
    α = [1 // 4; 3 // 4]

    β₁ = [1; 0]
    β₂ = [0; 0]
    β₃ = [0; 0]
    β₄ = [0; 0]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRDI2WM(T = Float64, T2 = Float64)
    c₀ = [0; 1; 0]
    c₁ = [0; 2 // 3; 2 // 3]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1 0 0
          0 0 0]
    A₁ = [0 0 0
          2//3 0 0
          2//3 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          1 0 0
          0 0 0]
    B₁ = [0 0 0
          sqrt(2 // 3) 0 0
          -sqrt(2 // 3) 0 0]
    B₂ = [0 0 0
          sqrt(2) 0 0
          -sqrt(2) 0 0]
    α = [1 // 2; 1 // 2; 0]

    β₁ = [1 // 4; 3 // 8; 3 // 8]
    β₂ = [0; sqrt(6) / 4; -sqrt(6) / 4]
    β₃ = [-1 // 4; 1 // 8; 1 // 8]
    β₄ = [0; sqrt(2) / 4; -sqrt(2) / 4]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRDI3WM(T = Float64, T2 = Float64)
    c₀ = [0; 1 // 2; 3 // 4]
    c₁ = [0; 2 // 3; 2 // 3]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1//2 0 0
          0 3//4 0]
    A₁ = [0 0 0
          2//3 0 0
          2//3 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          (9 - 2 * sqrt(15))/14 0 0
          (18 + 3 * sqrt(15))/28 0 0]
    B₁ = [0 0 0
          sqrt(2 // 3) 0 0
          -sqrt(2 // 3) 0 0]
    B₂ = [0 0 0
          sqrt(2) 0 0
          -sqrt(2) 0 0]
    α = [2 // 9; 1 // 3; 4 // 9]

    β₁ = [1 // 4; 3 // 8; 3 // 8]
    β₂ = [0; sqrt(6) / 4; -sqrt(6) / 4]
    β₃ = [-1 // 4; 1 // 8; 1 // 8]
    β₄ = [0; sqrt(2) / 4; -sqrt(2) / 4]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function constructRDI4WM(T = Float64, T2 = Float64)
    c₀ = [0; 1 // 2; 1]
    c₁ = [0; 2 // 3; 2 // 3]
    c₂ = [0; 0; 0]
    A₀ = [0 0 0
          1//2 0 0
          -1 2 0]
    A₁ = [0 0 0
          2//3 0 0
          2//3 0 0]
    A₂ = [0 0 0
          0 0 0
          0 0 0]
    B₀ = [0 0 0
          (6 - sqrt(6))/10 0 0
          (3 + 2 * sqrt(6))/5 0 0]
    B₁ = [0 0 0
          sqrt(2 // 3) 0 0
          -sqrt(2 // 3) 0 0]
    B₂ = [0 0 0
          sqrt(2) 0 0
          -sqrt(2) 0 0]
    α = [1 // 6; 2 // 3; 1 // 6]

    β₁ = [1 // 4; 3 // 8; 3 // 8]
    β₂ = [0; sqrt(6) / 4; -sqrt(6) / 4]
    β₃ = [-1 // 4; 1 // 8; 1 // 8]
    β₄ = [0; sqrt(2) / 4; -sqrt(2) / 4]
    RoesslerRI(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂),
        map(T, α), map(T, β₁), map(T, β₂),
        map(T, β₃), map(T, β₄),
        1 // 1, -0.9674215661017014)
end

function checkRIOrder(RI; tol = 1e-6, ps = 2)
    @unpack c₀, c₁, c₂, A₀, A₁, A₂, B₀, B₁, B₂, α, β₁, β₂, β₃, β₄ = RI
    e = ones(size(α))
    if ps == 2
        conditions = Vector{Bool}(undef, 59) # 9 conditions for first order, 59 in total
    elseif ps == 1
        conditions = Vector{Bool}(undef, 9)
    else
        error("Only the conditions for first and second order weak convergence are implemented. Choose ps=1 or ps=2.")
    end
    #first order weak sense
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₄, e) - 0) < tol
    conditions[3] = abs(dot(β₃, e) - 0) < tol
    conditions[4] = abs(dot(β₁, e)^2 - 1) < tol
    conditions[5] = abs(dot(β₂, e) - 0) < tol
    conditions[6] = abs(dot(β₁, (B₁ * e)) - 0) < tol
    conditions[7] = abs(dot(β₄, (A₂ * e)) - 0) < tol
    conditions[8] = abs(dot(β₃, (B₂ * e)) - 0) < tol
    conditions[9] = abs(dot(β₄, (B₂ * e) .^ 2) - 0) < tol

    #second order weak sense
    if ps == 2
        conditions[10] = abs(dot(α, (A₀ * e)) - 0.5) < tol
        conditions[11] = abs(dot(α, (B₀ * e) .^ 2) - 0.5) < tol
        conditions[12] = abs(dot(β₁, e) * dot(α, (B₀ * e)) - 0.5) < tol
        conditions[13] = abs(dot(β₁, e) * dot(β₁, (A₁ * e)) - 0.5) < tol
        conditions[14] = abs(dot(β₃, (A₂ * e)) - 0) < tol
        conditions[15] = abs(dot(β₂, (B₁ * e)) - 1) < tol
        conditions[16] = abs(dot(β₄, (B₂ * e)) - 1) < tol
        conditions[17] = abs(dot(β₁, e) * dot(β₁, (B₁ * e) .^ 2) - 0.5) < tol
        conditions[18] = abs(dot(β₁, e) * dot(β₃, (B₂ * e) .^ 2) - 0.5) < tol
        conditions[19] = abs(dot(β₁, B₁ * (B₁ * e)) - 0) < tol
        conditions[20] = abs(dot(β₃, B₂ * (B₁ * e)) - 0) < tol
        conditions[21] = abs(dot(β₃, B₂ * (B₁ * (B₁ * e))) - 0) < tol
        conditions[22] = abs(dot(β₁, A₁ * (B₀ * e)) - 0) < tol
        conditions[23] = abs(dot(β₃, A₂ * (B₀ * e)) - 0) < tol
        conditions[24] = abs(dot(β₄, (A₂ * e) .^ 2) - 0) < tol
        conditions[25] = abs(dot(β₄, A₂ * (A₀ * e)) - 0) < tol
        conditions[26] = abs(dot(α, B₀ * (B₁ * e)) - 0) < tol
        conditions[27] = abs(dot(β₂, A₁ * e) - 0) < tol
        conditions[28] = abs(dot(β₁, (A₁ * e) .* (B₁ * e)) - 0) < tol
        conditions[29] = abs(dot(β₃, (A₂ * e) .* (B₂ * e)) - 0) < tol
        conditions[30] = abs(dot(β₄, A₂ * (B₀ * e)) - 0) < tol
        conditions[31] = abs(dot(β₂, A₁ * (B₀ * e)) - 0) < tol
        conditions[32] = abs(dot(β₄, ((B₂ * e) .^ 2) .* (A₂ * e)) - 0) < tol
        conditions[33] = abs(dot(β₄, A₂ * (B₀ * e) .^ 2) - 0) < tol
        conditions[34] = abs(dot(β₂, A₁ * (B₀ * e) .^ 2) - 0) < tol
        conditions[35] = abs(dot(β₁, B₁ * (A₁ * e)) - 0) < tol
        conditions[36] = abs(dot(β₃, B₂ * (A₁ * e)) - 0) < tol
        conditions[37] = abs(dot(β₂, (B₁ * e) .^ 2) - 0) < tol
        conditions[38] = abs(dot(β₄, B₂ * (B₁ * e)) - 0) < tol
        conditions[39] = abs(dot(β₂, B₁ * (B₁ * e)) - 0) < tol
        conditions[40] = abs(dot(β₁, (B₁ * e) .^ 3) - 0) < tol
        conditions[41] = abs(dot(β₃, (B₂ * e) .^ 3) - 0) < tol
        conditions[42] = abs(dot(β₁, B₁ * ((B₁ * e) .^ 2)) - 0) < tol
        conditions[43] = abs(dot(β₃, B₂ * ((B₁ * e) .^ 2)) - 0) < tol
        conditions[44] = abs(dot(β₄, (B₂ * e) .^ 4) - 0) < tol
        conditions[45] = abs(dot(β₄, (B₂ * (B₁ * e)) .^ 2) - 0) < tol
        conditions[46] = abs(dot(β₄, (B₂ * e) .* (B₂ * (B₁ * e))) - 0) < tol
        conditions[47] = abs(dot(α, (B₀ * e) .* (B₀ * (B₁ * e))) - 0) < tol
        conditions[48] = abs(dot(β₁, (A₁ * (B₀ * e)) .* (B₁ * e)) - 0) < tol
        conditions[49] = abs(dot(β₃, (A₂ * (B₀ * e)) .* (B₂ * e)) - 0) < tol
        conditions[50] = abs(dot(β₁, A₁ * (B₀ * (B₁ * e))) - 0) < tol
        conditions[51] = abs(dot(β₃, A₂ * (B₀ * (B₁ * e))) - 0) < tol
        conditions[52] = abs(dot(β₄, (B₂ * (A₁ * e)) .* (B₂ * e)) - 0) < tol
        conditions[53] = abs(dot(β₁, B₁ * (A₁ * (B₀ * e))) - 0) < tol
        conditions[54] = abs(dot(β₃, B₂ * (A₁ * (B₀ * e))) - 0) < tol
        conditions[55] = abs(dot(β₁, (B₁ * e) .* (B₁ * (B₁ * e))) - 0) < tol
        conditions[56] = abs(dot(β₃, (B₂ * e) .* (B₂ * (B₁ * e))) - 0) < tol
        conditions[57] = abs(dot(β₁, B₁ * (B₁ * (B₁ * e))) - 0) < tol
        conditions[58] = abs(dot(β₄, (B₂ * e) .* (B₂ * ((B₁ * e) .^ 2))) - 0) < tol
        conditions[59] = abs(dot(β₄, (B₂ * e) .* (B₂ * (B₁ * (B₁ * e)))) - 0) < tol
    end
    return (conditions)
end

"""

RoesslerRS

Holds the Butcher tableaus for a Roessler RS method. (high weak order in Stratonovich sense)

"""

struct RoesslerRS{T, T2} <: Tableau
    c₀::Vector{T2}
    c₁::Vector{T2}
    c₂::Vector{T2}
    A₀::Matrix{T}
    A₁::Matrix{T}
    A₂::Matrix{T}
    B₀::Matrix{T}
    B₁::Matrix{T}
    B₂::Matrix{T}
    B₃::Matrix{T}
    α::Vector{T}
    β₁::Vector{T}
    β₂::Vector{T}
    order::Rational{Int}
    quantile::T
end

function constructRS1(T = Float64, T2 = Float64)
    c₀ = [0; 0; 1; 0]
    c₁ = [0; 0; 1; 1]
    c₂ = [0; 0; 0; 0]
    A₀ = [0 0 0 0
          0 0 0 0
          1 0 0 0
          0 0 0 0]
    A₁ = [0 0 0 0
          0 0 0 0
          1 0 0 0
          1 0 0 0]
    A₂ = [0 0 0 0
          0 0 0 0
          0 0 0 0
          0 0 0 0]
    B₀ = [0 0 0 0
          0 0 0 0
          1//4 3//4 0 0
          0 0 0 0]
    B₁ = [0 0 0 0
          2//3 0 0 0
          1//12 1//4 0 0
          -5//4 1//4 2 0]
    B₂ = [0 0 0 0
          1 0 0 0
          -1 0 0 0
          0 0 0 0]
    B₃ = [0 0 0 0
          0 0 0 0
          1//4 3//4 0 0
          1//4 3//4 0 0]
    α = [0; 0; 1 // 2; 1 // 2]

    β₁ = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    β₂ = [0; -1 // 4; 1 // 4; 0]

    RoesslerRS(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂), map(T, B₃),
        map(T, α), map(T, β₁), map(T, β₂),
        1 // 1, -0.9674215661017014)
end

function constructRS2(T = Float64, T2 = Float64)
    c₀ = [0; 2 // 3; 2 // 3; 0]
    c₁ = [0; 0; 1; 1]
    c₂ = [0; 0; 0; 0]
    A₀ = [0 0 0 0
          2//3 0 0 0
          1//6 1//2 0 0
          0 0 0 0]
    A₁ = [0 0 0 0
          0 0 0 0
          1 0 0 0
          1 0 0 0]
    A₂ = [0 0 0 0
          0 0 0 0
          0 0 0 0
          0 0 0 0]
    B₀ = [0 0 0 0
          0 0 0 0
          1//4 3//4 0 0
          0 0 0 0]
    B₁ = [0 0 0 0
          2//3 0 0 0
          1//12 1//4 0 0
          -5//4 1//4 2 0]
    B₂ = [0 0 0 0
          1 0 0 0
          -1 0 0 0
          0 0 0 0]
    B₃ = [0 0 0 0
          0 0 0 0
          1//4 3//4 0 0
          1//4 3//4 0 0]
    α = [1 // 4; 1 // 4; 1 // 2; 0]

    β₁ = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    β₂ = [0; -1 // 4; 1 // 4; 0]

    RoesslerRS(map(T2, c₀), map(T2, c₁), map(T2, c₂),
        map(T, A₀), map(T, A₁), map(T, A₂),
        map(T, B₀), map(T, B₁), map(T, B₂), map(T, B₃),
        map(T, α), map(T, β₁), map(T, β₂),
        1 // 1, -0.9674215661017014)
end

function checkRSOrder(RS; tol = 1e-6, ps = 2)
    @unpack c₀, c₁, c₂, A₀, A₁, A₂, B₀, B₁, B₂, B₂, B₃, α, β₁, β₂ = RS
    e = ones(size(α))
    if ps == 2
        conditions = Vector{Bool}(undef, 55) # 6 conditions for first order, 55 in total
    elseif ps == 1
        conditions = Vector{Bool}(undef, 6)
    else
        error("Only the conditions for first and second order weak convergence are implemented. Choose ps=1 or ps=2.")
    end
    #first order weak sense
    conditions[1] = abs(dot(α, e) - 1) < tol
    conditions[2] = abs(dot(β₁, e)^2 - 1) < tol
    conditions[3] = abs(dot(β₂, e) - 0) < tol
    conditions[4] = abs(dot(β₁, (B₁ * e)) - 0.5) < tol
    conditions[5] = abs(dot(β₂, (A₂ * e)) - 0) < tol
    conditions[6] = abs(dot(β₂, (B₂ * e) .^ 2) - 0) < tol

    #second order weak sense
    if ps == 2
        conditions[7] = abs(dot(α, (A₀ * e)) - 0.5) < tol
        conditions[8] = abs(dot(α, B₀ * (B₁ * e)) - 0.25) < tol
        conditions[9] = abs(dot(α, (B₀ * e) .^ 2) - 0.5) < tol
        conditions[10] = abs(dot(β₁, e) * dot(α, (B₀ * e)) - 0.5) < tol
        conditions[11] = abs(dot(β₁, e) * dot(β₁, (A₁ * e)) - 0.5) < tol
        conditions[12] = abs(dot(β₁, B₁ * (A₁ * e)) - 0.25) < tol
        conditions[13] = abs(dot(β₁, (B₁ * e) .* (A₁ * e)) - 0.25) < tol
        conditions[14] = abs(dot(β₁, (B₃ * e)) - 0.5) < tol
        conditions[15] = abs(dot(β₁, e) * dot(β₁, (B₁ * e) .^ 2) - 1 / 3) < tol
        conditions[16] = abs(dot(β₁, e) * dot(β₁, (B₃ * e) .^ 2) - 0.5) < tol
        conditions[17] = abs(dot(β₁, B₃ * (B₃ * e)) - 0) < tol
        conditions[18] = abs(dot(β₂, (B₂ * e))^2 - 0.25) < tol
        conditions[19] = abs(dot(β₁, (B₁ * e) .^ 3) - 0.25) < tol
        conditions[20] = abs(dot(β₁, B₁ * (B₁ * e) .^ 2) - 1 / 12) < tol
        conditions[21] = abs(dot(β₁, B₁ * (B₃ * e) .^ 2) - 0.25) < tol
        conditions[22] = abs(dot(β₁, A₁ * (B₀ * e)) - 0) < tol
        conditions[23] = abs(dot(β₂, (A₂ * e) .^ 2) - 0) < tol
        conditions[24] = abs(dot(β₂, A₂ * (A₀ * e)) - 0) < tol
        conditions[25] = abs(dot(β₁, B₁ * B₁ * (B₁ * e)) - 1 / 24) < tol
        conditions[26] = abs(dot(β₂, A₂ * (B₀ * e)) - 0) < tol
        conditions[27] = abs(dot(β₂, A₂ * (B₀ * e) .^ 2) - 0) < tol
        conditions[28] = abs(dot(β₂, (B₂ * e) .^ 4) - 0) < tol
        conditions[29] = abs(dot(β₂, (B₂ * (B₁ * e)) .^ 2) - 0) < tol
        conditions[30] = abs(dot(β₂, (B₂ * (B₃ * e)) .^ 2) - 0) < tol
        conditions[31] = abs(dot(β₁, (B₁ * e) .* (B₃ * e) .^ 2) - 0.25) < tol
        conditions[32] = abs(dot(β₂, (A₂ * e) .* (B₂ * e) .^ 2) - 0) < tol
        conditions[33] = abs(dot(β₁, B₁ * B₃ * (B₁ * e)) - 1 // 8) < tol
        conditions[34] = abs(dot(β₁, B₃ * B₃ * (B₃ * e)) - 0) < tol
        conditions[35] = abs(dot(β₁, B₃ * B₁ * (B₃ * e)) - 0) < tol
        conditions[36] = abs(dot(β₂, A₂ * B₀ * (B₁ * e)) - 0) < tol
        conditions[37] = abs(dot(β₁, e) * dot(β₁, (B₃ * e) .* (B₁ * e)) - 0.25) < tol
        conditions[38] = abs(dot(β₁, e) * dot(β₁, B₁ * (B₁ * e)) - 1 / 6) < tol
        conditions[39] = abs(dot(β₁, e) * dot(β₁, B₃ * (B₁ * e)) - 0.25) < tol
        conditions[40] = abs(dot(β₁, e) * dot(β₁, B₁ * (B₃ * e)) - 0.25) < tol
        conditions[41] = abs(dot(β₁, (B₁ * e) .* (B₁ * (B₁ * e))) - 1 / 8) < tol
        conditions[42] = abs(dot(β₁, (B₁ * e) .* (B₃ * (B₁ * e))) - 1 / 8) < tol
        conditions[43] = abs(dot(β₁, (B₃ * e) .* (B₁ * (B₃ * e))) - 1 / 4) < tol
        conditions[44] = abs(dot(β₁, (B₃ * e) .* (B₃ * (B₃ * e))) - 0) < tol
        conditions[45] = abs(dot(β₁, B₃ * ((B₃ * e) .* (B₁ * e))) - 0) < tol
        conditions[46] = abs(dot(β₂, (B₂ * (A₁ * e)) .* (B₂ * e)) - 0) < tol
        conditions[47] = abs(dot(β₂, (B₂ * e) .* (B₂ * (B₁ * e))) - 0) < tol
        conditions[48] = abs(dot(β₂, (B₂ * e) .* (B₂ * (B₃ * e))) - 0) < tol
        conditions[49] = abs(dot(β₂, (B₂ * e) .* (B₂ * ((B₁ * e) .^ 2))) - 0) < tol
        conditions[50] = abs(dot(β₂, (B₂ * e) .* (B₂ * ((B₃ * e) .^ 2))) - 0) < tol
        conditions[51] = abs(dot(β₂, (B₂ * e) .* (B₂ * ((B₁ * e) .* (B₃ * e)))) - 0) < tol
        conditions[52] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₁ * (B₁ * e))) - 0) < tol
        conditions[53] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₃ * (B₁ * e))) - 0) < tol
        conditions[54] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₃ * (B₃ * e))) - 0) < tol
        conditions[55] = abs(dot(β₂, (B₂ * e) .* (B₂ * B₁ * (B₃ * e))) - 0) < tol
    end
    return (conditions)
end

"""
NON

Holds the Butcher tableau for Komori's NON method. (second weak order in Stratonovich sense)
"""

struct KomoriNON{T} <: Tableau
    c0::Vector{T}
    cj::Vector{T}
    cjl::Vector{T}
    clj::Vector{T}

    α00::Matrix{T}
    α0j::Matrix{T}
    αj0::Matrix{T}
    αjj::Matrix{T}
    αjl::Matrix{T}

    αjljj::Matrix{T}
    αljjl::Matrix{T}

    order::Rational{Int}
    quantile::T
end

function constructNON(T = Float64)
    c0 = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    cj = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    cjl = [0; 1 // 2; -1 // 2; 0]
    clj = [0; -1 // 2; 1 // 2; 0]

    α00 = [0 0 0 0
           1//2 0 0 0
           0 1//2 0 0
           0 0 1 0]

    α0j = [0 0 0 0
           1 0 0 0
           -9//8 9//8 0 0
           1 0 0 0]

    αj0 = [0 0 0 0
           2 0 0 0
           0 0 0 0
           -2 0 0 0]

    αjj = [0 0 0 0
           2//3 0 0 0
           1//12 1//4 0 0
           -5//4 1//4 2 0]

    αjl = [0 0 0 0
           0 0 0 0
           1//4 3//4 0 0
           1//4 3//4 0 0]

    αjljj = [0 0 0 0
             0 0 0 0
             1 0 0 0
             0 0 0 0]

    αljjl = [0 0 0 0
             -1//2 0 0 0
             1//2 0 0 0
             0 0 0 0]

    KomoriNON(map(T, c0), map(T, cj), map(T, cjl),
        map(T, clj), map(T, α00), map(T, α0j),
        map(T, αj0), map(T, αjj), map(T, αjl), map(T, αjljj), map(T, αljjl),
        1 // 1, -0.9674215661017014)
end

function checkNONOrder(NON; tol = 1e-6)
    if NON isa KomoriNON
        @unpack c0, cj, cjl, clj, α00, α0j, αj0, αjj, αjl, αjljj, αljjl = NON
    elseif NON isa KomoriNON2
        @unpack c0, cj, ckj, α00, α0j, αj0, αjj, αjl, αkjjl = NON
    else
        @unpack c0, cj, α00, α0j, αj0, αjj, αjl = NON
    end
    e = ones(size(c0))

    conditions = Vector{Bool}(undef, 44) # 38 conditions for second order, 6 extra conditions for fourth deterministic order, 44 in total
    conditions[1] = abs(dot(c0, e) - 1) < tol
    conditions[2] = abs(dot(cj, e) - 1) < tol
    conditions[3] = abs(dot(cj, αjj * e) - 1 / 2) < tol
    conditions[4] = abs(dot(cj, α0j * e) - 1 / 2) < tol
    conditions[5] = abs(dot(cj, αj0 * e) - 1 / 2) < tol
    conditions[6] = abs(dot(c0, α00 * e) - 1 / 2) < tol
    conditions[7] = abs(dot(cj, αjj * (αj0 * e)) - 1 / 4) < tol

    conditions[8] = abs(dot(c0, α0j * (αjj * e)) - 1 / 4) < tol
    conditions[9] = abs(dot(cj, αj0 * (α0j * e)) - 0) < tol
    conditions[10] = abs(dot(c0, (α0j * e) .^ 2) - 1 / 2) < tol
    conditions[11] = abs(dot(cj, (αj0 * e) .* (αjj * e)) - 1 / 4) < tol
    conditions[12] = abs(dot(cj, αjj * αjj * (αjj * e)) - 1 / 24) < tol

    conditions[13] = abs(dot(cj, αjj * (αjj * e) .^ 2) - 1 / 12) < tol
    conditions[14] = abs(dot(cj .* (αjj * e), αjj * (αjj * e)) - 1 / 8) < tol
    conditions[15] = abs(dot(cj, (αjj * e) .^ 3) - 1 / 4) < tol
    conditions[16] = abs(dot(cj, αjj * (αjj * e)) - 1 / 6) < tol
    conditions[17] = abs(dot(cj, (αjj * e) .^ 2) - 1 / 3) < tol
    conditions[18] = abs(dot(cj, (αjl * e)) - 1 / 2) < tol
    conditions[19] = abs(dot(cj .* (αjl * e), αjl * (αjl * e)) - 0) < tol
    conditions[20] = abs(dot(cj, (αjl * e) .^ 2) - 1 / 2) < tol
    conditions[21] = abs(dot(cj, αjj * αjl * (αjj * e)) - 1 / 8) < tol
    conditions[22] = abs(dot(cj, αjl * αjl * (αjl * e)) - 0) < tol
    conditions[23] = abs(dot(cj, αjl * αjj * (αjl * e)) - 0) < tol
    conditions[24] = abs(dot(cj, αjj * (αjl * e) .^ 2) - 1 / 4) < tol
    conditions[25] = abs(dot(cj, αjl * ((αjj * e) .* (αjl * e))) - 0) < tol
    conditions[26] = abs(dot(cj .* (αjj * e), αjl * (αjj * e)) - 1 / 8) < tol
    conditions[27] = abs(dot(cj .* (αjl * e), αjj * (αjl * e)) - 1 / 4) < tol
    conditions[28] = abs(dot(cj .* (αjj * e), (αjl * e) .^ 2) - 1 / 4) < tol
    conditions[29] = abs(dot(cj, αjj * (αjl * e)) - 1 / 4) < tol
    conditions[30] = abs(dot(cj, αjl * (αjj * e)) - 1 / 4) < tol
    conditions[31] = abs(dot(cj, αjl * (αjl * e)) - 0) < tol
    conditions[32] = abs(dot(cj, (αjj * e) .* (αjl * e)) - 1 / 4) < tol

    if NON isa KomoriNON
        conditions[33] = abs(dot(clj, e) - 0) < tol
        conditions[34] = abs(dot(cjl, e) - 0) < tol
        conditions[35] = abs(dot(clj, αljjl * e) - 1 / 2) < tol
        conditions[36] = abs(dot(cjl, αjljj * e) + 1 / 2) < tol
        conditions[37] = abs(dot(clj .* (αljjl * e), αljjl * (αjljj * e)) - 0) < tol
        conditions[38] = abs(dot(clj, (αljjl * e) .^ 2) - 0) < tol
    elseif NON isa KomoriNON2
        conditions[33] = abs(ckj[4] + ckj[3] - 0) < tol
        conditions[34] = abs(ckj[2] - 0) < tol
        conditions[35] = abs(αkjjl[4, 3] - 0) < tol
        conditions[36] = true # we set alpha^(k(j),j,0,0) = 0
        conditions[37] = abs(ckj[3] * αkjjl[3, 2]^2 + ckj[4] * αkjjl[4, 2]^2 - 0) < tol
        conditions[38] = abs(ckj[3] * αkjjl[3, 2] + ckj[4] * αkjjl[4, 2] - 1 / 2) < tol
    end
    # deterministic fourth order
    conditions[39] = abs(dot(c0, α00 * α00 * (α00 * e)) - 1 / 24) < tol
    conditions[40] = abs(dot(c0, α00 * (α00 * e) .^ 2) - 1 / 12) < tol
    conditions[41] = abs(dot(c0 .* (α00 * e), α00 * (α00 * e)) - 1 / 8) < tol
    conditions[42] = abs(dot(c0, (α00 * e) .^ 3) - 1 / 4) < tol
    conditions[43] = abs(dot(c0, α00 * (α00 * e)) - 1 / 6) < tol
    conditions[44] = abs(dot(c0, (α00 * e) .^ 2) - 1 / 3) < tol

    return (conditions)
end

struct KomoriNON2{T} <: Tableau
    c0::Vector{T}
    cj::Vector{T}
    ckj::Vector{T}

    α00::Matrix{T}
    α0j::Matrix{T}
    αj0::Matrix{T}
    αjj::Matrix{T}
    αjl::Matrix{T}
    αkjjl::Matrix{T}

    order::Rational{Int}
    quantile::T
end

function constructNON2(T = Float64)
    # gamme is a free parameter
    γ = 1
    c0 = [1 // 6; 1 // 3; 1 // 3; 1 // 6]
    cj = [1 // 8; 3 // 8; 3 // 8; 1 // 8]
    ckj = [0; 0; γ; -γ]

    α00 = [0 0 0 0
           1//2 0 0 0
           0 1//2 0 0
           0 0 1 0]

    α0j = [0 0 0 0
           1 0 0 0
           -9//8 9//8 0 0
           1 0 0 0]

    αj0 = [0 0 0 0
           2 0 0 0
           0 0 0 0
           -2 0 0 0]

    αjj = [0 0 0 0
           2//3 0 0 0
           1//12 1//4 0 0
           -5//4 1//4 2 0]

    αjl = [0 0 0 0
           0 0 0 0
           1//4 3//4 0 0
           1//4 3//4 0 0]

    αkjjl = [0 0 0 0
             0 0 0 0
             0 1/(4 * γ) 0 0
             0 -1/(4 * γ) 0 0]

    KomoriNON2(map(T, c0), map(T, cj), map(T, ckj),
        map(T, α00), map(T, α0j),
        map(T, αj0), map(T, αjj), map(T, αjl), map(T, αkjjl),
        1 // 1, -0.9674215661017014)
end
