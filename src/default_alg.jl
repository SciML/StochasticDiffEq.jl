"""
Default algorithm selection for StochasticDiffEq.jl
Based on the logic from DifferentialEquations.jl but using a type-stable approach
similar to OrdinaryDiffEq.jl's DefaultODEAlgorithm
"""

EnumX.@enumx DefaultSDESolverChoice begin
    SOSRI = 1
    SOSRA = 2
    RKMilCommute = 3
    ImplicitRKMil = 4
    SKenCarp = 5
    LambaEM = 6
    ISSEM = 7
    LambaEulerHeun = 8
    ImplicitEulerHeun = 9
end

const NUM_NONSTIFF_SDE = 3  # SOSRI, SOSRA, RKMilCommute
const NUM_STIFF_SDE = 2     # ImplicitRKMil, SKenCarp

"""
    DefaultSDEAlgorithm(; kwargs...)

Constructs a default SDE algorithm that automatically switches between different
solvers based on problem characteristics. This provides a type-stable default
algorithm similar to OrdinaryDiffEq.jl's DefaultODEAlgorithm.

The algorithm selection is based on:

  - Stiffness detection
  - Noise characteristics (additive, commutative, general)
  - Problem interpretation (Ito vs Stratonovich)

## Keyword Arguments

  - `lazy::Bool=true`: Whether to use lazy interpolants
  - `stiffalgfirst::Bool=false`: Whether to start with the stiff algorithm
  - `autodiff::Union{Bool,ADTypes.AbstractADType}=true`: Automatic differentiation backend
  - `kwargs...`: Additional keyword arguments passed to the algorithms
"""
function DefaultSDEAlgorithm(; lazy = true, stiffalgfirst = false,
        autodiff = true, kwargs...)
    # For now, we use a simpler approach with just two algorithms
    # This can be expanded later to include more sophisticated selection
    nonstiff = SOSRI()
    stiff = ImplicitRKMil(autodiff = autodiff)

    AutoAlgSwitch(nonstiff, stiff; stiffalgfirst = stiffalgfirst)
end

"""
    DefaultAdaptiveSDEAlgorithm(; kwargs...)

Constructs a more sophisticated default SDE algorithm that adapts based on
problem characteristics including noise type and algorithm hints.

## Keyword Arguments

  - `alg_hints::Vector{Symbol}=Symbol[]`: Algorithm hints like :additive, :commutative, :stiff, :stratonovich
  - `autodiff::Union{Bool,ADTypes.AbstractADType}=true`: Automatic differentiation backend
  - `kwargs...`: Additional keyword arguments passed to the algorithms
"""
function DefaultAdaptiveSDEAlgorithm(; alg_hints = Symbol[],
        autodiff = true, kwargs...)
    is_stiff = :stiff ∈ alg_hints
    is_stratonovich = :stratonovich ∈ alg_hints
    is_additive = :additive ∈ alg_hints
    is_commutative = :commutative ∈ alg_hints

    # Select appropriate algorithms based on hints
    if is_additive
        nonstiff = SOSRA()
        stiff = ISSEM(autodiff = autodiff)  # Use ISSEM for stiff additive noise
    elseif is_commutative
        nonstiff = RKMilCommute()
        stiff = ImplicitRKMil(autodiff = autodiff)
    else
        nonstiff = SOSRI()
        if is_stratonovich
            stiff = ImplicitRKMil(autodiff = autodiff,
                interpretation = SciMLBase.AlgorithmInterpretation.Stratonovich)
        else
            stiff = ImplicitRKMil(autodiff = autodiff)
        end
    end

    AutoAlgSwitch(nonstiff, stiff; stiffalgfirst = is_stiff)
end

"""
    default_sde_alg_choice(prob::SDEProblem)

Determines the default algorithm choice for an SDE problem based on its characteristics.
Returns algorithm hints that can be used to select appropriate solvers.
"""
function default_sde_alg_choice(prob::SDEProblem)
    alg_hints = Symbol[]

    # Check for mass matrix (implies stiffness)
    if hasproperty(prob.f, :mass_matrix) && prob.f.mass_matrix !== I
        push!(alg_hints, :stiff)
    end

    # Check for additive noise (diagonal noise with no state dependence)
    if prob.noise === nothing || isa(prob.noise, DiffEqNoiseProcess.NoiseProcess)
        if DiffEqBase.is_diagonal_noise(prob)
            # Simple heuristic for additive noise
            push!(alg_hints, :additive)
        end
    end

    # Check problem interpretation
    if hasproperty(prob, :interpretation) &&
       prob.interpretation == SciMLBase.AlgorithmInterpretation.Stratonovich
        push!(alg_hints, :stratonovich)
    end

    return alg_hints
end

# Hook into the solve interface
function SciMLBase.__init(prob::SDEProblem, ::Nothing, args...; kwargs...)
    alg_hints = default_sde_alg_choice(prob)
    alg = if isempty(alg_hints)
        DefaultSDEAlgorithm(; kwargs...)
    else
        DefaultAdaptiveSDEAlgorithm(; alg_hints = alg_hints, kwargs...)
    end
    SciMLBase.__init(prob, alg, args...; kwargs...)
end

function SciMLBase.__solve(prob::SDEProblem, ::Nothing, args...; kwargs...)
    alg_hints = default_sde_alg_choice(prob)
    alg = if isempty(alg_hints)
        DefaultSDEAlgorithm(; kwargs...)
    else
        DefaultAdaptiveSDEAlgorithm(; alg_hints = alg_hints, kwargs...)
    end
    SciMLBase.__solve(prob, alg, args...; kwargs...)
end

# Mark default algorithms for special handling
function isdefaultalg(alg::StochasticCompositeAlgorithm)
    # Check if it's one of our default algorithm combinations
    if isa(alg.algs, Tuple{SOSRI, ImplicitRKMil}) ||
       isa(alg.algs, Tuple{SOSRA, ISSEM}) ||
       isa(alg.algs, Tuple{RKMilCommute, ImplicitRKMil})
        return true
    end
    return false
end

# Also provide a fallback for other algorithms
isdefaultalg(alg) = false

export DefaultSDEAlgorithm, DefaultAdaptiveSDEAlgorithm, default_sde_alg_choice,
       isdefaultalg
