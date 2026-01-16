@verbosity_specifier SDEVerbosity begin
    toggles = (
        :ode_verbosity,
        :noise_evaluation,
        :adaptive_timestepping,
        :dt_NaN,
        :function_NaN
    )

    presets = (
        None = (
            ode_verbosity = ODEVerbosity(None()),
            noise_evaluation = Silent(),
            adaptive_timestepping = Silent(),
            dt_NaN = Silent(),
            function_NaN = Silent()
        ),
        Minimal = (
            ode_verbosity = ODEVerbosity(Minimal()),
            noise_evaluation = Silent(),
            adaptive_timestepping = Silent(),
            dt_NaN = WarnLevel(),
            function_NaN = WarnLevel()
        ),
        Standard = (
            ode_verbosity = ODEVerbosity(),
            noise_evaluation = Silent(),
            adaptive_timestepping = Silent(),
            dt_NaN = WarnLevel(),
            function_NaN = WarnLevel()
        ),
        Detailed = (
            ode_verbosity = ODEVerbosity(Detailed()),
            noise_evaluation = InfoLevel(),
            adaptive_timestepping = InfoLevel(),
            dt_NaN = WarnLevel(),
            function_NaN = WarnLevel()
        ),
        All = (
            ode_verbosity = ODEVerbosity(All()),
            noise_evaluation = InfoLevel(),
            adaptive_timestepping = InfoLevel(),
            dt_NaN = WarnLevel(),
            function_NaN = WarnLevel()
        )
    )

    groups = (
        sde_specific = (
            :noise_evaluation,
            :adaptive_timestepping,
            :dt_NaN,
            :function_NaN
        ),
    )
end

function Base.getproperty(v::SDEVerbosity, s::Symbol)
    # For parametric types, we need to use the base type
    SDE_fields = fieldnames(typeof(v).name.wrapper)
    if s in SDE_fields
        return getfield(v, s)
    else
        # Try to delegate to ODE verbosity
        ode_v = getfield(v, :ode_verbosity)
        ODE_fields = fieldnames(typeof(ode_v).name.wrapper)
        if s in ODE_fields
            return getfield(ode_v, s)
        else
            return error("type SDEVerbosity has no field ", s)
        end
    end
end


"""
    SDEVerbosity <: AbstractVerbositySpecifier

Verbosity configuration for StochasticDiffEq.jl solvers, providing fine-grained control over
diagnostic messages, warnings, and errors during SDE solution.

# Fields

## ODE Verbosity
- `ode_verbosity`: Verbosity configuration for the underlying ODE solver

## SDE-Specific Group
- `noise_evaluation`: Messages about noise term evaluation
- `adaptive_timestepping`: Messages about adaptive timestepping decisions
- `dt_NaN`: Messages when automatic dt is NaN
- `function_NaN`: Messages when drift or diffusion functions return NaN

# Constructors

    SDEVerbosity(preset::AbstractVerbosityPreset)

Create a `SDEVerbosity` using a preset configuration:
- `SciMLLogging.None()`: All messages disabled
- `SciMLLogging.Minimal()`: Only critical errors and fatal issues
- `SciMLLogging.Standard()`: Balanced verbosity (default)
- `SciMLLogging.Detailed()`: Comprehensive debugging information
- `SciMLLogging.All()`: Maximum verbosity

    SDEVerbosity(; ode_verbosity=nothing, sde_specific=nothing, kwargs...)

Create a `SDEVerbosity` with group-level or individual field control.

# Examples

```julia
# Use a preset
verbose = SDEVerbosity(SciMLLogging.Standard())

# Set ODE verbosity and SDE-specific group
verbose = SDEVerbosity(
    ode_verbosity = ODEVerbosity(SciMLLogging.Detailed()),
    sde_specific = SciMLLogging.InfoLevel()
)

# Set individual fields
verbose = SDEVerbosity(
    noise_evaluation = SciMLLogging.InfoLevel(),
    dt_NaN = SciMLLogging.WarnLevel()
)

# Mix group and individual settings
verbose = SDEVerbosity(
    sde_specific = SciMLLogging.InfoLevel(),  # Set all SDE-specific to InfoLevel
    dt_NaN = SciMLLogging.WarnLevel()  # Override specific field
)
```
"""
function SDEVerbosity end