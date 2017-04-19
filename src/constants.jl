@inline ODE_DEFAULT_ISOUTOFDOMAIN(t,u) = false
@inline UNITLESS_ABS2(x) = abs2(x)/(typeof(x)(one(x))*typeof(x)(one(x)))
@inline ODE_DEFAULT_NORM(u) = sqrt(sum(UNITLESS_ABS2,u) / length(u))
@inline ODE_DEFAULT_PROG_MESSAGE(dt,t,u) = "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,t,u) = any(isnan,u)
