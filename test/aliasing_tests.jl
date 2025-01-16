using Test
using StochasticDiffEq
using SDEProblemLibrary: prob_sde_linear
using SciMLBase

# Test that the old keyword works, and that the new AliasSpecier works.
@test_warn x -> contains("`alias_u0`", x) solve(prob_sde_linear, EM(), dt=0.1, alias_u0=true)

@test_nowarn solve(prob_sde_linear, EM(), dt = 0.1, alias = SciMLBase.SDEAliasSpecifier(alias_u0 = true))



