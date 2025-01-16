using Test
using SDEProblemLibrary: prob_sde_linear

# Test that the old keyword works, and that the new AliasSpecier works.
@test_warn x -> contains(x, "`alias_u0` keyword argument is deprecated, to set `alias_u0`,
      please use an SDEAliasSpecifier or RODEAliasSpecifier, e.g. `solve(prob, alias = SDEAliasSpecifier(alias_u0 = true))`") solve(prob_sde_linear, EM(), alias_u0 = true, alias_jumps = true)
@test_nowarn solve(prob_sde_linear, EM(), alias = ODEAliasSpecifier(alias_u0 = true, alias_jumps = true))



