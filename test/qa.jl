using StochasticDiffEq, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(StochasticDiffEq)
    Aqua.test_ambiguities(StochasticDiffEq, recursive = false)
    Aqua.test_deps_compat(StochasticDiffEq)
    Aqua.test_piracies(StochasticDiffEq,
        treat_as_own = [StochasticDiffEq.JumpProcesses.JumpProblem,
            StochasticDiffEq.SciMLBase.AbstractRODEProblem])
    Aqua.test_project_extras(StochasticDiffEq)
    Aqua.test_stale_deps(StochasticDiffEq)
    Aqua.test_unbound_args(StochasticDiffEq)
    Aqua.test_undefined_exports(StochasticDiffEq)
end
