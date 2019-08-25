using StochasticDiffEq, Test

function f(du, u, p, t)
    du[1] = u[1]
end
function g(du, u, p, t)
    0.0
end

u0 = [1.0]
prob = SDEProblem{true}(f, g, u0, (0.0, 0.1))

sol_ito = solve(prob, RKMil{:Ito}())
@test length(sol_ito) < 100

sol_strato = solve(prob, RKMil{:Stratonovich}(); dt=1e-2)
@test length(sol_strato) < 100

sol_ito = solve(prob, RKMil())
@test length(sol_ito) < 100

sol_strato = solve(prob, RKMil(interpretation=:Stratonovich); dt=1e-2)
@test length(sol_strato) < 100

sol_leh = solve(prob, LambaEulerHeun())
@test length(sol_leh) < 100
