using StochasticDiffEq

const mm_A = [-2.0 1 4
            4 -2 1
            2 1 3]
const mm_b = mm_A*ones(3)
function mm_f(du,u,p,t)
      A_mul_B!(du,mm_A,u)
      tmp = t*mm_b
      du .+= tmp
end
function mm_f(::Type{Val{:analytic}},u0,p,t,W)
      @. 2ones(3)*exp(t) - t - 1
end
function mm_g(du,u,p,t)
      du .= u + t
end
function mm_g(::Type{Val{:analytic}},u0,p,t,W)
      @. 2ones(3)*exp(t) - t - 1
end
function g!(t, u, du)
    du .= 0.0
end
prob2 = SDEProblem(mm_g,g!,ones(3),(0.0,1.0))
prob = SDEProblem(mm_f,g!,ones(3),(0.0,1.0),mass_matrix=mm_A)

sol = solve(prob, ImplicitRKMil(theta=1), dt = 0.01)
sol2 = solve(prob2, ImplicitRKMil(theta=1), dt = 0.01)

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob, ImplicitEM(theta=1), dt = 0.01)
sol2 = solve(prob2, ImplicitEM(theta=1), dt = 0.01)

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob, ImplicitRKMil(symplectic=true), dt = 0.01)
sol2 = solve(prob2, ImplicitRKMil(symplectic=true), dt = 0.01)

@test norm(sol .- sol2) ≈ 0 atol=1e-11

sol = solve(prob, ImplicitEM(symplectic=true), dt = 0.01)
sol2 = solve(prob2, ImplicitEM(symplectic=true), dt = 0.01)

@test norm(sol .- sol2) ≈ 0 atol=1e-11

function mm_f2(du,u,p,t)
      A_mul_B!(du,mm_A,u)
end
function no_mm_f2(du,u,p,t)
      du .= u
end
function no_mm_g2(du,u,p,t)
      du .= u
end
function mm_g2(t, u, du)
    A_mul_B!(du,mm_A,u)
end
prob2 = SDEProblem(no_mm_f2,no_mm_g2,ones(3),(0.0,1.0))
prob = SDEProblem(mm_f2,no_mm_g2,ones(3),(0.0,1.0),mass_matrix=mm_A)

srand(1)
sol = solve(prob, ImplicitEM(theta=1), dt = 0.01)
srand(1)
sol2 = solve(prob2, ImplicitEM(theta=1), dt = 0.01)

@test norm(sol .- sol2) ≈ 0 atol=1e-11
