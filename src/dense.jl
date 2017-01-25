function sde_interpolant(Θ,integrator)
  (1-Θ).*integrator.uprev .+ Θ.*integrator.u
end

function current_interpolant(t::Number,integrator)
  Θ = (t-integrator.tprev)/integrator.dt
  sde_interpolant(Θ,integrator)
end

function current_interpolant(t::AbstractArray,integrator)
  Θ = (t.-integrator.tprev)./integrator.dt
  [sde_interpolant(ϕ,integrator) for ϕ in Θ]
end
