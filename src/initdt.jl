function sde_determine_initdt(u0,t,tdir,dtmax,abstol,reltol,internalnorm,prob,order)
  f = prob.f
  g = prob.g
  d₀ = internalnorm(u0./(abstol+abs.(u0)*reltol))
  if typeof(u0) <: Number
    f₀ = f(t,u0)
    g₀ = 3g(t,u0)
  else
    f₀ = similar(u0)
    g₀ = similar(u0)
    f(t,u0,f₀)
    g(t,u0,g₀); g₀.*=3
  end

  d₁ = internalnorm(max(abs.(f₀.+g₀),abs.(f₀-g₀))./(abstol+abs.(u0)*reltol))
  if d₀ < 1e-5 || d₁ < 1e-5
    dt₀ = 1e-6
  else
    dt₀ = 0.01*(d₀/d₁)
  end
  dt₀ = min(dt₀,tdir*dtmax)
  u₁ = u0 + tdir*dt₀*f₀
  if typeof(u0) <: Number
    f₁ = f(t+tdir*dt₀,u₁)
    g₁ = 3g(t+tdir*dt₀,u₁)
  else
    f₁ = similar(u0)
    g₁ = similar(u0)
    f(t,u0,f₁)
    g(t,u0,g₁); g₁.*=3
  end
  ΔgMax = max(abs.(g₀-g₁),abs.(g₀+g₁))
  d₂ = internalnorm(max(abs.(f₁.-f₀.+ΔgMax),abs.(f₁.-f₀.-ΔgMax))./(abstol+abs.(u0)*reltol))/dt₀
  if max(d₁,d₂)<=1e-15
    dt₁ = max(1e-6,dt₀*1e-3)
  else
    dt₁ = 10.0^(-(2+log10(max(d₁,d₂)))/(order+.5))
  end
  dt = tdir*min(100dt₀,dt₁)
end
