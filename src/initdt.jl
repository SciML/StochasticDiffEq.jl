function sde_determine_initdt{tType,uType}(u0::uType,t::tType,tdir,dtmax,abstol,reltol,internalnorm,prob,order)
  f = prob.f
  g = prob.g
  d₀ = internalnorm(u0./(abstol+abs.(u0)*reltol))
  if !isinplace(prob)
    f₀ = f(t,u0)
    g₀ = 3g(t,u0)
  else
    f₀ = similar(u0)
    g₀ = similar(u0)
    f(t,u0,f₀)
    g(t,u0,g₀); g₀.*=3
  end

  d₁ = internalnorm(max.(abs.(f₀.+g₀),abs.(f₀.-g₀))./(abstol.+abs.(u0).*reltol))
  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1e-6)
  else
    dt₀ = tType(0.01*(d₀/d₁))
  end
  dt₀ = min(dt₀,tdir*dtmax)
  u₁ = u0 .+ tdir.*dt₀.*f₀
  if !isinplace(prob)
    f₁ = f(t+tdir*dt₀,u₁)
    g₁ = 3g(t+tdir*dt₀,u₁)
  else
    f₁ = similar(u0)
    g₁ = similar(u0)
    f(t,u0,f₁)
    g(t,u0,g₁); g₁.*=3
  end
  ΔgMax = max.(abs.(g₀.-g₁),abs.(g₀.+g₁))
  d₂ = internalnorm(max.(abs.(f₁.-f₀.+ΔgMax),abs.(f₁.-f₀.-ΔgMax))./(abstol.+abs.(u0).*reltol))./dt₀
  if max(d₁,d₂)<=T1(1e-15)
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)))/(order+.5)))
  end
  dt = tdir*min(100dt₀,dt₁)
end
