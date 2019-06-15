function sde_determine_initdt(u0::uType,t::tType,tdir,dtmax,abstol,reltol,internalnorm,prob,order,integrator) where {tType,uType}
  f = prob.f
  g = prob.g
  p = prob.p
  d₀ = internalnorm(u0./(abstol.+internalnorm.(u0,t).*reltol),t)
  if !isinplace(prob)
    f₀ = f(u0,p,t)
    if integrator.opts.verbose && any(x->any(isnan,x),f₀)
      @warn("First function call for f produced NaNs. Exiting.")
    end
    g₀ = 3g(u0,p,t)
    if integrator.opts.verbose && any(x->any(isnan,x),g₀)
      @warn("First function call for g produced NaNs. Exiting.")
    end
  else
    f₀ = zero(u0)
    if prob.noise_rate_prototype !== nothing
      g₀ = zero(prob.noise_rate_prototype)
    else
      g₀ = zero(u0)
    end
    f(f₀,u0,p,t)
    if integrator.opts.verbose && any(x->any(isnan,x),f₀)
      @warn("First function call for f produced NaNs. Exiting.")
    end
    g(g₀,u0,p,t); g₀.*=3
    if integrator.opts.verbose && any(x->any(isnan,x),g₀)
      @warn("First function call for g produced NaNs. Exiting.")
    end
  end

  d₁ = internalnorm(max.(internalnorm.(f₀.+g₀,t),internalnorm.(f₀.-g₀,t))./(abstol.+internalnorm.(u0,t).*reltol),t)
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
    f₁ = f(u₁,p,t+tdir*dt₀)
    g₁ = 3g(u₁,p,t+tdir*dt₀)
  else
    f₁ = zero(u0)
    if prob.noise_rate_prototype !== nothing
      g₁ = zero(prob.noise_rate_prototype)
    else
      g₁ = zero(u0)
    end
    f(f₁,u0,p,t)
    g(g₁,u0,p,t); g₁.*=3
  end
  ΔgMax = max.(internalnorm.(g₀.-g₁,t),internalnorm.(g₀.+g₁,t))
  d₂ = internalnorm(max.(internalnorm.(f₁.-f₀.+ΔgMax,t),internalnorm.(f₁.-f₀.-ΔgMax,t))./(abstol.+internalnorm.(u0,t).*reltol),t)./dt₀
  if max(d₁,d₂)<=T1(1//Int64(10)^(15))
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)))/(order+.5)))
  end
  dt = tdir*min(100dt₀,dt₁,tdir*dtmax)
end
