@muladd function perform_step!(integrator,cache::SROCK1ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator

  maxeig!(integrator, cache)
  cache.mdeg = Int(floor(sqrt(2*dt*integrator.eigen_est)+1)) # this is the spectral radius estimate to choose optimal stage
  choose_deg!(integrator,cache)

  mdeg = cache.mdeg
  η  = cache.optimal_η
  ω₀ = 1.0 + (η/(mdeg^2))
  ωSq = (ω₀^2) - 1.0
  Sqrt_ω = sqrt(ωSq)
  cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
  ω₁ = (Sqrt_ω*cosh(mdeg*cosh_inv))/(mdeg*sinh(mdeg*cosh_inv))

  if alg_interpretation(integrator.alg) == :Stratonovich
    α  = cosh(mdeg*cosh_inv)/(2*ω₀*cosh((mdeg-1)*cosh_inv))
    γ  = 1/(2*α)
    β  = -γ
  end

  uᵢ₋₂ = copy(uprev)
  k = integrator.f(uprev,p,t)
  Tᵢ₋₂ = one(eltype(u))
  Tᵢ₋₁ = convert(eltype(u),ω₀)
  Tᵢ   = Tᵢ₋₁
  tᵢ₋₁ = t+dt*(ω₁/ω₀)
  tᵢ   = tᵢ₋₁
  tᵢ₋₂ = t
  gₘ₋₁ = zero(k)
  gₘ₋₂ = zero(k)

  #stage 1
  uᵢ₋₁ = uprev + (dt*ω₁/ω₀)*k

  for i in 2:mdeg
    Tᵢ = 2*ω₀*Tᵢ₋₁ - Tᵢ₋₂
    μ = 2*ω₁*(Tᵢ₋₁/Tᵢ)
    ν = 2*ω₀*(Tᵢ₋₁/Tᵢ)
    κ = (-Tᵢ₋₂/Tᵢ)
    k = integrator.f(uᵢ₋₁,p,tᵢ₋₁)

    u = dt*μ*k + ν*uᵢ₋₁ + κ*uᵢ₋₂
    if (i > mdeg - 2) && alg_interpretation(integrator.alg) == :Stratonovich
        (i == mdeg - 1) && (gₘ₋₂ = integrator.g(uᵢ₋₁,p,tᵢ₋₁); u += α*gₘ₋₂.*W.dW)
        (i == mdeg) && (gₘ₋₁ = integrator.g(uᵢ₋₁,p,tᵢ₋₁); u += (β*gₘ₋₂ + γ*gₘ₋₁) .* W.dW)
    elseif (i == mdeg) && alg_interpretation(integrator.alg) == :Ito
      if typeof(W.dW) <: Number || is_diagonal_noise(integrator.sol.prob)
        gₘ₋₂ = integrator.g(uᵢ₋₁,p,tᵢ₋₁)
        uᵢ₋₂ = uᵢ₋₁ + sqrt(dt)*gₘ₋₂
        gₘ₋₁ = integrator.g(uᵢ₋₂,p,tᵢ₋₁)
        u += gₘ₋₂ .* W.dW + 1/(2.0*sqrt(dt)) .* (gₘ₋₁ - gₘ₋₂) .* (W.dW^2 - dt)
      else
        gₘ₋₁ = integrator.g(uᵢ₋₁,p,tᵢ₋₁)
        u += gₘ₋₁ .* W.dW
      end
    end

    if i < mdeg
      tᵢ = μ*dt + ν*tᵢ₋₁ + κ*tᵢ₋₂
      uᵢ₋₂ = uᵢ₋₁
      uᵢ₋₁ = u
      tᵢ₋₂ = tᵢ₋₁
      tᵢ₋₁ = tᵢ
      Tᵢ₋₂ = Tᵢ₋₁
      Tᵢ₋₁ = Tᵢ
    end
  end
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SROCK1Cache,f=integrator.f)
  @unpack uᵢ₋₁,uᵢ₋₂,k, gₘ₋₁, gₘ₋₂ = cache
  @unpack t,dt,uprev,u,W,p = integrator
  ccache = cache.constantcache
  maxeig!(integrator, cache)
  ccache.mdeg = Int(floor(sqrt(2*dt*integrator.eigen_est)+1))   # this is the spectral radius estimate to choose optimal stage
  choose_deg!(integrator,cache)

  mdeg = ccache.mdeg
  η  = ccache.optimal_η
  ω₀ = 1 + η/(mdeg^2)
  ωSq = ω₀^2 - 1
  Sqrt_ω = sqrt(ωSq)
  cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
  ω₁ = (Sqrt_ω*cosh(mdeg*cosh_inv))/(mdeg*sinh(mdeg*cosh_inv))

  if alg_interpretation(integrator.alg) == :Stratonovich
    α  = cosh(mdeg*cosh_inv)/(2*ω₀*cosh((mdeg-1)*cosh_inv))
    γ  = 1/(2*α)
    β  = -γ
  end

  @.. uᵢ₋₂ = uprev
  integrator.f(k,uprev,p,t)
  Tᵢ₋₂ = one(eltype(u))
  Tᵢ₋₁ = convert(eltype(u),ω₀)
  Tᵢ   = Tᵢ₋₁
  tᵢ₋₁ = t + dt*(ω₁/ω₀)
  tᵢ   = tᵢ₋₁
  tᵢ₋₂ = t

  #stage 1
  @.. uᵢ₋₁ = uprev + (dt*ω₁/ω₀)*k

  for i in 2:mdeg
    Tᵢ = 2*ω₀*Tᵢ₋₁ - Tᵢ₋₂
    μ = 2*ω₁*Tᵢ₋₁/Tᵢ
    ν = 2*ω₀*Tᵢ₋₁/Tᵢ
    κ = - Tᵢ₋₂/Tᵢ
    integrator.f(k,uᵢ₋₁,p,tᵢ₋₁)
    @.. u = dt*μ*k + ν*uᵢ₋₁ + κ*uᵢ₋₂
    if (i > mdeg - 2) && alg_interpretation(integrator.alg) == :Stratonovich
        (i == mdeg - 1) && (integrator.g(gₘ₋₂,uᵢ₋₁,p,tᵢ₋₁); @.. u += α*gₘ₋₂*W.dW)
        (i == mdeg) && (integrator.g(gₘ₋₁,uᵢ₋₁,p,tᵢ₋₁); @.. u += (β*gₘ₋₂ + γ*gₘ₋₁)*W.dW)
    elseif (i == mdeg) && alg_interpretation(integrator.alg) == :Ito
      if typeof(W.dW) <: Number || is_diagonal_noise(integrator.sol.prob)
        integrator.g(gₘ₋₂,uᵢ₋₁,p,tᵢ₋₁)
        @.. uᵢ₋₂ = uᵢ₋₁ + sqrt(dt)*gₘ₋₂
        integrator.g(gₘ₋₁,uᵢ₋₂,p,tᵢ₋₁)
        @.. u += gₘ₋₂*W.dW + 1/(2.0*sqrt(dt))*(gₘ₋₁ - gₘ₋₂)*(W.dW^2 - dt)
      else
        integrator.g(gₘ₋₁,uᵢ₋₁,p,tᵢ₋₁)
        @.. u += gₘ₋₁*W.dW
      end
    end

    if i < mdeg
      tᵢ = dt*μ + ν*tᵢ₋₁ + κ*tᵢ₋₂
      @.. uᵢ₋₂ = uᵢ₋₁
      @.. uᵢ₋₁ = u
      tᵢ₋₂ = tᵢ₋₁
      tᵢ₋₁ = tᵢ
      Tᵢ₋₂ = Tᵢ₋₁
      Tᵢ₋₁ = Tᵢ
    end
  end
  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SROCK2ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator
  @unpack recf, recf2, mα, mσ, mτ = cache

  gen_prob = !((is_diagonal_noise(integrator.sol.prob)) || (typeof(W.dW) <: Number) || (length(W.dW) == 1))
  gen_prob && (vec_χ = zeros(eltype(W.dW),length(W.dW)))

  maxeig!(integrator, cache)
  cache.mdeg = Int(floor(sqrt((2*dt*integrator.eigen_est+1.5)/0.811)+1))
  cache.mdeg = max(3,min(cache.mdeg,200))-2
  choose_deg!(integrator,cache)

  mdeg      = cache.mdeg
  start     = cache.start
  deg_index = cache.deg_index
  α = mα[deg_index]
  σ = (1.0-α)*0.5 + α*mσ[deg_index]

  # I'm not sure about which one is correct τ
  τ = 0.5*((1.0-α)^2) + 2*α*(1.0-α)*mσ[deg_index] + (α^2.0)*(mσ[deg_index]*(mσ[deg_index]+mτ[deg_index]))
  # τ = 0.5*((1.0-α)^2) + 2*α*(1.0-α)*mσ[deg_index] + (α^2.0)*mτ[deg_index]

  sqrt_dt   = sqrt(dt)
  sqrt_3    = sqrt(3*one(eltype(W.dW)))

  if gen_prob
    for i in 1:length(W.dW)
      if rand() < 0.5
        vec_χ[i] = -one(eltype(W.dW))
      else
        vec_χ[i] = one(eltype(W.dW))
      end
    end
  end

  μ = recf[start]  # here κ = 0
  tᵢ = t + α*dt*μ
  tᵢ₋₁ = tᵢ
  tᵢ₋₂ = t

  # stage 1
  uᵢ₋₂ = uprev
  uᵢ = integrator.f(uprev,p,t)
  uᵢ₋₁ = uprev + α*dt*μ*uᵢ

  # stages 2 upto s-2
  for i in 2:mdeg
    μ, κ = recf[start + 2*(i-2) + 1], recf[start + 2*(i-2) + 2]
    ν    = 1.0 + κ
    uᵢ   = integrator.f(uᵢ₋₁,p,tᵢ₋₁)
    uᵢ   = α*dt*μ*uᵢ + ν*uᵢ₋₁ - κ*uᵢ₋₂
    uᵢ₋₂ = uᵢ₋₁
    uᵢ₋₁ = uᵢ
    tᵢ   = α*dt*μ + ν*tᵢ₋₁ - κ*tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
    tᵢ₋₁ = tᵢ
  end

  #2 stage-finishing procedure
  #stage s-1
  μ, κ = recf2[(deg_index-1)*4 + 1], recf2[(deg_index-1)*4 + 2]
  ν    = 1.0 + κ
  uᵢ   = integrator.f(uᵢ₋₁,p,tᵢ₋₁)

  tₓ   = tᵢ₋₁ + 2*dt*τ                    # So that we don't have to calculate f(uₛ₋₂) again
  uₓ   = uᵢ₋₁ + 2*dt*τ*uᵢ                 # uₓ and tₓ represent u_star
  u    = uᵢ₋₁ + (2*σ - 0.5)*dt*uᵢ

  uᵢ   = α*dt*μ*uᵢ + ν*uᵢ₋₁ - κ*uᵢ₋₂
  tᵢ   = α*dt*μ + ν*tᵢ₋₁ - κ*tᵢ₋₂
  tᵢ₋₂ = tᵢ₋₁
  tᵢ₋₁ = tᵢ
  uᵢ₋₂ = uᵢ₋₁
  uᵢ₋₁ = uᵢ

  #stage s
  μ, κ = recf2[(deg_index-1)*4 + 3], recf2[(deg_index-1)*4 + 4]
  ν    = 1.0 + κ
  uᵢ   = integrator.f(uᵢ₋₁,p,tᵢ₋₁)
  uᵢ   = α*dt*μ*uᵢ + ν*uᵢ₋₁ - κ*uᵢ₋₂
  tᵢ   = α*dt*μ + ν*tᵢ₋₁ - κ*tᵢ₋₂

  # Now uᵢ₋₂ = uₛ₋₂, uᵢ₋₁ = uₛ₋₁, uᵢ = uₛ
  # Similarly tᵢ₋₂ = tₛ₋₂, tᵢ₋₁ = tₛ₋₁, tᵢ = tₛ

  if (typeof(W.dW) <: Number) || (length(W.dW) == 1)
    Gₛ = integrator.g(uᵢ₋₁,p,tᵢ₋₁)
    u  += Gₛ*W.dW
    Gₛ = integrator.g(uᵢ,p,tᵢ)
    uₓ += Gₛ*W.dW

    uₓ = integrator.f(uₓ,p,tₓ)
    u  += (1//2*dt)*uₓ
    uₓ  = Gₛ*((W.dW^2-dt)/2)
    uᵢ₋₂ = uᵢ + uₓ
    Gₛ₁ = integrator.g(uᵢ₋₂,p,tᵢ)
    u   += (1//2)*Gₛ₁
    uᵢ₋₂ = uᵢ - uₓ
    Gₛ₁ = integrator.g(uᵢ₋₂,p,tᵢ)
    u   -= (1//2)*Gₛ₁

    uₓ  = Gₛ*sqrt_dt
    uᵢ₋₂ = uᵢ + uₓ
    Gₛ₁ = integrator.g(uᵢ₋₂,p,tᵢ)
    u += (1//4*W.dW)*(Gₛ₁ - 2*Gₛ)
    uᵢ₋₂ = uᵢ - uₓ
    Gₛ₁ = integrator.g(uᵢ₋₂,p,tᵢ)
    u  += (1//4*W.dW)*Gₛ₁
  elseif is_diagonal_noise(integrator.sol.prob)

    Gₛ = integrator.g(uᵢ₋₁,p,tᵢ₋₁)
    u += Gₛ .* W.dW
    Gₛ = integrator.g(uᵢ,p,tᵢ)
    uₓ += Gₛ .* W.dW

    uₓ = integrator.f(uₓ,p,tₓ)
    u  += (1//2)*dt*uₓ

    uₓ   = Gₛ .* ((W.dW .^ 2 .- dt)/2)
    uᵢ₋₂ = uᵢ + uₓ
    Gₛ₁  = integrator.g(uᵢ₋₂,p,tᵢ)
    u    += (1//2)*Gₛ₁
    uᵢ₋₂ = uᵢ - uₓ
    Gₛ₁  = integrator.g(uᵢ₋₂,p,tᵢ)
    u    -= (1//2)*Gₛ₁

    uₓ = sqrt_dt*Gₛ
    uᵢ₋₂ = uᵢ + uₓ
    Gₛ₁  = integrator.g(uᵢ₋₂,p,tᵢ)
    u += (1//4*W.dW) .* (Gₛ₁ .- 2*Gₛ)

    uᵢ₋₂ = uᵢ - uₓ
    Gₛ₁  = integrator.g(uᵢ₋₂,p,tᵢ)
    u += (1//4*W.dW) .* Gₛ₁
  else
      Gₛ = integrator.g(uᵢ₋₁,p,tᵢ₋₁)
      for i in 1:length(W.dW)
        u += @view(Gₛ[:,i])*(W.dW[i])
      end
      Gₛ = integrator.g(uᵢ,p,tᵢ)
      for i in 1:length(W.dW)
        uₓ += @view(Gₛ[:,i])*(W.dW[i])
      end

      uₓ = integrator.f(uₓ,p,tₓ)
      u  += (1//2)*dt*uₓ
      for i in 1:length(W.dW)
        for j in 1:length(W.dW)
          (i == j) && (Jᵢⱼ = (W.dW[i]^2 - dt)/2)
          (i < j) && (Jᵢⱼ = (W.dW[i]*W.dW[j])/2)
          (i > j) && (Jᵢⱼ = (W.dW[i]*W.dW[j])/2)

          (j == 1) && (uₓ = @view(Gₛ[:,j])*Jᵢⱼ)
          (j > 1) && (uₓ += @view(Gₛ[:,j])*Jᵢⱼ)
        end
        Gₛ₁ = integrator.g(uᵢ + uₓ,p,tᵢ)
        u   += (1//2)*@view(Gₛ₁[:,i])
        Gₛ₁ = integrator.g(uᵢ - uₓ,p,tᵢ)
        u   -= (1//2)*@view(Gₛ₁[:,i])
      end

      for i in 1:length(W.dW)
        (i == 1) && (uₓ = @view(Gₛ[:,i])*(vec_χ[i]*sqrt_dt))
        (i > 1) && (uₓ += @view(Gₛ[:,i])*(vec_χ[i]*sqrt_dt))
      end
      uᵢ₋₂ = uᵢ + uₓ
      Gₛ₁ = integrator.g(uᵢ₋₂,p,tᵢ)
      for i in 1:length(W.dW)
        u += (1//4*W.dW[i])*(@view(Gₛ₁[:,i]) - 2*@view(Gₛ[:,i]))
      end
      uᵢ₋₂ = uᵢ - uₓ
      Gₛ₁ = integrator.g(uᵢ₋₂,p,tᵢ)
      for i in 1:length(W.dW)
        u += (1//4*W.dW[i])*@view(Gₛ₁[:,i])
      end
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SROCK2Cache,f=integrator.f)
  @unpack uᵢ, uₓ, uᵢ₋₁, uᵢ₋₂, k, Gₛ, Gₛ₁, vec_χ = cache

  @unpack t,dt,uprev,u,W,p = integrator

  @unpack recf, recf2, mα, mσ, mτ = cache.constantcache
  ccache = cache.constantcache
  gen_prob = !((is_diagonal_noise(integrator.sol.prob)) || (typeof(W.dW) <: Number) || (length(W.dW) == 1))

  maxeig!(integrator, cache)
  ccache.mdeg = Int(floor(sqrt((2*dt*integrator.eigen_est+1.5)/0.811)+1))
  ccache.mdeg = max(3,min(ccache.mdeg,200))-2
  choose_deg!(integrator,cache)

  mdeg      = ccache.mdeg
  start     = ccache.start
  deg_index = ccache.deg_index
  α = mα[deg_index]
  σ = (1.0-α)*0.5 + α*mσ[deg_index]

  # I'm not sure about which one is correct τ
  τ = 0.5*((1.0-α)^2) + 2*α*(1.0-α)*mσ[deg_index] + (α^2.0)*(mσ[deg_index]*(mσ[deg_index]+mτ[deg_index]))
  # τ = 0.5*((1.0-α)^2) + 2*α*(1.0-α)*mσ[deg_index] + (α^2.0)*mτ[deg_index]

  sqrt_dt   = sqrt(dt)
  sqrt_3    = sqrt(3.0)
  if gen_prob
    for i in 1:length(W.dW)
      if rand() < 0.5
        vec_χ[i] = -sqrt_dt
      else
        vec_χ[i] = sqrt_dt
      end
    end
  end

  μ = recf[start]  # here κ = 0
  tᵢ = t + α*dt*μ
  tᵢ₋₁ = tᵢ
  tᵢ₋₂ = t

  # stage 1
  @.. uᵢ₋₂ = uprev
  integrator.f(k,uprev,p,t)
  @.. uᵢ₋₁ = uprev + α*dt*μ*k

  # stages 2 upto s-2
  for i in 2:mdeg
    μ, κ = recf[start + 2*(i-2) + 1], recf[start + 2*(i-2) + 2]
    ν    = 1.0 + κ
    integrator.f(k,uᵢ₋₁,p,t)
    @.. uᵢ   = α*dt*μ*k + ν*uᵢ₋₁ - κ*uᵢ₋₂
    @.. uᵢ₋₂ = uᵢ₋₁
    @.. uᵢ₋₁ = uᵢ
    tᵢ   = α*dt*μ + ν*tᵢ₋₁ - κ*tᵢ₋₂
    tᵢ₋₂ = tᵢ₋₁
    tᵢ₋₁ = tᵢ
  end

  #2 stage-finishing procedure
  #stage s-1
  μ, κ = recf2[(deg_index-1)*4 + 1], recf2[(deg_index-1)*4 + 2]
  ν    = 1.0 + κ
  integrator.f(k,uᵢ₋₁,p,tᵢ₋₁)

  tₓ   = tᵢ₋₁ + 2*dt*τ                    # So that we don't have to calculate f(uₛ₋₂) again
  @.. uₓ   = uᵢ₋₁ + 2*dt*τ*k                 # uₓ and tₓ represent u_star
  @.. u    = uᵢ₋₁ + (2*σ - 0.5)*dt*k

  @.. uᵢ   = α*dt*μ*k + ν*uᵢ₋₁ - κ*uᵢ₋₂
  tᵢ   = α*dt*μ + ν*tᵢ₋₁ - κ*tᵢ₋₂
  tᵢ₋₂ = tᵢ₋₁
  tᵢ₋₁ = tᵢ
  @.. uᵢ₋₂ = uᵢ₋₁
  @.. uᵢ₋₁ = uᵢ

  #stage s
  μ, κ = recf2[(deg_index-1)*4 + 3], recf2[(deg_index-1)*4 + 4]
  ν    = 1.0 + κ
  integrator.f(k,uᵢ₋₁,p,tᵢ₋₁)
  @.. uᵢ   = α*dt*μ*k + ν*uᵢ₋₁ - κ*uᵢ₋₂
  tᵢ   = α*dt*μ + ν*tᵢ₋₁ - κ*tᵢ₋₂

  # Now uᵢ₋₂ = uₛ₋₂, uᵢ₋₁ = uₛ₋₁, uᵢ = uₛ
  # Similarly tᵢ₋₂ = tₛ₋₂, tᵢ₋₁ = tₛ₋₁, tᵢ = tₛ

  if (typeof(W.dW) <: Number) || (length(W.dW) == 1)
    integrator.g(Gₛ,uᵢ₋₁,p,tᵢ₋₁)
    @.. u  += Gₛ*W.dW
    integrator.g(Gₛ,uᵢ,p,tᵢ)
    @.. uₓ += Gₛ*W.dW
    integrator.f(k,uₓ,p,tₓ)
    @.. u  += (1//2)*dt*k

    @.. uₓ  = Gₛ*((W.dW^2 - dt)/2)
    @.. uᵢ₋₂ = uᵢ + uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u   += (1//2)*Gₛ₁
    @.. uᵢ₋₂ = uᵢ - uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u   -= (1//2)*Gₛ₁

    @.. uₓ = sqrt_dt*Gₛ
    @.. uᵢ₋₂ = uᵢ + uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u += (1//4*W.dW)*(Gₛ₁ - 2*Gₛ)
    @.. uᵢ₋₂ = uᵢ - uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u += (1//4*W.dW)*Gₛ₁
  elseif is_diagonal_noise(integrator.sol.prob)
    integrator.g(Gₛ,uᵢ₋₁,p,tᵢ₋₁)
    @.. u += Gₛ*W.dW

    integrator.g(Gₛ,uᵢ,p,tᵢ)
    @.. uₓ += Gₛ*W.dW

    integrator.f(k,uₓ,p,tₓ)
    @.. u  += (1//2)*dt*k

    @.. uₓ   = Gₛ*((W.dW^2 - dt)/2)
    @.. uᵢ₋₂ = uᵢ + uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u += (1//2)*Gₛ₁
    @.. uᵢ₋₂ = uᵢ - uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u -= (1//2)*Gₛ₁

    @.. uₓ = Gₛ
    @.. uₓ   *= sqrt_dt
    @.. uᵢ₋₂ = uᵢ + uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u += (1//4)*W.dW*(Gₛ₁-2*Gₛ)
    @.. uᵢ₋₂ = uᵢ - uₓ
    integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
    @.. u += (1//4)*W.dW*Gₛ₁
  else
      integrator.g(Gₛ,uᵢ₋₁,p,tᵢ₋₁)
      for i in 1:length(W.dW)
        @.. u += @view(Gₛ[:,i])*W.dW[i]
      end
      integrator.g(Gₛ,uᵢ,p,tᵢ)
      for i in 1:length(W.dW)
        @.. uₓ += @view(Gₛ[:,i])*W.dW[i]
      end
      integrator.f(k,uₓ,p,tₓ)
      @.. u  += (1//2)*dt*k

      for i in 1:length(W.dW)
        for j in 1:length(W.dW)
          (i == j) && (Jᵢⱼ = (W.dW[i]^2 - dt)/2)
          (i < j) && (Jᵢⱼ = (W.dW[i]*W.dW[j] - dt*vec_χ[i])/2)
          (i > j) && (Jᵢⱼ = (W.dW[i]*W.dW[j] + dt*vec_χ[i])/2)

          (j == 1) && (@.. uₓ = @view(Gₛ[:,j])*Jᵢⱼ)
          (j > 1) && (@.. uₓ += @view(Gₛ[:,j])*Jᵢⱼ)
        end
        @.. uᵢ₋₂ = uᵢ + uₓ
        integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
        @.. u   += (1//2)*@view(Gₛ₁[:,i])
        @.. uᵢ₋₂ = uᵢ - uₓ
        integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
        @.. u   -= (1//2)*@view(Gₛ₁[:,i])
      end

      for i in 1:length(W.dW)
        (i == 1) && (@.. uₓ = @view(Gₛ[:,i])*(vec_χ[i]*sqrt_dt))
        (i > 1) && (@.. uₓ += @view(Gₛ[:,i])*(vec_χ[i]*sqrt_dt))
      end
      @.. uᵢ₋₂ = uᵢ + uₓ
      integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
      for i in 1:length(W.dW)
        @.. u += (1//4)*W.dW[i]*(@view(Gₛ₁[:,i]) - 2*@view(Gₛ[:,i]))
      end
      @.. uᵢ₋₂ = uᵢ - uₓ
      integrator.g(Gₛ₁,uᵢ₋₂,p,tᵢ)
      for i in 1:length(W.dW)
        @.. u += (1//4)*W.dW[i]*@view(Gₛ₁[:,i])
      end
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SROCKEMConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u,W,p = integrator

  maxeig!(integrator, cache)
  if integrator.alg.strong_order_1
    cache.mdeg = Int(floor(sqrt(2*dt*integrator.eigen_est/0.19)+1))
  else
    cache.mdeg = Int(floor(sqrt(2*dt*integrator.eigen_est/0.33)+1))
  end
  cache.mdeg = max(3,min(cache.mdeg,200))
  choose_deg!(integrator,cache)

  mdeg = cache.mdeg
  η  = cache.optimal_η
  ω₀ = 1.0 + (η/(mdeg^2))
  ωSq = (ω₀^2) - 1.0
  Sqrt_ω = sqrt(ωSq)
  cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
  ω₁ = (Sqrt_ω*cosh(mdeg*cosh_inv))/(mdeg*sinh(mdeg*cosh_inv))

  uᵢ₋₂ = copy(uprev)
  k = integrator.f(uprev,p,t)
  Tᵢ₋₂ = one(eltype(u))
  Tᵢ₋₁ = convert(eltype(u),ω₀)
  Tᵢ   = Tᵢ₋₁
  tᵢ₋₁ = t+dt*(ω₁/ω₀)
  tᵢ   = tᵢ₋₁
  tᵢ₋₂ = t

  #stage 1
  uᵢ₋₁ = uprev + (dt*ω₁/ω₀)*k

  for i in 2:mdeg
    Tᵢ = 2*ω₀*Tᵢ₋₁ - Tᵢ₋₂
    μ = 2*ω₁*(Tᵢ₋₁/Tᵢ)
    ν = 2*ω₀*(Tᵢ₋₁/Tᵢ)
    κ = (-Tᵢ₋₂/Tᵢ)
    k = integrator.f(uᵢ₋₁,p,tᵢ₋₁)

    u = dt*μ*k + ν*uᵢ₋₁ + κ*uᵢ₋₂
    tᵢ = μ*dt + ν*tᵢ₋₁ + κ*tᵢ₋₂

    if i < mdeg
      uᵢ₋₂ = uᵢ₋₁
      uᵢ₋₁ = u
      tᵢ₋₂ = tᵢ₋₁
      tᵢ₋₁ = tᵢ
      Tᵢ₋₂ = Tᵢ₋₁
      Tᵢ₋₁ = Tᵢ
    end
  end

  Gₛ = integrator.g(u,p,tᵢ)
  if (typeof(W.dW) <: Number) || (length(W.dW) == 1)
    uᵢ₋₁ = Gₛ*W.dW
  elseif is_diagonal_noise(integrator.sol.prob)
    uᵢ₋₁ = Gₛ .* W.dW
  else
    for i in 1:length(W.dW)
      (i == 1) && (uᵢ₋₁ = @view(Gₛ[:,i])*W.dW[i])
      (i > 1) && (uᵢ₋₁ += @view(Gₛ[:,i])*W.dW[i])
    end
  end

  if integrator.alg.strong_order_1
    if (typeof(W.dW) <: Number) || (length(W.dW) == 1) || (is_diagonal_noise(integrator.sol.prob))
      if (typeof(W.dW) <: Number) || (length(W.dW) == 1)
        uᵢ₋₂  = Gₛ*(W.dW^2 - dt)*0.5
      else
        uᵢ₋₂  = Gₛ .* (W.dW .^ 2 .- dt) .* 0.5
      end

      tmp = u + uᵢ₋₂
      Gₛ  = integrator.g(tmp,p,tᵢ)
      uᵢ₋₁ += 0.5*Gₛ
      tmp = u - uᵢ₋₂
      Gₛ  = integrator.g(tmp,p,tᵢ)
      uᵢ₋₁ -= 0.5*Gₛ
    else
      for i in 1:length(W.dW)
        for j in 1:length(W.dW)
          (i == j) && (Jᵢⱼ = (W.dW[i]*W.dW[j]-dt)*0.5)
          (i != j) && (Jᵢⱼ = (W.dW[i]*W.dW[j])*0.5)

          (j == 1) && (uᵢ₋₂ = @view(Gₛ[:,j])*Jᵢⱼ)
          (j > 1) && (uᵢ₋₂ += @view(Gₛ[:,j])*Jᵢⱼ)
        end
        tmp = u + uᵢ₋₂
        Gₛ₁ = integrator.g(tmp,p,tᵢ)
        uᵢ₋₁ += @view(Gₛ₁[:,i])*(0.5)
        tmp = u - uᵢ₋₂
        Gₛ₁ = integrator.g(tmp,p,tᵢ)
        uᵢ₋₁ -= @view(Gₛ₁[:,i])*(0.5)
      end
    end
  end

  integrator.u = u + uᵢ₋₁
end

@muladd function perform_step!(integrator,cache::SROCKEMCache,f=integrator.f)
  @unpack uᵢ₋₁,uᵢ₋₂,tmp,k,Gₛ,Gₛ₁ = cache
  @unpack t,dt,uprev,u,W,p = integrator
  ccache = cache.constantcache
  maxeig!(integrator, cache)
  if integrator.alg.strong_order_1
    ccache.mdeg = Int(floor(sqrt(2*dt*integrator.eigen_est/0.19)+1))
  else
    ccache.mdeg = Int(floor(sqrt(2*dt*integrator.eigen_est/0.33)+1))
  end
  ccache.mdeg = max(3,min(ccache.mdeg,200))
  choose_deg!(integrator,cache)

  mdeg = ccache.mdeg
  η  = ccache.optimal_η
  ω₀ = 1 + η/(mdeg^2)
  ωSq = ω₀^2 - 1
  Sqrt_ω = sqrt(ωSq)
  cosh_inv = log(ω₀ + Sqrt_ω)             # arcosh(ω₀)
  ω₁ = (Sqrt_ω*cosh(mdeg*cosh_inv))/(mdeg*sinh(mdeg*cosh_inv))

  @.. uᵢ₋₂ = uprev
  integrator.f(k,uprev,p,t)
  Tᵢ₋₂ = one(eltype(u))
  Tᵢ₋₁ = convert(eltype(u),ω₀)
  Tᵢ   = Tᵢ₋₁
  tᵢ₋₁ = t + dt*(ω₁/ω₀)
  tᵢ   = tᵢ₋₁
  tᵢ₋₂ = t

  #stage 1
  @.. uᵢ₋₁ = uprev + (dt*ω₁/ω₀)*k

  for i in 2:mdeg
    Tᵢ = 2*ω₀*Tᵢ₋₁ - Tᵢ₋₂
    μ = 2*ω₁*Tᵢ₋₁/Tᵢ
    ν = 2*ω₀*Tᵢ₋₁/Tᵢ
    κ = - Tᵢ₋₂/Tᵢ
    integrator.f(k,uᵢ₋₁,p,tᵢ₋₁)
    @.. u = dt*μ*k + ν*uᵢ₋₁ + κ*uᵢ₋₂
    tᵢ = dt*μ + ν*tᵢ₋₁ + κ*tᵢ₋₂

    if i < mdeg
      @.. uᵢ₋₂ = uᵢ₋₁
      @.. uᵢ₋₁ = u
      tᵢ₋₂ = tᵢ₋₁
      tᵢ₋₁ = tᵢ
      Tᵢ₋₂ = Tᵢ₋₁
      Tᵢ₋₁ = Tᵢ
    end
  end

  integrator.g(Gₛ,u,p,tᵢ)
  if (typeof(W.dW) <: Number) || (length(W.dW) == 1) || is_diagonal_noise(integrator.sol.prob)
    @.. uᵢ₋₁ = Gₛ*W.dW
  else
    for i in 1:length(W.dW)
      (i == 1) && (@.. uᵢ₋₁ = @view(Gₛ[:,i])*W.dW[i])
      (i > 1) && (@.. uᵢ₋₁ += @view(Gₛ[:,i])*W.dW[i])
    end
  end

  if integrator.alg.strong_order_1
    if (typeof(W.dW) <: Number) || (length(W.dW) == 1) || (is_diagonal_noise(integrator.sol.prob))
      @.. uᵢ₋₂  = Gₛ*(W.dW^2 - dt)*0.5
      @.. tmp = u + uᵢ₋₂
      integrator.g(Gₛ,tmp,p,tᵢ)
      @.. uᵢ₋₁ += 0.5*Gₛ
      @.. tmp = u - uᵢ₋₂
      integrator.g(Gₛ,tmp,p,tᵢ)
      @.. uᵢ₋₁ -= 0.5*Gₛ
    else
      for i in 1:length(W.dW)
        for j in 1:length(W.dW)
          (i == j) && (Jᵢⱼ = (W.dW[i]*W.dW[j]-dt)*0.5)
          (i != j) && (Jᵢⱼ = (W.dW[i]*W.dW[j])*0.5)

          (j == 1) && (@.. uᵢ₋₂ = @view(Gₛ[:,j])*Jᵢⱼ)
          (j > 1) && (@.. uᵢ₋₂ += @view(Gₛ[:,j])*Jᵢⱼ)
        end
        @.. tmp = u + uᵢ₋₂
        integrator.g(Gₛ₁,tmp,p,tᵢ)
        @.. uᵢ₋₁ += @view(Gₛ₁[:,i])*(0.5)
        @.. tmp = u - uᵢ₋₂
        integrator.g(Gₛ₁,tmp,p,tᵢ)
        @.. uᵢ₋₁ -= @view(Gₛ₁[:,i])*(0.5)
      end
    end
  end

  integrator.u = u + uᵢ₋₁
end
