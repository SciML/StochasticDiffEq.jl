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

  maxeig!(integrator, cache)
  cache.mdeg = Int(floor(sqrt((2*dt*integrator.eigen_est+1.5)/0.811)+1)) # this is the spectral radius estimate to choose optimal stage
  cache.mdeg = min(3,min(cache.mdeg,200))-2
  choose_deg!(integrator,cache)

  mdeg = cache.mdeg
  η  = cache.optimal_η
  start = cache.start
  deg_index = cache.deg_index
  sqrt_dt = sqrt(dt)
  sqrt_3  = sqrt(3.0)

  for i in 1:length(W.dW)
    winc = rand()*6
    if winc < 1.0
      vec_ξ[i] = -sqrt_dt*sqrt_3
    elseif winc < 2.0
      vec_ξ[i] = sqrt_dt*sqrt_3
    else
      vec_ξ[i] = 0*sqrt_dt
    end

    if !is_diagonal_noise(integrator.sol.prob)
      if rand() < 0.5
        vec_χ[i] = -sqdt
      else
        vec_χ[i] = sqdt
      end
    end
  end

  α = mα[deg_index]
  dt_α = α*dt
  sigma_alpha = (1-α)/2.0 + α*mσ[deg_index]
  tau_alpha   = ((α-1)^2)/2 + 2*α*(1-α)*mσ[deg_index] + (α^2)*mτ[deg_index]

  mu = recf[start]  # here kappa = 0
  t_i = t + dt_α*mu
  t_im1 = t_i
  t_im2 = t
  beta_sm1 = t
  beta_s = t

  # stage 1
  u_im2 = uprev
  k = integrator.f(uprev,p,t)
  u_im1 = uprev + dt_α*mu*k
  (mdeg < 2) && (u_i = u_im1)

  # stages 2 upto s-2
  for i in 2:(mdeg+2)
    mu, kappa = recf[start + 2*(i-2) + 1], recf[start + 2*(i-2) + 2]
    nu = 1.0 + kappa
    u_i = integrator.f(u_im1,p,t)
    u_i = dt_α*mu*u + nu*u_im1 - kappa*u_im2
    u_im2 = u_im1; u_im1 = u_i
    t_i = dt_α*mu + nu*t_im1 - kappa*t_im2
    t_im2 = t_im1
    t_im1 = t_i
  end



  #stage s-1
  mu, kappa = recf2[(deg_index-1)*4 + 1], recf2[(deg_index-1)*4 + 2]
  nu = 1.0 + kappa
  u_i = integrator.f(u_im1,p,t_im1)

  t_star = t_im1 + 2*dt*tau_alpha
  u_star_sm1 = u_im1 + 2*dt*tau_alpha*u_i                   # So that we don't have to calculate f(uₛ₋₂) again
  u = u_im1 + (2*sigma_alpha - 0.5)*dt*u_i

  u_i = dt_α*mu*u_i + nu*u_im1 - kappa*u_im2
  t_im2 = t_im1; t_im1 = t_i; u_im2 = u_im1; uim1 = u_i
  beta_sm1 = t_i

  #stage s
  mu, kappa = recf2[(deg_index-1)*4 + 3], recf2[(deg_index-1)*4 + 4]
  nu = 1.0 + kappa
  u_i = integrator.f(u_im1,p,t_im1)
  u_i = dt_α*mu*u_i + nu*u_im1 - kappa*u_im2
  t_im2 = t_im1; t_im1 = t_i;
  beta_s = t_i

  # Now u_im2 is u_sm2, u_im1 = u_sm1 and u_i = u_s
  # Similarly t_im2 = t_sm2, t_im1 = t_sm1 and t_i = t_s


  g_s = intgrator.g(u_im1,p,t_im1)
  for i in 1:length(W.dW)
    u += @view(g_s[:,i])*vec_ξ[i]
  end

  g_s = integrator.g(u_i,p,t_i)
  for i in 1:length(W.dW)
    u_star_sm1 += @view(g_s[:,i])*vec_ξ[i]
  end

  u_star_sm1 = integrator.f(u_star_sm1,p,t_star)
  u += (1//2)*dt*u_star_sm1

  for i in 1:length(W.dW)
    for j in 1:length(W.dW)
      if i == j
        J_qr = (vec_ξ[i]^2 - 1)/2
      elseif i < j
        J_qr = (vec_ξ[i]*vec_ξ[j] - abs(vec_χ[j])*vec_χ[j])/2
      else
        J_qr = (vec_ξ[i]*vec_ξ[j] + abs(vec_χ[i])*vec_χ[i])/2
      end

      if j == 1
        u_star_sm1 = @view(g_s[:,j])*J_qr
      else
        u_star_sm1 += @view(g_s[:,j])*J_qr
      end
    end

    g_s1 = integrator.g(u_i + u_star_sm1,p,t_i)
    u += (1//2)*@view(g_s1[:,i])
    g_s1 = integrator.g(u_i - u_star_sm1,p,t_i)
    u -= (1//2)*@view(g_s1[:,i])
  end

  for i in 1:length(W.dW)
    if i == 1
      u_star_sm1 = @view(g_s[:,i])*vec_χ[i]
    else
      u_star_sm1 += @view(g_s[:,i])*vec_χ[i]
    end
  end

  g_s1 = integrator.g(u_i + u_star_sm1,p,t_i)
  for i in 1:length(W.dW)
    u += (1//4)*@view(g_s1[:,i])*vec_ξ[i]
  end

  g_s1 = integrator.g(u_i - u_star_sm1,p,t_i)
  for i in 1:length(W.dW)
    u += (1//4)*@view(g_s1[:,i])*vec_ξ[i]
  end

  for i in 1:lenth(W.dW)
    u -= (1//2)*@view(g_s[:,i])*vec_ξ[i]
  end

  integrator.u = u
end

@muladd function perform_step!(integrator,cache::SROCK2Cache,f=integrator.f)
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
