# This function calculates WikJ, a mxm Array for a m dimensional general noise problem,
# which is a approximation to the second order iterated integrals
# this is implementation of the section 4 of the paper doi:10.1016/j.cam.2006.05.037

function get_iterated_I!(integrator, cache::StochasticDiffEqConstantCache)
    @unpack dt, u, uprev, t, p, W = integrator
    @unpack WikJ = cache
    dW     = W.dW

    if typef(dW) <: Number || is_diagonal_noise(integrator.sol.prob)
        cache.WikJ = 1//2 .* dW .^ 2
    else
        m      = length(dW)
        M      = m*(m-1)/2
        m_seq  = cache.m_seq
        # sum_dW² = zero(eltype(dW))
        sum_dW² = dW'*dW

        WikJ = dW*dW'
        Gp1 = randn(M)
        α = sqrt(1 + sum_dW²/dt)
        Gp2 = Gp1/(sqrt(2)*(1+α)*dt)

        #operator (Iₘ² - Pₘ)Kₘᵀ
        for i in 1:M
            WikJ2[m_seq[i,1], m_seq[i,2]] = Gp2[i]
            WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp2[i]
        end

        #operator (Iₘ X W*Wᵀ)
        WikJ2 = WikJ*WikJ2

        #operator Kₘ(Iₘ² - Pₘ)
        WikJ2 = WikJ2 - WikJ2'
        for i in 1:M
            Gp2[i] = WikJ2[m_seq[i,1], m_seq[i,2]]
        end
        Gp = Gp/sqrt(2) + Gp2

        #operator (Iₘ² - Pₘ)Kₘᵀ
        for i in 1:M
            WikJ2[m_seq[i,1], m_seq[i,2]] = Gp[i]
            WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp[i]
        end

        WikJ *= 1//2
        a2ₚ = (π^2)/6
        p = Int(floor((1/π)*sqrt(M/(24*dt))*sqrt(m + 4*sum_dW²/dt) + 1))
        for i in 1:p
            a2ₚ -= (1/i^2)
            var = sqrt(dt/(2*π*i))
            vec_ζ = randn(m)*var
            vec_η = randn(m)*var
            WikJ += (vec_ζ*vec_η' - vec_η*vec_ζ')
            Aₚ -= (2/sqrt(π*i))*vec_ζ
        end

        WikJ -= 1//2*(dW*Aₚ' - Aₚ*dW')
        WikJ += (sqrt(a2ₚ)*dt/π)*WikJ2

        cache.WikJ = WikJ
    end
    return false
end

function get_iterated_I!(integrator, cache::StochasticDiffEqMutableCache)
    @unpack dt, u, uprev, t, p, W = integrator
    @unpack WikJ = cache
    dW = W.dW

    if typef(dW) <: Number || is_diagonal_noise(integrator.sol.prob)
        @.. cache.WikJ = 1//2 .* dW .^ 2
    else
        m      = length(dW)
        M      = m*(m-1)/2
        m_seq  = cache.m_seq; WikJ2 = cache.WikJ2; WikJ3 = cache.WikJ3;
        Gp1    = cache.Gp1; Gp2 = cache.Gp2
        vec_ζ  = cache.vec_ζ; vec_η = cache.vec_η

        sum_dW² = zero(eltype(dW))
        mul!(sum_dW²,dW', dW)

        @.. Gp1 = randn(M)
        α = sqrt(1 + sum_dW²/dt)
        @.. Gp2 = Gp1/(sqrt(2)*(1+α)*dt)

        #operator (Iₘ² - Pₘ)Kₘᵀ
        for i in 1:M
            WikJ2[m_seq[i,1], m_seq[i,2]] = Gp2[i]
            WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp2[i]
        end

        #operator (Iₘ X W*Wᵀ)
        mul!(WikJ,dW,dW')
        mul!(WikJ3,WikJ,WikJ2)

        #operator Kₘ(Iₘ² - Pₘ)
        @.. WikJ2 = WikJ3 - WikJ3'
        for i in 1:M
            Gp2[i] = WikJ2[m_seq[i,1], m_seq[i,2]]
        end
        Gp1 = Gp1/sqrt(2) + Gp2

        #operator (Iₘ² - Pₘ)Kₘᵀ
        for i in 1:M
            WikJ2[m_seq[i,1], m_seq[i,2]] = Gp1[i]
            WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp1[i]
        end

        @.. WikJ *= 1//2
        a2ₚ = (π^2)/6
        p = Int(floor((1/π)*sqrt(M/(24*dt))*sqrt(m + 4*sum_dW²/dt) + 1))
        for i in 1:p
            a2ₚ -= (1/i^2)
            var = sqrt(dt/(2*π*i))
            @.. vec_ζ = randn(m)*var
            @.. vec_η = randn(m)*var
            mul!(WikJ3, vec_ζ, vec_η')
            @.. WikJ += WikJ3 - WikJ3'
            @.. Aₚ -= (2/sqrt(π*i))*vec_ζ
        end
        mul!(WikJ3, dW, Aₚ')
        @.. WikJ -= 1//2*(WikJ3 - WikJ3')
        @.. WikJ += (sqrt(a2ₚ)*dt/π)*WikJ2

        @.. cache.WikJ = WikJ
    end
    return false
end
