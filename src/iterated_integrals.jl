# This function calculates WikJ, a mxm Array for a m dimensional general noise problem,
# which is a approximation to the second order iterated integrals
# this is implementation of the section 4 of the paper doi:10.1016/j.cam.2006.05.037

function get_iterated_I!(integrator, cache::StochasticDiffEqConstantCache)
    @unpack dt, u, uprev, t, p, W = integrator
    @unpack m_seq, WikJ = cache

    dW     = W.dW
    m      = length(dW)
    M      = m*(m-1)/2

    if typef(dW) <: Number || is_diagonal_noise(integrator.sol.prob)
        cache.WikJ = 1//2 .* J .^ 2
    else
        sum_dW² = zero(eltype(dW))
        for i in 1:length(dW)
            sum_dW² += dW[i]^2
        end

        Gp1 = randn(M)
        α = sqrt(1 + sum_dW²/dt)
        Gp2 = Gp1/(sqrt(2)*(1+α)*dt)

        #operator (Iₘ² - Pₘ)Kₘᵀ
        for i in 1:M
            WikJ2[m_seq[i,1], m_seq[i,2]] = Gp2[i]
            WikJ2[m_seq[i,2], m_seq[i,1]] = -Gp2[i]
        end

        #operator (Iₘ X W*Wᵀ)
        WikJ2 = dW*(dW'*WikJ2)

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

        WikJ = 1//2*(dW*dW')
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
    end
    cache.WikJ = WikJ
    return false
end
