# Shared derivative kernels for the ForwardDiff, Zygote and Enzyme extensions.

@inline function _∂ellipke(K, E, m)
    iszero(m) && return (oftype(K, π) / 8, -oftype(E, π) / 8)
    return ((E - (1 - m) * K) / (2 * (1 - m) * m), (E - K) / (2 * m))
end

@inline function _∂F_E(F, E, ϕ, m)
    sqrt_term = √(1 - m * sin(ϕ)^2)
    ∂ϕF = 1 / sqrt_term
    ∂ϕE = sqrt_term
    sin_2ϕ = sin(2 * ϕ)
    iszero(m) && return (∂ϕF, (2 * ϕ - sin_2ϕ) / 8, ∂ϕE, (sin_2ϕ - 2 * ϕ) / 8)
    ∂mE = (E - F) / (2m)
    if isone(m)
        sec_ϕ = sec(ϕ)
        tan_ϕ = tan(ϕ)
        return (∂ϕF, (sec_ϕ * tan_ϕ - log(abs(sec_ϕ + tan_ϕ))) / 4, ∂ϕE, ∂mE)
    end
    ∂mF = E / (2 * m * (1 - m)) - F / 2 / m - sin_2ϕ / (4 * (1 - m) * sqrt_term)
    return (∂ϕF, ∂mF, ∂ϕE, ∂mE)
end

function _sn_parameter_derivative(E_am, cd, u, m, s, c, d)
    if iszero(m)
        return -c * (u - s * c) / 4
    elseif isone(m)
        return -(s - u * c^2) / 4
    end

    return d * c * ((1 - m) * u - E_am + m * cd * s) / (2 * (1 - m) * m)
end

function _cn_parameter_derivative(E_am, cd, u, m, s, c, d)
    if iszero(m)
        return s * (u - s * c) / 4
    elseif isone(m)
        return s * (s - u * c^2) / (4c)
    end

    return d * s * ((m - 1) * u + E_am - m * cd * s) / (2 * (1 - m) * m)
end

function _dn_parameter_derivative(E_am, cd, u, m, s, c, d)
    iszero(m) && return -s^2 / 2
    ∂m_sn = _sn_parameter_derivative(E_am, cd, u, m, s, c, d)
    return -(s^2 + 2m * s * ∂m_sn) / (2d)
end

function _sc_parameter_derivative(E_am, cd, u, m, s, c, d)
    ∂m_sn = _sn_parameter_derivative(E_am, cd, u, m, s, c, d)
    ∂m_cn = _cn_parameter_derivative(E_am, cd, u, m, s, c, d)
    return (c * ∂m_sn - s * ∂m_cn) / c^2
end

function _complete_j_derivatives(E, K, Pi, J, n, m)
    j = J(n, m)
    k = K(m)
    e = E(m)
    ∂m_k = iszero(m) ? oftype(k, π) / 8 : (e - (1 - m) * k) / (2 * (1 - m) * m)
    ∂m_e = iszero(m) ? -oftype(e, π) / 8 : (e - k) / (2m)

    if iszero(n)
        ∂n_j = iszero(m) ? oftype(j, 3π) / 16 : ((2 + m) * k - 2 * (1 + m) * e) / (3m^2)
        ∂m_j = iszero(m) ? ∂n_j / 2 : (∂m_k - ∂m_e - j) / m
        return j, ∂n_j, ∂m_j
    end

    pi_value = Pi(n, m)
    ∂n_pi = (n * e + (m - n) * k + (n^2 - m) * pi_value) / (2 * (m - n) * (n - 1) * n)
    ∂m_pi = (e / (m - 1) + pi_value) / (2 * (n - m))
    return j, (∂n_pi - j) / n, (∂m_pi - ∂m_k) / n
end

function _incomplete_j_derivatives(E, F, Pi, J, n, φ, m)
    j = J(n, φ, m)
    f = F(φ, m)
    e = E(φ, m)
    sin_φ, cos_φ = sincos(φ)
    sin_φ_squared = sin_φ^2
    sqrt_term = √(1 - m * sin_φ_squared)
    ∂φ_j = sin_φ_squared / ((1 - n * sin_φ_squared) * sqrt_term)
    ∂m_f =
        iszero(m) ? (2φ - sin(2φ)) / 8 :
        e / (2 * m * (1 - m)) - f / (2m) - sin(2φ) / (4 * (1 - m) * sqrt_term)
    ∂m_e = iszero(m) ? (sin(2φ) - 2φ) / 8 : (e - f) / (2m)

    if iszero(n)
        ∂n_j =
            iszero(m) ? 3φ / 8 - sin(2φ) / 4 + sin(4φ) / 32 :
            (m * sin_φ * cos_φ * sqrt_term - m * f + 2 * (1 + m) * (f - e)) / (3m^2)
        ∂m_j = iszero(m) ? ∂n_j / 2 : (∂m_f - ∂m_e - j) / m
        return j, ∂n_j, ∂φ_j, ∂m_j
    end

    pi_value = Pi(n, φ, m)
    ∂n_pi =
        (
            e + (m - n) * f / n + (n^2 - m) * pi_value / n -
            n * sqrt_term * sin(2φ) / (2 * (1 - n * sin_φ_squared))
        ) / (2 * (m - n) * (n - 1))
    ∂m_pi =
        (e / (m - 1) + pi_value - m * sin(2φ) / (2 * (m - 1) * sqrt_term)) / (2 * (n - m))
    return j, (∂n_pi - j) / n, ∂φ_j, (∂m_pi - ∂m_f) / n
end

@inline function _∂Pi(K, E, Pi, n, m)
    iszero(n) && iszero(m) && return (oftype(Pi, π) / 4, oftype(Pi, π) / 8)
    ∂m = (E / (m - 1) + Pi) / (2 * (n - m))
    iszero(n) && return ((K - E) / m, ∂m)
    ∂n = (E + (m - n) * K / n + (n^2 - m) * Pi / n) / (2 * (m - n) * (n - 1))
    return (∂n, ∂m)
end

@inline function _∂Pi(F, E, Pi, n, ϕ, m)
    sin_ϕ = sin(ϕ)
    sin_2ϕ = sin(2 * ϕ)
    sqrt_term = √(1 - m * sin_ϕ^2)
    one_minus_n_sin2 = 1 - n * sin_ϕ^2
    ∂ϕ = 1 / (sqrt_term * one_minus_n_sin2)
    iszero(n) && iszero(m) && return ((2 * ϕ - sin_2ϕ) / 4, ∂ϕ, (2 * ϕ - sin_2ϕ) / 8)
    ∂m = (E / (m - 1) + Pi - m * sin_2ϕ / (2 * (m - 1) * sqrt_term)) / (2 * (n - m))
    iszero(n) && return ((F - E) / m, ∂ϕ, ∂m)
    ∂n =
        (
            E + (m - n) * F / n + (n^2 - m) * Pi / n -
            n * sqrt_term * sin_2ϕ / (2 * one_minus_n_sin2)
        ) / (2 * (m - n) * (n - 1))
    return (∂n, ∂ϕ, ∂m)
end

# sn, cn, dn and E(am(u, m) | m) at one argument. `scd_am` shares a single amplitude
# evaluation between the three Jacobi functions and hands back the amplitude that E
# needs; reaching for the four public entry points re-derives it every time.
@inline function _jacobi_data(scd_am::S, Efun::F, u, m) where {S,F}
    s, c, d, φ = scd_am(u, m)
    return (s, c, d, Efun(φ, m))
end

# ∂/∂u of sn, cn, dn and sc at the same argument, from one (s, c, d) triple.
@inline _∂u_sn(s, c, d, m) = d * c
@inline _∂u_cn(s, c, d, m) = -d * s
@inline _∂u_dn(s, c, d, m) = -m * s * c
@inline _∂u_sc(s, c, d, m) = d / c^2
