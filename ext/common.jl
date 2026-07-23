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
    sqrt_term = sqrt(1 - m * sin_φ_squared)
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
