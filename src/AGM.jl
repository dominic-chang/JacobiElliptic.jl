module ArithmeticGeometricMeanAlg
export K, E

# Arithmetic Geometric mean implementation
# https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
# https://www.math.emory.edu/~gliang7/AGM.pdf
# Theorem 2.1' and Equation 6
@inline function _ellipke_agm_base(m::T) where {T}
    oneT = one(T)
    twoT = T(2)
    halfpi = T(π / 2)
    ERRTOL = eps(T)

    a = oneT
    b = sqrt(oneT - m)
    sum_c2 = m / twoT
    pow2 = oneT

    while true
        c = (a - b) / twoT
        abs(c) < ERRTOL && break
        next_a = (a + b) / twoT
        next_b = sqrt(a * b)
        sum_c2 = muladd(pow2, c * c, sum_c2)
        pow2 += pow2
        a = next_a
        b = next_b
    end

    k = halfpi / a
    e = k * (oneT - sum_c2)
    return (k, e)
end

"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
function ellipke(m::T) where {T}
    oneT = one(T)
    if isnan(m)
        return (T(NaN), T(NaN))
    elseif m < zero(T)
        # Transformation for m < 0
        m_abs = -m
        one_plus_m_abs = m_abs + oneT
        m_transformed = m_abs / one_plus_m_abs
        k_trans, e_trans = ellipke(m_transformed)
        sqrt_factor = sqrt(one_plus_m_abs)
        k_result = k_trans / sqrt_factor
        e_result = e_trans * sqrt_factor
        return (k_result, e_result)
    elseif m > oneT
        # Reciprocal modulus transformation for m > 1
        k = sqrt(m)
        k_inv = inv(k)
        m_inv = inv(m)
        k_complete, e_complete = ellipke(m_inv)
        k_transformed = k_inv * k_complete
        e_transformed = k * (e_complete - (1 - m_inv) * k_complete)
        return (k_transformed, e_transformed)
    elseif m == oneT
        return (T(Inf), oneT)
    else
        # 0 ≤ m < 1
        return _ellipke_agm_base(m)
    end
end
"""
``K(m) = \\int_0^{\\pi/2}\\frac{d\\theta}{\\sqrt{1-k^2\\sin(\\theta)^2}}.``

Returns the complete elliptic integral of the first kind.
    
# Arguments

- `m` : Elliptic modulus
"""
K(m) = ellipke(m)[1]

"""
``E(m) = \\int_0^{\\pi/2}\\sqrt{1-k^2\\sin(\\theta)^2}d\\theta.``

Returns the complete elliptic integral of the second kind.

# Arguments

- `m` : Elliptic modulus
"""
E(m) = ellipke(m)[2]
end
