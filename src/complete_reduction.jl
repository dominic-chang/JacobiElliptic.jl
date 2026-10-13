# Argument reduction shared by every complete-integral implementation.
#
# The AGM and Carlson kernels are only valid on `0 ≤ m < 1`. Two standard transformations
# cover the rest of the real line, and because each one lands directly inside `[0, 1)`, a
# single step always suffices:
#
# - `m < 0` : imaginary-modulus transformation, A&S 17.4.15 and 17.4.16
# - `m > 1` : reciprocal-modulus transformation, A&S 17.4.17
#
# `base` evaluates the kernel on the reduced parameter and returns the pair `(K, E)`.
# Callers that only want `K` pass a base whose second component is a placeholder; the
# compiler drops the resulting dead arithmetic.
@inline function _ellipke_reduced(base::F, m::T) where {F,T}
    oneT = one(T)
    if isnan(m)
        return (T(NaN), T(NaN))
    elseif m < zero(T)
        m_abs = -m
        one_plus_m_abs = m_abs + oneT
        k, e = _ellipke_base_guarded(base, m_abs / one_plus_m_abs)
        sqrt_factor = _sqrt(one_plus_m_abs)
        return (k / sqrt_factor, e * sqrt_factor)
    elseif m > oneT
        modulus = _sqrt(m)
        m_inv = inv(m)
        k, e = _ellipke_base_guarded(base, m_inv)
        return (inv(modulus) * k, modulus * (e - (oneT - m_inv) * k))
    elseif m == oneT
        return (T(Inf), oneT)
    else
        return base(m)
    end
end

# Reduction does not always land strictly inside `[0, 1)`: `-m / (1 - m)` rounds to exactly
# 1 once `m ≤ -9.1e15`, and `m == -Inf` produces `Inf / Inf == NaN`. Both would send the
# kernel off to a degenerate parameter, so they get the same answers the unreduced cases do.
@inline function _ellipke_base_guarded(base::F, m::T) where {F,T}
    if isnan(m)
        return (T(NaN), T(NaN))
    elseif m == one(T)
        return (T(Inf), one(T))
    else
        return base(m)
    end
end
