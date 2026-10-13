module ArithmeticGeometricMeanAlg
export K, E

include(joinpath(@__DIR__, "complete_reduction.jl"))

# The shared reduction is written against the `_sqrt` seam that `CarlsonAlg` points at an
# unchecked intrinsic. Nothing here needs that, so forward to the generic square root.
_sqrt(x) = √(x)

# Arithmetic Geometric mean implementation
# https://en.wikipedia.org/wiki/Arithmetic%E2%80%93geometric_mean
# https://www.math.emory.edu/~gliang7/AGM.pdf
# Theorem 2.1' and Equation 6
@inline function _ellipke_agm_base(m::T) where {T}
    oneT = one(T)
    twoT = T(2)
    halfπ = T(π / 2)
    ERRTOL = eps(T)

    a = oneT
    b = √(oneT - m)
    sum_c2 = m / twoT
    pow2 = oneT

    while true
        c = (a - b) / twoT
        abs(c) < ERRTOL && break
        next_a = (a + b) / twoT
        next_b = √(a * b)
        sum_c2 = muladd(pow2, c * c, sum_c2)
        pow2 += pow2
        a = next_a
        b = next_b
    end

    k = halfπ / a
    e = k * (oneT - sum_c2)
    return (k, e)
end

"""
`ellipke(m::Real)`
returns `(K(m), E(m))` for scalar `0 ≤ m ≤ 1`
"""
ellipke(m::T) where {T} = _ellipke_reduced(_ellipke_agm_base, m)

"""
``K(m) = \\int_0^{\\pi/2}\\frac{d\\theta}{\\sqrt{1-m\\sin(\\theta)^2}}.``

Returns the complete elliptic integral of the first kind.

# Arguments

- `m` : Elliptic parameter
"""
K(m) = ellipke(m)[1]

"""
``E(m) = \\int_0^{\\pi/2}\\sqrt{1-m\\sin(\\theta)^2}d\\theta.``

Returns the complete elliptic integral of the second kind.

# Arguments

- `m` : Elliptic parameter
"""
E(m) = ellipke(m)[2]
end
