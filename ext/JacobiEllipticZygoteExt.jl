module JacobiEllipticZygoteExt

using JacobiElliptic
using Zygote: @adjoint

#----------------------------------------------------------------------------------------
# Elliptic E(a, b)
#----------------------------------------------------------------------------------------

@adjoint JacobiElliptic.CarlsonAlg.E(a, b) = JacobiElliptic.CarlsonAlg.E(a, b),
(
    c̄ ->
        c̄ .* (
            (sqrt(1 - b * sin(a)^2)),
            iszero(b) ? -π / 8 :
            (JacobiElliptic.CarlsonAlg.E(a, b) - JacobiElliptic.CarlsonAlg.F(a, b)) / (2b),
        )
)

#----------------------------------------------------------------------------------------
# Elliptic CN(a, b)
#----------------------------------------------------------------------------------------

@adjoint JacobiElliptic.cn(a, b) = JacobiElliptic.cn(a, b),
(
    c̄ ->
        c̄ .* (
            (-JacobiElliptic.dn(a, b) * JacobiElliptic.sn(a, b)),
            inv(2 * (1 - b) * b) *
            JacobiElliptic.dn(a, b) *
            JacobiElliptic.sn(a, b) *
            (
                (b - 1) * a + JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(a, b), b) -
                b * JacobiElliptic.cd(a, b) * JacobiElliptic.sn(a, b)
            ),
        )
)

#----------------------------------------------------------------------------------------
# Elliptic SN(a, b)
#----------------------------------------------------------------------------------------

@adjoint JacobiElliptic.sn(a, b) = JacobiElliptic.sn(a, b),
(
    c̄ ->
        c̄ .* (
            (JacobiElliptic.dn(a, b) * JacobiElliptic.cn(a, b)),
            inv(2 * (1 - b) * b) *
            JacobiElliptic.dn(a, b) *
            JacobiElliptic.cn(a, b) *
            (
                (1 - b) * a - JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(a, b), b) +
                b * JacobiElliptic.cd(a, b) * JacobiElliptic.sn(a, b)
            ),
        )
)


end # module
