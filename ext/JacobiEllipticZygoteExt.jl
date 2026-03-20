module JacobiEllipticZygoteExt

using JacobiElliptic
using Zygote: @adjoint


@adjoint JacobiElliptic.CarlsonAlg._sqrt(a) = JacobiElliptic.CarlsonAlg._sqrt(a),
(c̄ -> (c̄ * inv(2 * JacobiElliptic.CarlsonAlg._sqrt(a)),))

for alg in (JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg)
    @eval begin
        #------------------------------------------------------------------------------------
        # Elliptic E(a, b)
        #------------------------------------------------------------------------------------

        @adjoint ($alg).E(a, b) = ($alg).E(a, b),
        (
            c̄ -> (
                c̄ * sqrt(1 - b * sin(a)^2),
                c̄ * (
                    iszero(b) ? (sin(2 * a) - 2 * a) / 8 :
                    (($alg).E(a, b) - ($alg).F(a, b)) / (2b)
                ),
            )
        )

        #------------------------------------------------------------------------------------
        # Elliptic Pi(a, b)
        #------------------------------------------------------------------------------------

        @adjoint ($alg).Pi(a, b) = ($alg).Pi(a, b),
        (
            c̄ -> (
                c̄ * (
                    iszero(a) ?
                    (iszero(b) ? π / 4 : (($alg).K(b) - ($alg).E(b)) / b) :
                    (
                        ($alg).E(b) +
                        (b - a) * ($alg).K(b) / a +
                        (a^2 - b) * ($alg).Pi(a, b) / a
                    ) / (2 * (b - a) * (a - 1))
                ),
                c̄ * ((($alg).E(b) / (b - 1) + ($alg).Pi(a, b)) / (2 * (a - b))),
            )
        )

        #------------------------------------------------------------------------------------
        # Elliptic Pi(a, b, c)
        #------------------------------------------------------------------------------------

        @adjoint ($alg).Pi(a, b, c) = ($alg).Pi(a, b, c),
        (
            c̄ -> (
                c̄ * (
                    iszero(a) ?
                    (
                        iszero(c) ? (2 * b - sin(2 * b)) / 4 :
                        (($alg).F(b, c) - ($alg).E(b, c)) / c
                    ) :
                    (
                        ($alg).E(b, c) +
                        (c - a) * ($alg).F(b, c) / a +
                        (a^2 - c) * ($alg).Pi(a, b, c) / a -
                        a * sqrt(1 - c * sin(b)^2) * sin(2 * b) / (2 * (1 - a * sin(b)^2))
                    ) / (2 * (c - a) * (a - 1))
                ),
                c̄ * (1 / (sqrt(1 - c * sin(b)^2) * (1 - a * sin(b)^2))),
                c̄ * (
                    (
                        ($alg).E(b, c) / (c - 1) + ($alg).Pi(a, b, c) -
                        c * sin(2 * b) / (2 * (c - 1) * sqrt(1 - c * sin(b)^2))
                    ) / (2 * (a - c))
                ),
            )
        )

        #------------------------------------------------------------------------------------
        # Elliptic CN(a, b)
        #------------------------------------------------------------------------------------

        @adjoint ($alg).cn(a, b) = ($alg).cn(a, b),
        (
            c̄ -> (
                c̄ * (-($alg).dn(a, b) * ($alg).sn(a, b)),
                c̄ *
                inv(2 * (1 - b) * b) *
                ($alg).dn(a, b) *
                ($alg).sn(a, b) *
                (
                    (b - 1) * a + ($alg).E(($alg).am(a, b), b) -
                    b * ($alg).cd(a, b) * ($alg).sn(a, b)
                ),
            )
        )

        #------------------------------------------------------------------------------------
        # Elliptic SN(a, b)
        #------------------------------------------------------------------------------------

        @adjoint ($alg).sn(a, b) = ($alg).sn(a, b),
        (
            c̄ -> (
                c̄ * (($alg).dn(a, b) * ($alg).cn(a, b)),
                c̄ *
                inv(2 * (1 - b) * b) *
                ($alg).dn(a, b) *
                ($alg).cn(a, b) *
                (
                    (1 - b) * a - ($alg).E(($alg).am(a, b), b) +
                    b * ($alg).cd(a, b) * ($alg).sn(a, b)
                ),
            )
        )
    end
end


end # module
