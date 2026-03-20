module JacobiEllipticZygoteExt

using JacobiElliptic
using Zygote: @adjoint


@adjoint JacobiElliptic.CarlsonAlg._sqrt(a) = JacobiElliptic.CarlsonAlg._sqrt(a),
(c̄ -> (c̄ * inv(2 * JacobiElliptic.CarlsonAlg._sqrt(a)),))

for alg in (JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg)
    @eval begin
        #------------------------------------------------------------------------------------
        # Elliptic F(a, b)
        #------------------------------------------------------------------------------------

        @adjoint function ($alg).F(a, b)
            F = ($alg).F(a, b)
            E = ($alg).E(a, b)
            sin_a = sin(a)
            sin_2a = sin(2 * a)
            sqrt_term = sqrt(1 - b * sin_a^2)
            return F,
            c̄ -> (
                c̄ / sqrt_term,
                c̄ * (
                    iszero(b) ? (2 * a - sin_2a) / 8 :
                    E / (2 * b * (1 - b)) - F / (2 * b) -
                    sin_2a / (4 * (1 - b) * sqrt_term)
                ),
            )
        end

        #------------------------------------------------------------------------------------
        # Elliptic E(a, b)
        #------------------------------------------------------------------------------------

        @adjoint function ($alg).E(a, b)
            E = ($alg).E(a, b)
            F = ($alg).F(a, b)
            sin_a = sin(a)
            return E,
            c̄ -> (
                c̄ * sqrt(1 - b * sin_a^2),
                c̄ * (iszero(b) ? (sin(2 * a) - 2 * a) / 8 : (E - F) / (2 * b)),
            )
        end

        #------------------------------------------------------------------------------------
        # Elliptic Pi(a, b)
        #------------------------------------------------------------------------------------

        @adjoint function ($alg).Pi(a, b)
            Pi = ($alg).Pi(a, b)
            K = ($alg).K(b)
            E_b = ($alg).E(b)
            return Pi,
            c̄ -> (
                c̄ * (
                    iszero(a) ? (iszero(b) ? π / 4 : (K - E_b) / b) :
                    (E_b + (b - a) * K / a + (a^2 - b) * Pi / a) / (2 * (b - a) * (a - 1))
                ),
                c̄ * ((E_b / (b - 1) + Pi) / (2 * (a - b))),
            )
        end

        #------------------------------------------------------------------------------------
        # Elliptic Pi(a, b, c)
        #------------------------------------------------------------------------------------

        @adjoint function ($alg).Pi(a, b, c)
            Pi = ($alg).Pi(a, b, c)
            F_bc = ($alg).F(b, c)
            E_bc = ($alg).E(b, c)
            sin_b = sin(b)
            sin_2b = sin(2 * b)
            sqrt_term = sqrt(1 - c * sin_b^2)
            return Pi,
            c̄ -> (
                c̄ * (
                    iszero(a) ? (iszero(c) ? (2 * b - sin_2b) / 4 : (F_bc - E_bc) / c) :
                    (
                        E_bc + (c - a) * F_bc / a + (a^2 - c) * Pi / a -
                        a * sqrt_term * sin_2b / (2 * (1 - a * sin_b^2))
                    ) / (2 * (c - a) * (a - 1))
                ),
                c̄ * (1 / (sqrt_term * (1 - a * sin_b^2))),
                c̄ * (
                    (E_bc / (c - 1) + Pi - c * sin_2b / (2 * (c - 1) * sqrt_term)) /
                    (2 * (a - c))
                ),
            )
        end

        #------------------------------------------------------------------------------------
        # Elliptic CN(a, b)
        #------------------------------------------------------------------------------------

        @adjoint function ($alg).cn(a, b)
            cn = ($alg).cn(a, b)
            dn = ($alg).dn(a, b)
            sn = ($alg).sn(a, b)
            am = ($alg).am(a, b)
            cd = ($alg).cd(a, b)
            E_am = ($alg).E(am, b)
            return cn,
            c̄ -> (
                c̄ * (-dn * sn),
                c̄ *
                inv(2 * (1 - b) * b) *
                dn *
                sn *
                ((b - 1) * a + E_am - b * (cn / dn) * sn),
            )
        end

        #------------------------------------------------------------------------------------
        # Elliptic SN(a, b)
        #------------------------------------------------------------------------------------

        @adjoint function ($alg).sn(a, b)
            sn = ($alg).sn(a, b)
            dn = ($alg).dn(a, b)
            cn = ($alg).cn(a, b)
            am = ($alg).am(a, b)
            cd = ($alg).cd(a, b)
            E_am = ($alg).E(am, b)
            return sn,
            c̄ -> (
                c̄ * (dn * cn),
                c̄ *
                inv(2 * (1 - b) * b) *
                dn *
                cn *
                ((1 - b) * a - E_am + b * (cn / dn) * sn),
            )
        end
    end
end


end # module
