module JacobiEllipticForwardDiffExt

using JacobiElliptic, ForwardDiff
include("common.jl")

function JacobiElliptic.CarlsonAlg._zero(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    zero(V)
end

function JacobiElliptic.CarlsonAlg._one(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    one(V)
end

_types = [ForwardDiff.Dual, Real]
for (T1, T2) in zip(_types, _types)
    function JacobiElliptic.CarlsonAlg._isequals(a::T1, b::T2)
        return ForwardDiff.value(a) == ForwardDiff.value(b)
    end
end

#ForwardDiff.DiffRules.@define_diffrule JacobiElliptic.CarlsonAlg._sqrt(x) = :(inv(2 * JacobiElliptic.CarlsonAlg._sqrt($x)))
function JacobiElliptic.CarlsonAlg._sqrt(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg._sqrt(xval)
    ∂xf = inv(2 * fval)
    return ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

begin
    alg = JacobiElliptic.ArithmeticGeometricMeanAlg
    #----------------------------------------------------------------------------------------
    # Elliptic K(ϕ)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).K(x::ForwardDiff.Dual{T}) where {T}
        xval = x.value
        fval = ($alg).K(xval)
        ∂xf = (($alg).E(xval) - (1 - xval) * ($alg).K(xval)) / (2 * (1 - xval) * xval)
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
    end

    #----------------------------------------------------------------------------------------
    # Elliptic E(ϕ)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).E(x::ForwardDiff.Dual{T}) where {T}
        xval = x.value
        fval = ($alg).E(xval)
        ∂xf = (($alg).E(xval) - ($alg).K(xval)) / (2xval)
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
    end

end

for alg in [JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg]
    #----------------------------------------------------------------------------------------
    # Elliptic K(ϕ)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).K(x::ForwardDiff.Dual{T}) where {T}
        xval = x.value
        fval = ($alg).K(xval)
        ∂xf = (($alg).E(xval) - (1 - xval) * ($alg).K(xval)) / (2 * (1 - xval) * xval)
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
    end

    #----------------------------------------------------------------------------------------
    # Elliptic E(ϕ)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).E(x::ForwardDiff.Dual{T}) where {T}
        xval = x.value
        fval = ($alg).E(xval)
        ∂xf = (($alg).E(xval) - ($alg).K(xval)) / (2xval)
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
    end
    #----------------------------------------------------------------------------------------
    # Elliptic F(ϕ, m)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).F(x::ForwardDiff.Dual{T}, y::U) where {T,U}
        xval = x.value
        fval = ($alg).F(xval, y)
        ∂xf = inv(sqrt(1 - y * sin(xval)^2))
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
    end

    @eval function ($alg).F(x::U, y::ForwardDiff.Dual{T}) where {T,U}
        yval = y.value
        fval = ($alg).F(x, yval)
        if iszero(yval)
            ∂yf = (2x - sin(2x)) / 8
        elseif isone(yval)
            sec_x = sec(x)
            tan_x = tan(x)
            ∂yf = 0.25 * (sec_x * tan_x - log(abs(sec_x + tan_x)))
        else
            ∂yf =
                ($alg).E(x, yval) / (2 * yval * (1 - yval)) - fval / (2 * yval) -
                sin(2 * x) / (4 * (1 - yval) * sqrt(1 - yval * sin(x)^2))
        end
        ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
    end

    @eval function ($alg).F(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = x.value
        yval = y.value
        fval = ($alg).F(xval, yval)
        sin_x = sin(xval)
        sin_2x = sin(2 * xval)
        sqrt_term = sqrt(1 - yval * sin_x^2)
        ∂xf = inv(sqrt_term)
        if iszero(yval)
            ∂yf = (2xval - sin_2x) / 8
        elseif isone(yval)
            sec_x = sec(xval)
            tan_x = tan(xval)
            ∂yf = 0.25 * (sec_x * tan_x - log(abs(sec_x + tan_x)))
        else
            ∂yf =
                ($alg).E(xval, yval) / (2 * yval * (1 - yval)) - fval / (2 * yval) -
                sin_2x / (4 * (1 - yval) * sqrt_term)
        end
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials + ∂yf * y.partials)
    end

    #----------------------------------------------------------------------------------------
    # Elliptic E(ϕ, m)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).E(x::ForwardDiff.Dual{T}, y::U) where {T,U}
        xval = x.value
        fval = ($alg).E(xval, y)
        ∂xf = sqrt(1 - y * sin(xval)^2)
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
    end

    @eval function ($alg).E(x::U, y::ForwardDiff.Dual{T}) where {T,U}
        yval = y.value
        fval = ($alg).E(x, yval)
        ∂yf = (fval - ($alg).F(x, yval)) / (2yval)
        ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
    end

    @eval function ($alg).E(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = x.value
        yval = y.value

        fval = ($alg).E(xval, yval)
        sin_x = sin(xval)
        sin_2x = sin(2 * xval)
        sqrt_term = sqrt(1 - yval * sin_x^2)
        ∂xf = sqrt_term
        ∂yf = iszero(yval) ? (sin_2x - 2xval) / 8 : (fval - ($alg).F(xval, yval)) / (2yval)
        ForwardDiff.Dual{T}(fval, ∂xf * x.partials + ∂yf * y.partials)
    end

    #----------------------------------------------------------------------------------------
    # Associate complete elliptic integral J(n, m)
    #----------------------------------------------------------------------------------------
    for mask ∈ 1:3
        n_is_dual = mask & 1 != 0
        m_is_dual = mask & 2 != 0
        n_type = n_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        m_type = m_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        n_value = n_is_dual ? :(ForwardDiff.value(n)) : :n
        m_value = m_is_dual ? :(ForwardDiff.value(m)) : :m
        tangent_terms = Any[]
        n_is_dual && push!(tangent_terms, :(∂n_j * ForwardDiff.partials(n)))
        m_is_dual && push!(tangent_terms, :(∂m_j * ForwardDiff.partials(m)))
        tangent = reduce((a, b) -> :($a + $b), tangent_terms)

        @eval function ($alg).J(n::$n_type, m::$m_type) where {T}
            n_value = $n_value
            m_value = $m_value
            j, ∂n_j, ∂m_j = _complete_j_derivatives(
                ($alg).E,
                ($alg).K,
                ($alg).Pi,
                ($alg).J,
                n_value,
                m_value,
            )
            return ForwardDiff.Dual{T}(j, $tangent)
        end
    end

    #----------------------------------------------------------------------------------------
    # Associate incomplete elliptic integral J(n, φ, m)
    #----------------------------------------------------------------------------------------
    for mask ∈ 1:7
        n_is_dual = mask & 1 != 0
        φ_is_dual = mask & 2 != 0
        m_is_dual = mask & 4 != 0
        n_type = n_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        φ_type = φ_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        m_type = m_is_dual ? :(ForwardDiff.Dual{T}) : :Real
        n_value = n_is_dual ? :(ForwardDiff.value(n)) : :n
        φ_value = φ_is_dual ? :(ForwardDiff.value(φ)) : :φ
        m_value = m_is_dual ? :(ForwardDiff.value(m)) : :m
        tangent_terms = Any[]
        n_is_dual && push!(tangent_terms, :(∂n_j * ForwardDiff.partials(n)))
        φ_is_dual && push!(tangent_terms, :(∂φ_j * ForwardDiff.partials(φ)))
        m_is_dual && push!(tangent_terms, :(∂m_j * ForwardDiff.partials(m)))
        tangent = reduce((a, b) -> :($a + $b), tangent_terms)

        @eval function ($alg).J(n::$n_type, φ::$φ_type, m::$m_type) where {T}
            n_value = $n_value
            φ_value = $φ_value
            m_value = $m_value
            j, ∂n_j, ∂φ_j, ∂m_j = _incomplete_j_derivatives(
                ($alg).E,
                ($alg).F,
                ($alg).Pi,
                ($alg).J,
                n_value,
                φ_value,
                m_value,
            )
            return ForwardDiff.Dual{T}(j, $tangent)
        end
    end

    #----------------------------------------------------------------------------------------
    # Elliptic cn(u, m)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).cn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
        xval = ForwardDiff.value(x)
        fval = ($alg).cn(xval, y)
        ∂xf = -($alg).dn(xval, y) * ($alg).sn(xval, y)
        ForwardDiff.Dual{T}(fval, ∂xf * ForwardDiff.partials(x))
    end

    @eval function ($alg).cn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
        yval = ForwardDiff.value(y)
        fval = ($alg).cn(x, yval)
        s = ($alg).sn(x, yval)
        d = ($alg).dn(x, yval)
        am = ($alg).am(x, yval)
        E_am = ($alg).E(am, yval)
        cd = fval/d
        ∂yf = _cn_parameter_derivative(E_am, cd, x, yval, s, fval, d)
        ForwardDiff.Dual{T}(fval, ∂yf * ForwardDiff.partials(y))
    end

    @eval function ($alg).cn(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = ForwardDiff.value(x)
        yval = ForwardDiff.value(y)
        fval = ($alg).cn(xval, yval)
        s = ($alg).sn(xval, yval)
        d = ($alg).dn(xval, yval)
        am = ($alg).am(xval, yval)
        E_am = ($alg).E(am, yval)
        cd = fval/d
        ∂xf = -d * s
        ∂yf = _cn_parameter_derivative(E_am, cd, xval, yval, s, fval, d)
        ForwardDiff.Dual{T}(
            fval,
            ∂xf * ForwardDiff.partials(x) + ∂yf * ForwardDiff.partials(y),
        )
    end

    #----------------------------------------------------------------------------------------
    # Elliptic sn(u, m)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).sn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
        xval = ForwardDiff.value(x)
        fval = ($alg).sn(xval, y)
        ∂xf = ($alg).dn(xval, y) * ($alg).cn(xval, y)
        ForwardDiff.Dual{T}(fval, ∂xf * ForwardDiff.partials(x))
    end

    @eval function ($alg).sn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
        yval = ForwardDiff.value(y)
        fval = ($alg).sn(x, yval)
        c = ($alg).cn(x, yval)
        d = ($alg).dn(x, yval)
        am = ($alg).am(x, yval)
        E_am = ($alg).E(am, yval)
        cd = c/d
        ∂yf = _sn_parameter_derivative(E_am, cd, x, yval, fval, c, d)
        ForwardDiff.Dual{T}(fval, ∂yf * ForwardDiff.partials(y))
    end

    @eval function ($alg).sn(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = ForwardDiff.value(x)
        yval = ForwardDiff.value(y)
        fval = ($alg).sn(xval, yval)
        c = ($alg).cn(xval, yval)
        d = ($alg).dn(xval, yval)
        am = ($alg).am(xval, yval)
        E_am = ($alg).E(am, yval)
        cd = c/d
        ∂xf = d * c
        ∂yf = _sn_parameter_derivative(E_am, cd, xval, yval, fval, c, d)
        ForwardDiff.Dual{T}(
            fval,
            ∂xf * ForwardDiff.partials(x) + ∂yf * ForwardDiff.partials(y),
        )
    end

    #----------------------------------------------------------------------------------------
    # Elliptic dn(u, m)
    #----------------------------------------------------------------------------------------
    @eval function ($alg).dn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
        xval = ForwardDiff.value(x)
        fval = ($alg).dn(xval, y)
        ∂xf = -y * ($alg).sn(xval, y) * ($alg).cn(xval, y)
        ForwardDiff.Dual{T}(fval, ∂xf * ForwardDiff.partials(x))
    end

    @eval function ($alg).dn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
        yval = ForwardDiff.value(y)
        fval = ($alg).dn(x, yval)
        s = ($alg).sn(x, yval)
        c = ($alg).cn(x, yval)
        am = ($alg).am(x, yval)
        E_am = ($alg).E(am, yval)
        cd = c/fval
        ∂yf = _dn_parameter_derivative(E_am, cd, x, yval, s, c, fval)
        ForwardDiff.Dual{T}(fval, ∂yf * ForwardDiff.partials(y))
    end

    @eval function ($alg).dn(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = ForwardDiff.value(x)
        yval = ForwardDiff.value(y)
        fval = ($alg).dn(xval, yval)
        s = ($alg).sn(xval, yval)
        c = ($alg).cn(xval, yval)
        am = ($alg).am(xval, yval)
        E_am = ($alg).E(am, yval)
        cd = c/fval
        ∂xf = -yval * s * c
        ∂yf = _dn_parameter_derivative(E_am, cd, xval, yval, s, c, fval)
        ForwardDiff.Dual{T}(
            fval,
            ∂xf * ForwardDiff.partials(x) + ∂yf * ForwardDiff.partials(y),
        )
    end


end
end # module
