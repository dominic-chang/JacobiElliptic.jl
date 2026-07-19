module JacobiEllipticForwardDiffExt

using JacobiElliptic, ForwardDiff

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

function _sn_parameter_derivative(E, am, cd, u, m, s, c, d)
    if iszero(m)
        return -c * (u - s * c) / 4
    elseif isone(m)
        return -(s - u * c^2) / 4
    end

    return d * c * ((1 - m) * u - E(am(u, m), m) + m * cd(u, m) * s) / (2 * (1 - m) * m)
end

function _cn_parameter_derivative(E, am, cd, u, m, s, c, d)
    if iszero(m)
        return s * (u - s * c) / 4
    elseif isone(m)
        return s * (s - u * c^2) / (4c)
    end

    return d * s * ((m - 1) * u + E(am(u, m), m) - m * cd(u, m) * s) / (2 * (1 - m) * m)
end

function _dn_parameter_derivative(E, am, cd, u, m, s, c, d)
    iszero(m) && return -s^2 / 2
    ∂m_sn = _sn_parameter_derivative(E, am, cd, u, m, s, c, d)
    return -(s^2 + 2m * s * ∂m_sn) / (2d)
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
        ∂yf =
            iszero(yval) ? (2x - sin(2x)) / 8 :
            ($alg).E(x, yval) / (2 * yval * (1 - yval)) - fval / (2 * yval) -
            sin(2 * x) / (4 * (1 - yval) * sqrt(1 - yval * sin(x)^2))
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
        ∂yf =
            iszero(yval) ? (2xval - sin_2x) / 8 :
            ($alg).E(xval, yval) / (2 * yval * (1 - yval)) - fval / (2 * yval) -
            sin_2x / (4 * (1 - yval) * sqrt_term)
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
        ∂yf = _cn_parameter_derivative(($alg).E, ($alg).am, ($alg).cd, x, yval, s, fval, d)
        ForwardDiff.Dual{T}(fval, ∂yf * ForwardDiff.partials(y))
    end

    @eval function ($alg).cn(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = ForwardDiff.value(x)
        yval = ForwardDiff.value(y)
        fval = ($alg).cn(xval, yval)
        s = ($alg).sn(xval, yval)
        d = ($alg).dn(xval, yval)
        ∂xf = -d * s
        ∂yf =
            _cn_parameter_derivative(($alg).E, ($alg).am, ($alg).cd, xval, yval, s, fval, d)
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
        ∂yf = _sn_parameter_derivative(($alg).E, ($alg).am, ($alg).cd, x, yval, fval, c, d)
        ForwardDiff.Dual{T}(fval, ∂yf * ForwardDiff.partials(y))
    end

    @eval function ($alg).sn(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = ForwardDiff.value(x)
        yval = ForwardDiff.value(y)
        fval = ($alg).sn(xval, yval)
        c = ($alg).cn(xval, yval)
        d = ($alg).dn(xval, yval)
        ∂xf = d * c
        ∂yf =
            _sn_parameter_derivative(($alg).E, ($alg).am, ($alg).cd, xval, yval, fval, c, d)
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
        ∂yf = _dn_parameter_derivative(($alg).E, ($alg).am, ($alg).cd, x, yval, s, c, fval)
        ForwardDiff.Dual{T}(fval, ∂yf * ForwardDiff.partials(y))
    end

    @eval function ($alg).dn(x::ForwardDiff.Dual{T}, y::ForwardDiff.Dual{T}) where {T}
        xval = ForwardDiff.value(x)
        yval = ForwardDiff.value(y)
        fval = ($alg).dn(xval, yval)
        s = ($alg).sn(xval, yval)
        c = ($alg).cn(xval, yval)
        ∂xf = -yval * s * c
        ∂yf =
            _dn_parameter_derivative(($alg).E, ($alg).am, ($alg).cd, xval, yval, s, c, fval)
        ForwardDiff.Dual{T}(
            fval,
            ∂xf * ForwardDiff.partials(x) + ∂yf * ForwardDiff.partials(y),
        )
    end


end
end # module
