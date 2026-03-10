module JacobiEllipticForwardDiffExt

using JacobiElliptic, ForwardDiff

function JacobiElliptic.CarlsonAlg._zero(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    zero(V)
end
function JacobiElliptic.CarlsonAlg._one(::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N}
    one(V)
end

#ForwardDiff.DiffRules.@define_diffrule JacobiElliptic.CarlsonAlg._sqrt(x) = :(inv(2 * JacobiElliptic.CarlsonAlg._sqrt($x)))
function JacobiElliptic.CarlsonAlg._sqrt(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg._sqrt(xval)
    ∂xf = inv(2 * fval)
    return ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic K(ϕ)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.K(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg.K(xval)
    ∂xf =
        (
            JacobiElliptic.CarlsonAlg.E(xval) -
            (1 - xval) * JacobiElliptic.CarlsonAlg.K(xval)
        ) / (2 * (1 - xval) * xval)
    ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic E(ϕ)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.E(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg.E(xval)
    ∂xf = iszero(xval) ? -T(π) / 8 :
        (fval - JacobiElliptic.CarlsonAlg.K(xval)) / (2xval)
    ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.E(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg.E(xval, y)
    ∂xf = sqrt(1 - y * sin(xval)^2)
    ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

function JacobiElliptic.CarlsonAlg.E(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    fval = JacobiElliptic.CarlsonAlg.E(x, yval)
    ∂yf = iszero(yval) ? (sin(2x) - 2x) / 8 :
        (fval - JacobiElliptic.CarlsonAlg.F(x, yval)) / (2yval)
    ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
end

function JacobiElliptic.CarlsonAlg.E(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = x.value
    yval = y.value

    fval = JacobiElliptic.CarlsonAlg.E(xval, yval)
    ∂xf = sqrt(1 - yval * sin(xval)^2)
    ∂yf =
        iszero(yval) ? (sin(2xval) - 2xval) / 8 :
        (fval - JacobiElliptic.CarlsonAlg.F(xval, yval)) / (2yval)
    ForwardDiff.Dual{T}(
        fval,
        ForwardDiff.Partials((∂xf * x.partials[1], ∂yf * y.partials[2])),
    )
end

#----------------------------------------------------------------------------------------
# Fukushima incomplete Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.FukushimaAlg.E(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value
    fval = JacobiElliptic.FukushimaAlg.E(xval, y)
    ∂xf = sqrt(1 - y * sin(xval)^2)

    ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

function JacobiElliptic.FukushimaAlg.E(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    fval = JacobiElliptic.FukushimaAlg.E(x, yval)
    ∂yf =
        iszero(yval) ? (sin(2x) - 2x) / 8 :
        (fval - JacobiElliptic.FukushimaAlg.F(x, yval)) / (2yval)

    ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
end

function JacobiElliptic.FukushimaAlg.E(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = x.value
    yval = y.value
    fval = JacobiElliptic.FukushimaAlg.E(xval, yval)
    ∂xf = sqrt(1 - yval * sin(xval)^2)
    ∂yf =
        iszero(yval) ? (sin(2xval) - 2xval) / 8 :
        (fval - JacobiElliptic.FukushimaAlg.F(xval, yval)) / (2yval)

    ForwardDiff.Dual{T}(fval, ∂xf * x.partials + ∂yf * y.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic cn(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.cn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg.cn(xval, y)
    ∂xf = -JacobiElliptic.dn(xval, y) * JacobiElliptic.sn(xval, y)

    ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

function JacobiElliptic.CarlsonAlg.cn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    fval = JacobiElliptic.CarlsonAlg.cn(x, yval)
    ∂yf =
        inv(2 * (1 - yval) * yval) *
        JacobiElliptic.dn(x, yval) *
        JacobiElliptic.sn(x, yval) *
        (
            (yval - 1) * x + JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(x, yval), yval) -
            yval * JacobiElliptic.cd(x, yval) * JacobiElliptic.sn(x, yval)
        )
    ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
end

function JacobiElliptic.CarlsonAlg.cn(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = ForwardDiff.value(x)
    yval = ForwardDiff.value(y)

    fval = JacobiElliptic.cn(xval, yval)
    ∂xf = -JacobiElliptic.dn(xval, yval) * JacobiElliptic.sn(xval, yval)
    ∂yf =
        inv(2 * (1 - yval) * yval) *
        JacobiElliptic.dn(xval, yval) *
        JacobiElliptic.sn(xval, yval) *
        (
            (yval - 1) * xval +
            JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(xval, yval), yval) -
            yval * JacobiElliptic.cd(xval, yval) * JacobiElliptic.sn(xval, yval)
        )
    ForwardDiff.Dual{T}(fval, ∂xf * ForwardDiff.partials(x) + ∂yf * ForwardDiff.partials(y))
end

#----------------------------------------------------------------------------------------
# Elliptic sn(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.sn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value
    fval = JacobiElliptic.CarlsonAlg.sn(xval, y)
    ∂xf = JacobiElliptic.dn(xval, y) * JacobiElliptic.cn(xval, y)

    ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
end

function JacobiElliptic.CarlsonAlg.sn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    fval = JacobiElliptic.CarlsonAlg.sn(x, yval)
    ∂yf =
        inv(2 * (1 - yval) * yval) *
        JacobiElliptic.dn(x, yval) *
        JacobiElliptic.cn(x, yval) *
        (
            (1 - yval) * x - JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(x, yval), yval) +
            yval * JacobiElliptic.cd(x, yval) * JacobiElliptic.sn(x, yval)
        )
    ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
end

function JacobiElliptic.CarlsonAlg.sn(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = x.value
    yval = y.value

    fval = JacobiElliptic.CarlsonAlg.sn(xval, yval)
    ∂xf = JacobiElliptic.dn(xval, yval) * JacobiElliptic.cn(xval, yval)
    ∂yf =
        inv(2 * (1 - yval) * yval) *
        JacobiElliptic.dn(xval, yval) *
        JacobiElliptic.cn(xval, yval) *
        (
            (1 - yval) * xval -
            JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(xval, yval), yval) +
            yval * JacobiElliptic.cd(xval, yval) * JacobiElliptic.sn(xval, yval)
        )

    ForwardDiff.Dual{T}(fval, ∂xf * x.partials + ∂yf * y.partials)
end

end # module
