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
    return ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg._sqrt(xval), inv(2 * JacobiElliptic.CarlsonAlg._sqrt(xval)) * x.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic K(ϕ)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.K(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value

    ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg.K(xval), (JacobiElliptic.CarlsonAlg.E(xval) - (1-xval)*JacobiElliptic.CarlsonAlg.K(xval)) / (2*(1-xval)*xval) * x.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic E(ϕ)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.E(x::ForwardDiff.Dual{T}) where {T}
    xval = x.value

    ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg.E(xval), (JacobiElliptic.CarlsonAlg.E(xval) - JacobiElliptic.CarlsonAlg.K(xval)) / (2xval) * x.partials)
end

#----------------------------------------------------------------------------------------
# Elliptic E(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.E(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value

    ForwardDiff.Dual{T}(JacobiElliptic.CarlsonAlg.E(xval, y), (sqrt(1 - y * sin(xval)^2)) * x.partials)
end

function JacobiElliptic.CarlsonAlg.E(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.E(x, yval),
        (JacobiElliptic.CarlsonAlg.E(x, yval) - JacobiElliptic.CarlsonAlg.F(x, yval)) /
        (2yval) * y.partials,
    )
end

function JacobiElliptic.CarlsonAlg.E(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = x.value
    yval = y.value

    ∂yE =
        iszero(yval) ? -π / 8 :
        (
            JacobiElliptic.CarlsonAlg.E(xval, yval) -
            JacobiElliptic.CarlsonAlg.F(xval, yval)
        ) / (2yval)
    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.E(xval, yval),
        ForwardDiff.Partials(((sqrt(1 - yval * sin(xval)^2)) .* x.partials[1], ∂yE * y.partials[2],))
    )
end

#----------------------------------------------------------------------------------------
# Elliptic cn(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.cn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value

    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.cn(xval, y),
        -JacobiElliptic.dn(xval, y) * JacobiElliptic.sn(xval, y) * x.partials,
    )
end

function JacobiElliptic.CarlsonAlg.cn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.cn(x, yval),
        inv(2 * (1 - yval) * yval) *
        JacobiElliptic.dn(x, yval) *
        JacobiElliptic.sn(x, yval) *
        (
            (yval - 1) * x + JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(x, yval), yval) -
            yval * JacobiElliptic.cd(x, yval) * JacobiElliptic.sn(x, yval)
        ) * y.partials,
    )
end

function JacobiElliptic.CarlsonAlg.cn(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = x.value
    yval = y.value

    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.cn(xval, yval),
        ForwardDiff.Partials((
            -JacobiElliptic.dn(xval, m) * JacobiElliptic.sn(xval, m) * x.partials[1],
            inv(2 * (1 - m) * m) *
            JacobiElliptic.dn(xval, m) *
            JacobiElliptic.sn(xval, m) *
            (
                (yval - 1) * xval +
                JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(xval, yval), yval) -
                yval * JacobiElliptic.cd(xval, yval) * JacobiElliptic.sn(xval, yval)
            ) * y.partials[2],
        )),
    )
end

#----------------------------------------------------------------------------------------
# Elliptic sn(ϕ, m)
#----------------------------------------------------------------------------------------
function JacobiElliptic.CarlsonAlg.sn(x::ForwardDiff.Dual{T}, y::U) where {T,U}
    xval = x.value

    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.sn(xval, y),
        JacobiElliptic.dn(xval, y) * JacobiElliptic.cn(xval, y) * x.partials,
    )
end

function JacobiElliptic.CarlsonAlg.sn(x::U, y::ForwardDiff.Dual{T}) where {T,U}
    yval = y.value
    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.sn(x, yval),
        inv(2 * (1 - yval) * yval) *
        JacobiElliptic.dn(x, yval) *
        JacobiElliptic.cn(x, yval) *
        (
            (1 - yval) * x - JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(x, yval), yval) +
            yval * JacobiElliptic.cd(x, yval) * JacobiElliptic.sn(x, yval)
        ) * y.partials,
    )
end

function JacobiElliptic.CarlsonAlg.sn(
    x::ForwardDiff.Dual{T},
    y::ForwardDiff.Dual{T},
) where {T}
    xval = x.value
    yval = y.value

    ForwardDiff.Dual{T}(
        JacobiElliptic.CarlsonAlg.sn(xval, yval),
        ForwardDiff.Partials((
            JacobiElliptic.dn(xval, m) * JacobiElliptic.cn(xval, m) * x.partials[1],
            inv(2 * (1 - m) * m) *
            JacobiElliptic.dn(xval, m) *
            JacobiElliptic.cn(xval, m) *
            (
                (1 - yval) * xval -
                JacobiElliptic.CarlsonAlg.E(JacobiElliptic.am(xval, yval), yval) +
                yval * JacobiElliptic.cd(xval, yval) * JacobiElliptic.sn(xval, yval)
            ) * y.partials[2],
        )),
    )
end

end # module
