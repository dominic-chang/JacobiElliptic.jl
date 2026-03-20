module JacobiEllipticForwardDiffExt

using JacobiElliptic, ForwardDiff

function JacobiElliptic.CarlsonAlg._zero(::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N}
	zero(V)
end

function JacobiElliptic.CarlsonAlg._one(::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N}
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

for alg in [JacobiElliptic.CarlsonAlg, JacobiElliptic.FukushimaAlg]
	#----------------------------------------------------------------------------------------
	# Elliptic K(ϕ)
	#----------------------------------------------------------------------------------------
	@eval function ($alg).K(x::ForwardDiff.Dual{T}) where {T}
		xval = x.value
		fval = ($alg).K(xval)
		∂xf =
			(
				($alg).E(xval) -
				(1 - xval) * ($alg).K(xval)
			) / (2 * (1 - xval) * xval)
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
	@eval function ($alg).F(x::ForwardDiff.Dual{T}, y::U) where {T, U}
		xval = x.value
		fval = ($alg).F(xval, y)
		∂xf = inv(sqrt(1 - y * sin(xval)^2))
		ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
	end

	@eval function ($alg).F(x::U, y::ForwardDiff.Dual{T}) where {T, U}
		yval = y.value
		fval = ($alg).F(x, yval)
		∂yf =
			iszero(yval) ? (2x - sin(2x)) / 8 :
			($alg).E(x, yval) / (2 * yval * (1 - yval)) - fval / (2 * yval) -
			sin(2 * x) / (4 * (1 - yval) * sqrt(1 - yval * sin(x)^2))
		ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
	end

	@eval function ($alg).F(
		x::ForwardDiff.Dual{T},
		y::ForwardDiff.Dual{T},
	) where {T}
		xval = x.value
		yval = y.value
		fval = ($alg).F(xval, yval)
		∂xf = inv(sqrt(1 - yval * sin(xval)^2))
		∂yf =
			iszero(yval) ? (2xval - sin(2xval)) / 8 :
			($alg).E(xval, yval) / (2 * yval * (1 - yval)) - fval / (2 * yval) -
			sin(2 * xval) / (4 * (1 - yval) * sqrt(1 - yval * sin(xval)^2))
		ForwardDiff.Dual{T}(fval, ∂xf * x.partials + ∂yf * y.partials)
	end

	#----------------------------------------------------------------------------------------
	# Elliptic E(ϕ, m)
	#----------------------------------------------------------------------------------------
	@eval function ($alg).E(x::ForwardDiff.Dual{T}, y::U) where {T, U}
		xval = x.value
		fval = ($alg).E(xval, y)
		∂xf = sqrt(1 - y * sin(xval)^2)
		ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
	end

	@eval function ($alg).E(x::U, y::ForwardDiff.Dual{T}) where {T, U}
		yval = y.value
		fval = ($alg).E(x, yval)
		∂yf = (fval - ($alg).F(x, yval)) / (2yval)
		ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
	end

	@eval function ($alg).E(
		x::ForwardDiff.Dual{T},
		y::ForwardDiff.Dual{T},
	) where {T}
		xval = x.value
		yval = y.value

		fval = ($alg).E(xval, yval)
		∂xf = sqrt(1 - yval * sin(xval)^2)
		∂yf =
			iszero(yval) ? (sin(2xval) - 2xval) / 8 :
			(fval - ($alg).F(xval, yval)) / (2yval)
		ForwardDiff.Dual{T}(
			fval,
			ForwardDiff.Partials((∂xf * x.partials[1], ∂yf * y.partials[2])),
		)
	end

	#----------------------------------------------------------------------------------------
	# Elliptic cn(ϕ, m)
	#----------------------------------------------------------------------------------------
	@eval function ($alg).cn(x::ForwardDiff.Dual{T}, y::U) where {T, U}
		xval = x.value
		fval = ($alg).cn(xval, y)
		∂xf = -JacobiElliptic.dn(xval, y) * JacobiElliptic.sn(xval, y)

		ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
	end

	@eval function ($alg).cn(x::U, y::ForwardDiff.Dual{T}) where {T, U}
		yval = y.value
		fval = ($alg).cn(x, yval)
		∂yf =
			inv(2 * (1 - yval) * yval) *
			JacobiElliptic.dn(x, yval) *
			JacobiElliptic.sn(x, yval) *
			(
				(yval - 1) * x + ($alg).E(JacobiElliptic.am(x, yval), yval) -
				yval * JacobiElliptic.cd(x, yval) * JacobiElliptic.sn(x, yval)
			)
		ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
	end

	@eval function ($alg).cn(
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
				($alg).E(JacobiElliptic.am(xval, yval), yval) -
				yval * JacobiElliptic.cd(xval, yval) * JacobiElliptic.sn(xval, yval)
			)
		ForwardDiff.Dual{T}(fval, ∂xf * ForwardDiff.partials(x) + ∂yf * ForwardDiff.partials(y))
	end

	#----------------------------------------------------------------------------------------
	# Elliptic sn(ϕ, m)
	#----------------------------------------------------------------------------------------
	@eval function ($alg).sn(x::ForwardDiff.Dual{T}, y::U) where {T, U}
		xval = x.value
		fval = ($alg).sn(xval, y)
		∂xf = JacobiElliptic.dn(xval, y) * JacobiElliptic.cn(xval, y)

		ForwardDiff.Dual{T}(fval, ∂xf * x.partials)
	end

	@eval function ($alg).sn(x::U, y::ForwardDiff.Dual{T}) where {T, U}
		yval = y.value
		fval = ($alg).sn(x, yval)
		∂yf =
			inv(2 * (1 - yval) * yval) *
			JacobiElliptic.dn(x, yval) *
			JacobiElliptic.cn(x, yval) *
			(
				(1 - yval) * x - ($alg).E(JacobiElliptic.am(x, yval), yval) +
				yval * JacobiElliptic.cd(x, yval) * JacobiElliptic.sn(x, yval)
			)
		ForwardDiff.Dual{T}(fval, ∂yf * y.partials)
	end

	@eval function ($alg).sn(
		x::ForwardDiff.Dual{T},
		y::ForwardDiff.Dual{T},
	) where {T}
		xval = x.value
		yval = y.value

		fval = ($alg).sn(xval, yval)
		∂xf = JacobiElliptic.dn(xval, yval) * JacobiElliptic.cn(xval, yval)
		∂yf =
			inv(2 * (1 - yval) * yval) *
			JacobiElliptic.dn(xval, yval) *
			JacobiElliptic.cn(xval, yval) *
			(
				(1 - yval) * xval -
				($alg).E(JacobiElliptic.am(xval, yval), yval) +
				yval * JacobiElliptic.cd(xval, yval) * JacobiElliptic.sn(xval, yval)
			)

		ForwardDiff.Dual{T}(fval, ∂xf * x.partials + ∂yf * y.partials)
	end


end
end # module
