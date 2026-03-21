export am

# Abramowitz & Stegun, section 16.4, p571
#const _ambuf = Array{Float64}(undef, 10)
function _am(u::A, m::B, tol::C) where {A,B,C}

    T = promote_type(A, B, C)

    #Pre-allocate with SVector to avoid GPU allocations
    _ambuf = StaticArrays.@MVector[
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
    ]
    u == 0 && return zero(T)

    sqrt_tol = sqrt(tol)
    if m < sqrt_tol
        # A&S 16.13.4
        return u - m * (u - sin(2 * u) / 2) / 4
    end
    m1 = 1 - m
    if m1 < sqrt_tol
        # A&S 16.15.4
        t = tanh(u)
        return asin(t) + m1 * (t - u * (1 - t^2)) * cosh(u) / 4
    end

    a, b, c, n = 1, sqrt(m1), sqrt(m), 0
    while abs(c) > tol
        @assert n < 10
        a, b, c, n = ((a + b) / 2, sqrt(a * b), (a - b) / 2, n + 1)
        @inbounds _ambuf[n] = c / a
    end

    #phi = ldexp(a*u, n) # Doesn't work with Enzyme
    phi = a * u * (2^n)
    for i = n:-1:1
        @inbounds phi = (phi + asin(_ambuf[i] * sin(phi))) / 2
    end
    phi
end
function _am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    return _am(u, m, eps(T))
end

"""
    am(u::Real, m::Real)

Returns amplitude, φ, such that u = F(φ | m)
"""
function am(u::A, m::B) where {A,B}
    T = promote_type(A, B)
    if m < 0
        mu1 = inv(1 - m)
        mu = -m * mu1
        sqrtmu1 = sqrt(mu1)
        v = u / sqrtmu1
        phi = _am(v, mu)
        s = sin(phi)
        t = floor((phi + A(π / 2)) / A(π))

        return t * π + ((-1)^t) * asin(sqrtmu1 * s / sqrt(1 - mu * s^2))
    elseif m <= 1 # 0 <= m <= 1
        return _am(u, m, eps(T))
    else # m > 1
        # Reciprocal Modulus Transformation
        # P.F. Byrd & M.D. Friedman: Handbook of Elliptic Integrals for Engineers and Scientists, 1971, p. 38, eq. 162.01
        # am(u,k) = arcsin(sin(am(u,k))) = arcsin(sn(u,k)) = arcsin(1 / k * sn(k * u, 1 / k))
        k = sqrt(m)
        return asin(inv(k) * sn(k * u, inv(m)))
    end
end

@inline function custom_atanh(a::T) where {T}
    oneT = one(T)
    arg1 = abs(oneT + a)
    arg2 = abs(oneT - a)

    ans = (log(arg1 / arg2)) / 2
    return ans
end

@inline function FukushimaT(t::A, h::B) where {A, B}
	T = promote_type(A, B)
	if h > zero(T)
		sqrt_h = √h
		return atan(t * sqrt_h) / sqrt_h
	elseif h == zero(T)
		return t
	else
		sqrt_neg_h = √(-h)
		arg = t * sqrt_neg_h
		ans = abs(arg) < one(T) ? atanh(arg) : custom_atanh(arg)
		return ans / sqrt_neg_h
	end
end

#https://link-springer-com.ezp-prod1.hul.harvard.edu/article/T(10).1007/BF02165405
function cel(kc::A, p::B, a::C, b::D) where {A, B, C, D}
	T = promote_type(A, B, C, D)
	#ca = T(1e-6)
	ca = eps(T)
	oneT = one(T)
	twoT = T(2)
	pi_over_2 = T(π / 2)
	kc = abs(kc)
	e = kc
	m = oneT

	f, g, q = T(0), T(0), T(0)
	if p > T(0)
		p = √p
		b = b / p
	else
		f = kc^2
		q = oneT - f
		g = oneT - p
		f = f - p
		q = (b - a * p) * q
		p = √(f / g)
		a = (a - b) / g
		b = -q * (g^2 * p) + a * p
	end
	count = 0
	while count < 1000
		count+=1
		f = a
		invp = inv(p)
		a = muladd(invp, b, a)
		g = e * invp
		b = twoT * muladd(f, g, b)
		p = g + p
		g = m
		m = kc + m
		if abs(g - kc) < g * ca
			break
		end
		kc = twoT * √e
		e = kc * m
	end
	return pi_over_2 * muladd(a, m, b) / (m * (m + p))
end

