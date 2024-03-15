using Revise
using JacobiElliptic
using BenchmarkTools

@btime JacobiElliptic.Pi(JacobiElliptic.Fukushima(), rand(), rand(), rand())
@btime JacobiElliptic.Pi(JacobiElliptic.Carlson(), rand(), rand(), rand())

@btime JacobiElliptic.F(JacobiElliptic.Fukushima(), rand(), rand())
@btime JacobiElliptic.F(JacobiElliptic.Carlson(), rand(), rand())

@btime JacobiElliptic.E(JacobiElliptic.Fukushima(), rand(), rand())
@btime JacobiElliptic.E(JacobiElliptic.Carlson(), rand(), rand())

@btime JacobiElliptic.K(JacobiElliptic.Fukushima(), rand())
@btime JacobiElliptic.K(JacobiElliptic.Carlson(), rand())