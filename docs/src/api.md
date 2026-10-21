# JacobiElliptic API

## Index
```@index
Pages = ["api.md"]
```

## Algorithm Selection

Every function below dispatches to a default algorithm. Pass an
[`JacobiElliptic.AbstractAlgorithm`](@ref) instance as the first argument to choose one
explicitly, for example `F(Carlson(), φ, m)` or `K(Fukushima(), m)`.

```@docs
JacobiElliptic.default_algorithm
JacobiElliptic.set_default_algorithm!
JacobiElliptic.AbstractAlgorithm
JacobiElliptic.AGM
JacobiElliptic.Carlson
JacobiElliptic.Fukushima
```

## Elliptic Integrals

```@docs
JacobiElliptic.K
JacobiElliptic.E
JacobiElliptic.F
JacobiElliptic.Pi
JacobiElliptic.J
```

### Per-algorithm implementations

```@docs
JacobiElliptic.ArithmeticGeometricMeanAlg.K
JacobiElliptic.ArithmeticGeometricMeanAlg.E
JacobiElliptic.ArithmeticGeometricMeanAlg.ellipke
JacobiElliptic.CarlsonAlg.ellipke
JacobiElliptic.CarlsonAlg.F
JacobiElliptic.CarlsonAlg.E
JacobiElliptic.CarlsonAlg.Pi
JacobiElliptic.FukushimaAlg.K
JacobiElliptic.FukushimaAlg.E
JacobiElliptic.FukushimaAlg.J
```

## Jacobi Elliptic Functions

```@docs
JacobiElliptic.am
JacobiElliptic.sn
JacobiElliptic.cn
JacobiElliptic.dn
```

Note that `cd` is not exported, to avoid clashing
with `Base.cd`.

```@docs
JacobiElliptic.sd
JacobiElliptic.sc
JacobiElliptic.cd
JacobiElliptic.cs
JacobiElliptic.dc
JacobiElliptic.ds
JacobiElliptic.nd
JacobiElliptic.nc
JacobiElliptic.ns
JacobiElliptic.nn
JacobiElliptic.dd
JacobiElliptic.cc
JacobiElliptic.ss
```

### Per-algorithm implementations

```@docs
JacobiElliptic.CarlsonAlg.am
JacobiElliptic.FukushimaAlg.am
JacobiElliptic.FukushimaAlg.sn
JacobiElliptic.FukushimaAlg.cn
JacobiElliptic.FukushimaAlg.dn
```

## Inverse Jacobi Elliptic Functions

Only [`JacobiElliptic.Fukushima`](@ref) provides these.

```@docs
JacobiElliptic.FukushimaAlg.asn
JacobiElliptic.FukushimaAlg.acn
```

## Internals

Not part of the public API — documented because the GPU extensions depend on them, and
they may change without a breaking release.

```@docs
JacobiElliptic.AlgorithmCall
JacobiElliptic.FukushimaAlg._F_E
```
