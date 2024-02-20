# JacobiElliptic
JacobiElliptic is an implementation of Toshio Fukushima's algorithms for calculating [Elliptic Integrals and Jacobi Elliptic Functions](https://ieeexplore.ieee.org/document/7203795). 

## Features
  - Type stable and preserving
  - [Metal.jl](https://github.com/JuliaGPU/Metal.jl) and [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) compatible
## Documentation
[![Dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://dchang10.github.io/JacobiElliptic.jl/dev/)

## Repo Status
[![Build Status](https://github.com/dchang10/JacobiElliptic.jl/workflows/CI/badge.svg)](https://github.com/dchang10/JacobiElliptic.jl/actions)
[![Coverage](https://codecov.io/gh/dchang10/JacobiElliptic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dchang10/JacobiElliptic.jl)

## Incomplete Elliptic Integrals
|Function | Definition |
| --- | --- |
| `F(φ, m)` | $F(\phi\|m)$: Incomplete elliptic integral of the first kind|
| `E(φ, m)` |  $E(\phi\|m)$: Incomplete elliptic integral of the second kind |
| `Pi(n, φ, m)` | $\Pi(n;\\,\phi\|\\, m)$: Incomplete elliptic integral of the third kind|
| `J(n, φ, m)` | $J (n;\\, \phi \|\\,m)=\frac{\Pi(n;\\,\phi\|\\, m) - F(\phi,m)}{n}$: Associated incomplete elliptic integral of the third kind|

## Complete Elliptic Integrals
|Function | Definition |
| --- | --- |
| `K(m)` | $K(m)$: Complete elliptic integral of the first kind|
| `E(m)` |  $E(m)$: Complete elliptic integral of the second kind |
| `Pi(n, m)` | $\Pi(n\|\\, m)$: Complete elliptic integral of the third kind|
| `J(n, m)` | $J (n\|\\,m)=\frac{\Pi(n\|\\, m) - K(m)}{n}$: Associated incomplete elliptic integral of the third kind|

## Jacobi Elliptic Functions
|Function | Definition |
| --- | --- |
| `sn(u, m)` | $\text{sn}(u \| m) = \sin(\text{am}(u \| m))$ |
| `cn(u, m)` | $\text{cn}(u \| m) = \cos(\text{am}(u \| m))$ |
| `asn(u, m)` | $\text{asn}(u \| m) = \text{sn}^{-1}(u \| m)$ |
| `acn(u, m)` | $\text{cn}(u \| m) = \text{cn}^{-1}(u \| m)$ |
