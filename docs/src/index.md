```@meta
CurrentModule = JacobiElliptic
```

# JacobiElliptic
JacobiElliptic is an implementation of Toshio Fukushima's & Billie C. Carlson's algorithms for calculating [Elliptic Integrals and Jacobi Elliptic Functions](https://ieeexplore.ieee.org/document/7203795). 

### Features
  - Type stable and preserving
  - [Metal.jl](https://github.com/JuliaGPU/Metal.jl) and [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) compatible
  - Automatic Differentiable with ForwardDiff, Zygote and Enzyme