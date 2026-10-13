using Test
using JacobiElliptic
using ArbNumerics
using DelimitedFiles: readdlm
using Enzyme, EnzymeTestUtils
using Reactant

setprecision(ArbNumerics.ArbFloat, 200)

include("default_algorithm_tests.jl")
include("default_tests.jl")
include("edge_case_tests.jl")
include("identity_tests.jl")
include("carlson_tests.jl")
include("jacobi_tests.jl")
include("autodiff.jl")
include("metal_tests.jl")
include("reactant_tests.jl")
