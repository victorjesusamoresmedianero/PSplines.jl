module PSplines

using LinearAlgebra: diagind

# buildmatrices.jl
export D1matrix, D2matrix, D3matrix, Ω1matrix, Ω2matrix, Ω3matrix, Wmatrix
include("buildMatrices.jl")

# types.jl
export BSpline, knots, vertices, buildKnots
include("types.jl")

#utils.jl
export equiSpacedPosition, evalBSplineBasis, evalBSpline
include("utils.jl")
end
