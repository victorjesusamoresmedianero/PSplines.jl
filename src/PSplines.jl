module PSplines

# buildmatrices.jl
export D1matrix, D2matrix, D3matrix, buildOm1, buildOm2, buildOm3, buildW
include("buildMatrices.jl")

# types.jl
export BSpline, knots, vertices, buildKnots
include("types.jl")

#utils.jl
export equiSpacedPosition, evalBSplineBasis
include("utils.jl")
end
