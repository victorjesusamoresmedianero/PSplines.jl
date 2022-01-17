module PSplines

using LinearAlgebra: diagind, I
using Kronecker: kronecker

# buildmatrices.jl
export D1matrix, D2matrix, D3matrix, Ω1matrix, Ω2matrix, Ω3matrix, Wmatrix
export PD2matrices2dim
include("buildMatrices.jl")

# types.jl
export BSpline, BSpline2dim, knots, vertices, buildKnots
include("types.jl")

#utils.jl
export equiSpacedPosition, buildSystemMatrix, buildSystemIndVect
export evalBSplineBasis, evalBSpline
export evalBSplineBasis2dim, evalBSpline2dim
export evalBSplineBasisD1, evalBSplineD1
export evalBSplineBasis2dimD1x, evalBSpline2dimD1x, evalBSplineBasis2dimD1y, evalBSpline2dimD1y 
export expandBSpline, expandBSpline2dim
include("utils.jl")
end
