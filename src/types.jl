"""
    BSpline{T}
bspline is a structures which includes the knots and the 
vertices of a BSpline.
"""
struct BSpline{T}
    knots::Array{T}
    vertices::Array{T}
end

"""
    knots(bs::BSpline)
Returns the field knots for a certain `BSpline`.
"""
knots(bs::BSpline) = bs.knots

"""
    vertices(bs::BSpline)
Returns the field vertices for a certain `BSpline`.
"""
vertices(bs::BSpline) = bs.vertices

"""
    buildKnots(xin, xfin, nvertices)
Returns the knots of a `BSpline` with beginning xin, end xfin
with nvertices components. The knots returned by this function
are equispaced between xin and xfin.

#Examples
```julia julia-repl
julia> buildKnots(1.,5.,10)
10-element Vector{Float64}:
 0.4285714285714286
 1.0
 1.5714285714285714
 2.142857142857143
 2.7142857142857144
 3.2857142857142856
 3.857142857142857
 4.428571428571429
 5.0
 5.571428571428571
```
"""
function buildKnots(xin, xfin, nvertices)
    knotsIn = LinRange(xin, xfin, nvertices-2)
    knots = vcat(knotsIn[1]-(knotsIn[2]-knotsIn[1]), knotsIn,  knotsIn[end]+(knotsIn[2]-knotsIn[1]))
    return knots
end