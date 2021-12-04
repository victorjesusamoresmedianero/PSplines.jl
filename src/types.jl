"""
    BSpline{T}
bspline is a structures which includes the knots and the 
vertices of a BSpline.
"""
struct BSpline{T}
    knots::Vector{T}
    vertices::Vector{T}
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
    BSpline2dim{T}
bspline2dim is a structures which includes the knots in x and y direction 
and the vertices (in vector form) of a 2 dimensional BSpline.
"""
struct BSpline2dim{T}
    knotsx::Vector{T}
    knotsy::Vector{T}
    vertices::Vector{T}
end

knotsx(bsdim2::BSpline2dim) = bsdim2.knotsx
knotsy(bsdim2::BSpline2dim) = bsdim2.knotsy
vertices(bsdim2::BSpline2dim) = bsdim2.vertices

"""
    buildKnots(xin, xend, nvertices)
Returns the knots of a `BSpline` with beginning xin, end xend
with nvertices components. The knots returned by this function
are equispaced between xin and xend.

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
function buildKnots(xin, xend, nvertices)
    knotsIn = LinRange(xin, xend, nvertices-2)
    knots = vcat(knotsIn[1]-(knotsIn[2]-knotsIn[1]), knotsIn,  knotsIn[end]+(knotsIn[2]-knotsIn[1]))
    return knots
end