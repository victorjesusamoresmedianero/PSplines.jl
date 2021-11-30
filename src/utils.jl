"""
    equiSpacedPosition(x, knotsIn)
Returns the bin in which x is contained for the bins defined by the 
limits in knotsIn. By default if the value is in the limit it is pPosition
at the right bin except for the las limit which is positioned to the left.
#Examples
```julia julia-repl
julia> equiSpacedPosition(2, [1,2,3,4,5,6,7,8])
2
```
"""
function equiSpacedPosition(x::Float64, knotsIn::Vector{Float64})

    @assert all(diff(knotsIn).>=0)    # bin limits are monotonically increaseing

    nbins = length(knotsIn) - 1
    lsubd = knotsIn[2] - knotsIn[1]

    proximityTol = 1e-4
    xinDist = abs(knotsIn[1]-x)
    xfinDist = abs(x - knotsIn[end])

    if x <= prevfloat(knotsIn[1])

        xinDist/lsubd < proximityTol ? s = 1 : throw(DomainError(x, "The value is out of the domain (left)")) 

    elseif x >= nextfloat(knotsIn[end])

        xfinDist/lsubd < proximityTol ? s = nbins : throw(DomainError(x, "The value is out of the domain (right)"))

    else
        s = floor(Int,(x - knotsIn[1])/lsubd) + 1

        s > nbins ? s = nbins : nothing
    end
    return s
end

equiSpacedPosition(x, knotsIn) = equiSpacedPosition(Float64(x), Float64.(knotsIn))



"""
    evalBSplineBasis(x, knots)
Builds the evaluation of the cubic BSpline basis corresponding to
the `knots` at `x`.
#Examples
```julia julia-repl
julia> evalBSplineBasis(10.,buildKnots(1.,10., 12))
12-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.16666666666666663
 0.6666666666666666
 0.16666666666666666
```  
"""
function evalBSplineBasis(x::Float64, knots::Vector{Float64})
    transfMatrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                  3. -6. 3. 0.;
                                 -3. 0. 3. 0.;
                                  1. 4. 1. 0.]
    nvertices = length(knots)
    lsubd = knots[2] - knots[1]
    knotsIn = knots[2:end-1] # the internal vertices are all except extrema

    bsplineBasis = zeros(Float64, nvertices)

    s = equiSpacedPosition(x, knotsIn)

    ξ = (x - knotsIn[s])/lsubd
    ξvlocal = [ξ^3, ξ^2, ξ, 1.0]
    
    bsplineBasis[s:s+3] = ξvlocal'*transfMatrix

    return bsplineBasis::Vector{Float64}
end

evalBSplineBasis(x, knots) = evalBSplineBasis(Float64(x), Float64.(knots)) 


"""
evalBSpline(x, knots, vertices)
Evaluation the cubic BSpline at `x` for vertices values `vertices`
and the BSpline basis defined by `evalBSplineBasis(x, knots)`.
#Examples
```julia julia-repl
julia> evalBSpline(9.,buildKnots(1.,10., 12), ones(12))
0.9999999999999999
```  
"""
function evalBSpline(x, knots, vertices)
    return evalBSplineBasis(x, knots)'*vertices
end

evalBSpline(x, bs::BSpline) = evalBSpline(x, bs.knots, bs.vertices) 


function expandBSpline(knots, vertices, xin, xend)
    knotsInmin = knots[2]
    knotsInmax = knots[end-1]

    lsubd = knots[2] - knots[1]

    intervLeft = (knotsInmin - xin)/lsbud
    intervRight = (xend - knotsInmax)/lsubd

    intervLeft > 0 ? nKnotsLeft = ceil(intervLeft) : nKnotsLeft = 0
    intervRight > 0 ? nKnotsRight = ceil(intervRight) : nKnotsRight = 0

    addKnotsLeft = LinRange(knots[1]-lsubd*nKnotsLeft, knots[1]-lsubd, nKnotsLeft)
    addKnotsRight = LinRange(knots[end]+lsubd, knots[end]+lsubd*nKnotsRight, nKnotsRight)

    knotsNew = vcat(addKnotsLeft, knots, addKnotsRight)
end