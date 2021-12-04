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
function equiSpacedPosition(x::T, knotsIn::Vector{T}) where {T<:Real}

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
function evalBSplineBasis(x::T, knots::Vector{T}) where {T<:Real}
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

    return bsplineBasis
end

function evalBSplineBasis2dim(x::T, y::T, knotsx::Vector{T}, knotsy::Vector{T}) where {T<:Real}

    bsplineBasisx = evalBSplineBasis(x, knotsx)
    bsplineBasisy = evalBSplineBasis(y, knotsy)
    return vec(bsplineBasisx*bsplineBasisy')
end
 
"""
    evalBSplineBasisD1(x, knots)
Builds the evaluation of the derivative of the cubic BSpline basis corresponding to
the `knots` at `x`.
#Examples
```julia julia-repl
julia> evalBSplineBasisD1(10.,buildKnots(1.,10., 12))
12-element Vector{Float64}:
  0.0
  0.0
  0.0
  0.0
  ⋮
 -0.5
  0.0
  0.5
```  
"""
function evalBSplineBasisD1(x::T, knots::Vector{T}) where {T<:Real}

    transfMatrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                   3. -6. 3. 0.;
                                   -3. 0. 3. 0.;
                                    1. 4. 1. 0.]
    nvertices = length(knots)
    lsubd = knots[2] - knots[1]
    knotsIn = knots[2:end-1] # the internal vertices are all except extrema
    bsplineBasisD1 = zeros(Float64, nvertices)

    s = equiSpacedPosition(x, knotsIn)

    ξ = (x - knotsIn[s])/lsubd
    ξvlocalD1 = [3. * ξ^2, 2. * ξ, 1., 0.]


    bsplineBasisD1[s:s+3] = ξvlocalD1'*transfMatrix*(1/lsubd)
    return bsplineBasisD1
end

"""
    evalBSpline(x, knots, vertices)
Evaluation of the cubic BSpline at `x` for vertices values `vertices`
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


function evalBSpline2dim(x, y, knotsx, knotsy, vertices)
    return evalBSplineBasis2dim(x, y, knotsx, knotsy)'*vertices
end

evalBSpline2dim(x, y, bs2dim::BSpline2dim) = evalBSpline2dim(x, y, bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices)
"""
    evalBSplineD1(x, knots, vertices)
Evaluating the derivative ofthe cubic BSpline at `x` for vertices 
values `vertices` and the BSpline basis derivative defined by 
`evalBSplineBasisD1(x, knots)`.
#Examples
```julia julia-repl
julia> evalBSplineD1(9.,buildKnots(1.,10., 12), ones(12))
0.0
```  
"""
function evalBSplineD1(x, knots, vertices)
    return evalBSplineBasisD1(x, knots)'*vertices
end

evalBSplineD1(x, bs::BSpline) = evalBSplineD1(x, bs.knots, bs.vertices)



"""
    expandBSpline(knots, vertices, xin, xend)
Given the knots and the vertices of a given BSpline with a certain range of validity,
it returns a new BSpline which extends the range of validity of the previous BSpline
betwenn `xin` and `xend`.
#Examples
```julia julia-repl
julia> bs2 = expandBSpline(bs1, 0., 5.);
julia> bs2.knots
20-element Vector{Float64}:
 -0.6666666666666663
 -0.33333333333333304
  2.220446049250313e-16
  ⋮
  5.333333333333333
  5.666666666666666
julia> bs2.vertices
20-element Vector{Float64}:
  -1.8150147193435557
  -1.2594147415453065
  -0.7038147637472744
   ⋮
  27.07391859179092
  29.851651947327774
```  
"""
function expandBSpline(knots::Vector{T}, vertices::Vector{T}, 
                       xin::T, xend::T) where {T<:Real}
    knotsInmin = knots[2]
    knotsInmax = knots[end-1]

    nvertices = length(knots)

    lsubd = knots[2] - knots[1]

    intervLeft = (knotsInmin - xin)/lsubd
    intervRight = (xend - knotsInmax)/lsubd

    intervLeft > 0 ? nKnotsLeft = Int64(ceil(intervLeft)) : nKnotsLeft = 0
    intervRight > 0 ? nKnotsRight = Int64(ceil(intervRight)) : nKnotsRight = 0

    addKnotsLeft = range(knots[1]-lsubd*nKnotsLeft, knots[1]-lsubd, length = nKnotsLeft)
    addKnotsRight = range(knots[end]+lsubd, knots[end]+lsubd*nKnotsRight, length = nKnotsRight)

    knotsNew = vcat(addKnotsLeft, knots, addKnotsRight)
    nverticesNew = length(knotsNew)
    A = zeros(nverticesNew, nverticesNew)

    D2 = D2matrix(nverticesNew)
    Ω2 = Ω2matrix(nverticesNew)
    W = Wmatrix(nverticesNew)

    A[nKnotsLeft+1 : (nKnotsLeft + nvertices), nKnotsLeft+1 : (nKnotsLeft + nvertices)] =Float64.(I(nvertices)) 
    
    for i = nKnotsLeft:(nKnotsLeft + nvertices)
        W[i,i] = 1e4
    end

    Ω2[diagind(Ω2)] .= 1. 

    Asys = A'*W*A + D2'*Ω2*D2
    bsys = A'*W*vcat(zeros(nKnotsLeft), vertices, zeros(nKnotsRight))

    verticesNew = Asys\bsys

    return BSpline(knotsNew, verticesNew)
end

expandBSpline(knots, vertices, xin, xend) = expandBSpline(Float64.(knots), Float64.(vertices), 
                                                          Float64(xin), Float64(xend))
expandBSpline(bs::BSpline, xin, xend) = expandBSpline(bs.knots, bs.vertices, xin, xend)


function buildSystemMatrix(N::Matrix{T};
                            W = Wmatrix(size(N,1))::Matrix{T}, 
                           Ω1 = Ω1matrix(size(N,2))::Matrix{T}, 
                           Ω2 = Ω2matrix(size(N,2))::Matrix{T}, 
                           Ω3 = Ω3matrix(size(N,2))::Matrix{T}) where {T<:Real}
    nvertices = size(N,2)
    D1 = D1matrix(nvertices)
    D2 = D2matrix(nvertices)
    D3 = D3matrix(nvertices)

    Asys = N'*W*N + D1'*Ω1*D1 + D2'*Ω2*D2 + D3'*Ω3*D3
    
    return Asys
end

function buildSystemIndVect(N::Matrix{T}, y::Vector{T};
                            W = Wmatrix(size(N,1))::Matrix{T}) where {T <: Real}
    bsys = N'*W*y
    return bsys
end
