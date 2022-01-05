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
function equiSpacedPosition(x, knotsIn)

    @assert all(diff(knotsIn).>=0)    # bin limits are monotonically increaseing

    nbins = length(knotsIn) - 1
    lsubd = knotsIn[2] - knotsIn[1]

    proximityTol = 1e-4
    xinDist = abs(knotsIn[1]-x)
    xfinDist = abs(x - knotsIn[end])

    if x <= prevfloat(Float64(knotsIn[1]))

        xinDist/lsubd < proximityTol ? s = 1 : throw(DomainError(x, "The value is out of the domain (left)")) 

    elseif x >= nextfloat(Float64(knotsIn[end]))

        xfinDist/lsubd < proximityTol ? s = nbins : throw(DomainError(x, "The value is out of the domain (right)"))

    else
        s = floor(Int,(x - knotsIn[1])/lsubd) + 1

        s > nbins ? s = nbins : nothing
    end
    return s
end
# equiSpacedPosition(x, knotsIn) = equiSpacedPosition(Float64(x), Float64.(knotsIn))

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
function evalBSplineBasis(x, knots) 
    transfMatrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                  3. -6. 3. 0.;
                                 -3. 0. 3. 0.;
                                  1. 4. 1. 0.]
    nvertices = length(knots)
    lsubd = knots[2] - knots[1]
    knotsIn = knots[2:end-1] # the internal vertices are all except extrema

    bsplineBasis = zeros(typeof(x), nvertices)

    s = equiSpacedPosition(x, knotsIn)

    ξ = (x - knotsIn[s])/lsubd
    ξvlocal = [ξ^3, ξ^2, ξ, 1.0]
    
    bsplineBasis[s:s+3] = ξvlocal'*transfMatrix

    return bsplineBasis
end

#evalBSplineBasis(x, knots) = evalBSplineBasis(Float64(x), Float64.(knots))

## Document
function evalBSplineBasis2dim(x, y, knotsx, knotsy)

    bsplineBasisx = evalBSplineBasis(x, knotsx)
    bsplineBasisy = evalBSplineBasis(y, knotsy)
    return vec(bsplineBasisx*bsplineBasisy')
end

#evalBSplineBasis2dim(x, y, knotsx, knotsy) = evalBSplineBasis2dim(Float64(x), Float64(y), Float64.(knotsx), Float64.(knotsy))

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
function evalBSplineBasisD1(x, knots)

    transfMatrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                   3. -6. 3. 0.;
                                   -3. 0. 3. 0.;
                                    1. 4. 1. 0.]
    nvertices = length(knots)
    lsubd = knots[2] - knots[1]
    knotsIn = knots[2:end-1] # the internal vertices are all except extrema
    bsplineBasisD1 = zeros(typeof(x), nvertices)

    s = equiSpacedPosition(x, knotsIn)

    ξ = (x - knotsIn[s])/lsubd
    ξvlocalD1 = [3. * ξ^2, 2. * ξ, 1., 0.]


    bsplineBasisD1[s:s+3] = ξvlocalD1'*transfMatrix*(1/lsubd)
    return bsplineBasisD1
end
#evalBSplineBasisD1(x, knots) = evalBSplineBasisD1(Float64(x), Float64.(knots))


function evalBSplineBasis2dimD1x(x, y, knotsx, knotsy)
    dbsplineBasisx = evalBSplineBasisD1(x, knotsx)
    bsplineBasisy = evalBSplineBasis(y, knotsy)
    return vec(dbsplineBasisx*bsplineBasisy')
end

function evalBSplineBasis2dimD1y(x, y, knotsx, knotsy)
    bsplineBasisx = evalBSplineBasis(x, knotsx)
    dbsplineBasisy = evalBSplineBasisD1(y, knotsy)
    return vec(bsplineBasisx*dbsplineBasisy')
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

## Document
function evalBSpline2dim(x, y, knotsx, knotsy, vertices)
    return evalBSplineBasis2dim(x, y, knotsx, knotsy)'*vertices
end

evalBSpline2dim(x, y, bs2dim::BSpline2dim) = evalBSpline2dim(x, y, bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices)


function evalBSpline2dimD1x(x, y, knotsx, knotsy, vertices)
    return evalBSplineBasis2dimD1x(x, y, knotsx, knotsy)'*vertices
end

evalBSpline2dimD1x(x, y, bs2dim::BSpline2dim) = evalBSpline2dimD1x(x, y, bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices)


function evalBSpline2dimD1y(x, y, knotsx, knotsy, vertices)
    return evalBSplineBasis2dimD1y(x, y, knotsx, knotsy)'*vertices
end

evalBSpline2dimD1y(x, y, bs2dim::BSpline2dim) = evalBSpline2dimD1y(x, y, bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices)

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




# Document
function nVerticesNeededLeftRight(knots, xin, xend)
    knotsInmin = knots[2]
    knotsInmax = knots[end-1]

    lsubd = knots[2] - knots[1]

    intervLeft = (knotsInmin - xin)/lsubd
    intervRight = (xend - knotsInmax)/lsubd

    intervLeft > 0 ? nKnotsLeft = Int64(ceil(intervLeft)) : nKnotsLeft = 0
    intervRight > 0 ? nKnotsRight = Int64(ceil(intervRight)) : nKnotsRight = 0
    return nKnotsLeft, nKnotsRight
end


# Document
function extendedKnots(knots, nKnotsLeft, nKnotsRight)

    lsubd = knots[2] - knots[1]

    addKnotsLeft = range(knots[1]-lsubd*nKnotsLeft, knots[1]-lsubd, length = nKnotsLeft)
    addKnotsRight = range(knots[end]+lsubd, knots[end]+lsubd*nKnotsRight, length = nKnotsRight)

    knotsNew = vcat(addKnotsLeft, knots, addKnotsRight)

    return knotsNew
end

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
function expandBSpline(knots, vertices, 
                       xin, xend)
    
    nvertices = length(knots)

    nKnotsLeft, nKnotsRight = nVerticesNeededLeftRight(knots, xin, xend)
    knotsNew = extendedKnots(knots, nKnotsLeft, nKnotsRight)
    nverticesNew = length(knotsNew)

    A = zeros(nverticesNew, nverticesNew)

    D2 = D2matrix(nverticesNew)
    Ω2 = Ω2matrix(nverticesNew)
    W = Wmatrix(nverticesNew)

    A[nKnotsLeft+1 : (nKnotsLeft + nvertices), nKnotsLeft+1 : (nKnotsLeft + nvertices)] =Float64.(I(nvertices)) 
    
    for i = nKnotsLeft+1:(nKnotsLeft + nvertices)
        W[i,i] = 1e4
    end

    Ω2[diagind(Ω2)] .= 1. 

    Asys = A'*W*A + D2'*Ω2*D2
    bsys = A'*W*vcat(zeros(nKnotsLeft), vertices, zeros(nKnotsRight))

    verticesNew = Asys\bsys

    return BSpline(knotsNew, verticesNew)
end

# TO BE CHECKED AND DOCUMENTED
function expandBSpline2dim(knotsx, knotsy, vertices, 
                           xin, xend, yin, yend)
    nverticesx = length(knotsx)
    nverticesy = length(knotsy)

    nKnotsLeftx, nKnotsRightx = nVerticesNeededLeftRight(knotsx, xin, xend)
    knotsNewx = extendedKnots(knotsx, nKnotsLeftx, nKnotsRightx)

    nverticesxNew = length(knotsNewx)

    nKnotsLefty, nKnotsRighty = nVerticesNeededLeftRight(knotsy, yin, yend)
    knotsNewy = extendedKnots(knotsy, nKnotsLefty, nKnotsRighty)

    nverticesyNew = length(knotsNewy)

    indicesVector = reshape(1:nverticesxNew*nverticesyNew, nverticesxNew, nverticesyNew)

    yindVerticesNew = zeros(nverticesxNew*nverticesyNew)
    NNew = zeros(nverticesxNew*nverticesyNew, nverticesxNew*nverticesyNew)

    # Identify a priori known vertices

    verticesMatrix = reshape(vertices, nverticesx, nverticesy)

    wi = 1e4

    for i = (nKnotsLeftx+1):(nverticesx+nKnotsLeftx)
        for j = (nKnotsLefty+1):(nvertivesy+nKnotsLefty)
            yindVerticesNew[indicesVector[i,j]] = verticesMatrix[i,j]
            NNew[indicesVector[i,j],indicesVector[i,j]] = wi
        end
    end

    Px, Py = PD2matrices2dim(nverticesxNew, nverticesyNew)
    λx = 1.
    λy = 1.

    Asys = NNew'*NNew + λx*Px + λy*Py
    bsys = NNew'*yindVerticesNew

    verticesNew =  Asys\bsys    

    return BSpline2dim(knotsNewx, knotsNewy, verticesNew)
end

#expandBSpline2dim(knotsx, knotsy, vertices, xin, xend, yin, yend) = expandBSpline2dim(Float64.(knotsx), Float64.(knotsy), Float64.(vertices), Float64(xin), Float64(xend), Float64(yin), Float64(yend))

expandBSpline2dim(bs2dim::BSpline2dim, xin, xend, yin, yend) = expandBSpline2dim(bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices, xin, xend, yin, yend)

## Document
#expandBSpline(knots, vertices, xin, xend) = expandBSpline(Float64.(knots), Float64.(vertices), 
                                                          #Float64(xin), Float64(xend))
expandBSpline(bs::BSpline, xin, xend) = expandBSpline(bs.knots, bs.vertices, xin, xend)

## Document
function buildSystemMatrix(N;
                            W = Wmatrix(size(N,1)), 
                           Ω1 = Ω1matrix(size(N,2)), 
                           Ω2 = Ω2matrix(size(N,2)), 
                           Ω3 = Ω3matrix(size(N,2)))
    nvertices = size(N,2)
    D1 = D1matrix(nvertices)
    D2 = D2matrix(nvertices)
    D3 = D3matrix(nvertices)

    Asys = N'*W*N + D1'*Ω1*D1 + D2'*Ω2*D2 + D3'*Ω3*D3
    
    return Asys
end
## Document
function buildSystemIndVect(N, y;
                            W = Wmatrix(size(N,1)))
    bsys = N'*W*y
    return bsys
end
