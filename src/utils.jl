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

    if x <= knotsIn[1]

        xinDist/lsubd < proximityTol ? s = 1 : throw(DomainError(x, "The value is out of the domain (left)")) 

    elseif x >= knotsIn[end]

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
function evalBSplineBasis(x, knots)
    transfMatrix = (1. / 6.) * [ -1. 3. -3. 1.;
                                  3. -6. 3. 0.;
                                 -3. 0. 3. 0.;
                                  1. 4. 1. 0.]
    nvertices = length(knots)
    lsubd = knots[2] - knots[1]
    knotsIn = knots[2:end-1] # the internal vertices are all except extrema

    bsplineBasis = zeros(promote_type(typeof(x), eltype(knots)), nvertices)

    s = equiSpacedPosition(x, knotsIn)

    ξ = (x - knotsIn[s])/lsubd
    ξvlocal = [ξ^3, ξ^2, ξ, one(promote_type(typeof(x), eltype(knots)))]
    
    bsplineBasis[s:s+3] = ξvlocal'*transfMatrix

    return bsplineBasis
end


"""
    evalBSplineBasis2dim(x, y, knotsx, knotsy)
Builds the evaluation of the dyadic product of 2 cubic BSplines 
basis, the first corresponding to `knotsx` at `x` and the second 
corresponding to `knotsy`at `y`. The output is already vectorized
instead of being a matrix.
#Examples
```julia julia-repl
julia> evalBSplineBasis2dim(3.,2.,buildKnots(1.,10., 12),buildKnots(0.,5.,14))
168-element Vector{Float64}:
 0.0
 0.0
 0.0
 0.0
 ⋮
 0.0
 0.0
 0.0
```  
"""
function evalBSplineBasis2dim(x, y, knotsx, knotsy)

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

"""
    evalBSpline2dim(x, y, knotsx, knotsy, vertices)
Evaluation of the cubic bidimensional BSpline with knots `knotsx`, `knotsy`
at `(x,y)` for vertices `vertices`.

#Examples
```julia julia-repl
julia> knotsx = knotsy = 0:0.1:10;
julia> vertices = rand(length(knotsx)*length(knotsy ));
julia> evalBSpline2dim(2., 3., knotsx, knotsy, vertices)
0.3800038047730809
```  
"""
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
    evalBSpline2dimD1x(x, y, knotsx, knotsy, vertices)
Evaluation of the derivative with respect to `x`of the 
cubic bidimensional BSpline with knots `knotsx`, `knotsy`
at `(x,y)` for vertices `vertices`.

#Examples
```julia julia-repl
julia> knotsx = knotsy = 0:0.1:10;
julia> vertices = rand(length(knotsx)*length(knotsy ));
julia> evalBSpline2dimD1x(2., 3., knotsx, knotsy, vertices)
-2.0863387254122316
```  
"""
function evalBSpline2dimD1x(x, y, knotsx, knotsy, vertices)
    return evalBSplineBasis2dimD1x(x, y, knotsx, knotsy)'*vertices
end

evalBSpline2dimD1x(x, y, bs2dim::BSpline2dim) = evalBSpline2dimD1x(x, y, bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices)

"""
    evalBSpline2dimD1y(x, y, knotsx, knotsy, vertices)
Evaluation of the derivative with respect to `y`of the 
cubic bidimensional BSpline with knots `knotsx`, `knotsy`
at `(x,y)` for vertices `vertices`.

#Examples
```julia julia-repl
julia> knotsx = knotsy = 0:0.1:10;
julia> vertices = rand(length(knotsx)*length(knotsy ));
julia> evalBSpline2dimD1y(2., 3., knotsx, knotsy, vertices)
2.345106484778582
```  
"""
function evalBSpline2dimD1y(x, y, knotsx, knotsy, vertices)
    return evalBSplineBasis2dimD1y(x, y, knotsx, knotsy)'*vertices
end

evalBSpline2dimD1y(x, y, bs2dim::BSpline2dim) = evalBSpline2dimD1y(x, y, bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices)



"""
    nVerticesNeededLeftRight(knots, xin, xend)
Given the knots vector `knots` and the limits `xin` and `xend`,
determines the number of vertices that knots needs  on the left 
and on the right needed to reach the values `xin` and `xend`.

#Examples
```julia julia-repl
julia> knots = 0.:1.:10.;
julia> nVerticesNeededLeftRight(knots, -2, 14)
(3, 5)
```  
"""
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


"""
    extendedKnots(knots, nKnotsLeft, nKnotsRight)
Given the knots vector `knots` and `nKnotsLeft` and `nKnotsRight`
extends the `knots` on the left `nKnotsLeft` and on the right `nKnotsRight`.

#Examples
```julia julia-repl
julia> knots = 0.:1.:10.;
julia> extendedKnots(knots, 2, 2)
15-element Vector{Float64}:
 -2.0
 -1.0
  ⋮
 11.0
 12.0
```  
"""
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

expandBSpline(bs::BSpline, xin, xend) = expandBSpline(bs.knots, bs.vertices, xin, xend)

"""
    expandBSpline2dim(knotsx, knotsy, vertices, xin, xend, yin, yend)
Given a 2D BSpline by `knotsx` , `knotsy` and `vetices`
compute the expanded 2D BSpline to the limits `xin, `xend`, `yin`, `yend`.

#Examples
```julia julia-repl
julia> knotsx = knotsy = 0:0.5:10;
julia> vertices = rand(length(knotsx)*length(knotsy ));
julia> expandBSpline2dim(knotsx, knotsy, vertices, -1., 11., -1., 11.)
BSpline2dim{Float64}([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0  …  7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5], 
[-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0  …  7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5], 
[-1.0029300060952653, -0.4790706367303448, 0.03153695697931029, 0.48370119176194715, 0.8342290770083228, 1.1030443114892081, 1.2623356760827873, 1.2585422976688359, 1.131589871095645, 0.9860138921208538  
…  0.8984551404840291, 0.7049426036058777, 0.4552726974231055, 0.2849091663041328, 0.22446630431980158, 0.22686278460500844, 0.3058408589841374, 0.4651356516733504, 0.6732304620674567, 0.8966217568026448])
```  
"""
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

    w = 1e10

    for i = (nKnotsLeftx+1):(nverticesx+nKnotsLeftx)
        for j = (nKnotsLefty+1):(nverticesy+nKnotsLefty)
            locali = i - nKnotsLeftx  ;localj = j - nKnotsLefty
            yindVerticesNew[indicesVector[i,j]] = verticesMatrix[locali,localj]
            NNew[indicesVector[i,j],indicesVector[i,j]] = 1
        end
    end

    Wmat = Diagonal(fill(w, size(NNew,1)))

    Px, Py = PD2matrices2dim(nverticesxNew, nverticesyNew)
    λx = 1.
    λy = 1.

    Asys = NNew'*Wmat*NNew + λx*Px + λy*Py
    bsys = NNew'*Wmat*yindVerticesNew

    verticesNew =  Asys\bsys    

    return BSpline2dim(knotsNewx, knotsNewy, verticesNew)
end


expandBSpline2dim(bs2dim::BSpline2dim, xin, xend, yin, yend) = expandBSpline2dim(bs2dim.knotsx, bs2dim.knotsy, bs2dim.vertices, xin, xend, yin, yend)

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
