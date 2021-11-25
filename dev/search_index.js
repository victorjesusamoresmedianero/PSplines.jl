var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PSplines","category":"page"},{"location":"#PSplines","page":"Home","title":"PSplines","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for PSplines.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [PSplines]","category":"page"},{"location":"#PSplines.BSpline","page":"Home","title":"PSplines.BSpline","text":"BSpline{T}\n\nbspline is a structures which includes the knots and the  vertices of a BSpline.\n\n\n\n\n\n","category":"type"},{"location":"#PSplines.D1matrix-Tuple{Int64}","page":"Home","title":"PSplines.D1matrix","text":"D1matrix(nvertices)\n\nReturns the first difference of dimension nvercies-1 x nvertices. It is important to note that the numerical first derivative of a  vector produces a vector with 1 component less than the input vector.\n\n#Examples\n\njulia> D1matrix(3)\n2×3 Matrix{Float64}:\n -1.0   1.0  0.0\n  0.0  -1.0  1.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.D2matrix-Tuple{Int64}","page":"Home","title":"PSplines.D2matrix","text":"D2matrix(nvertices)\n\nReturns the second difference of dimension nvercies-2 x nvertices. It is important to note that the numerical second derivative of a  vector produces a vector with 2 component less than the input vector.\n\n#Examples\n\njulia> D2matrix(4)\n2×4 Matrix{Float64}:\n 1.0  -2.0   1.0  0.0\n 0.0   1.0  -2.0  1.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.D3matrix-Tuple{Int64}","page":"Home","title":"PSplines.D3matrix","text":"D3matrix(nvertices)\n\nReturns the second difference of dimension nvercies-3 x nvertices. It is important to note that the numerical second derivative of a  vector produces a vector with 2 component less than the input vector.\n\n#Examples\n\njulia> D3matrix(5)\n2×5 Matrix{Float64}:\n -1.0   3.0  -3.0   1.0  0.0\n  0.0  -1.0   3.0  -3.0  1.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.buildKnots-Tuple{Any, Any, Any}","page":"Home","title":"PSplines.buildKnots","text":"buildKnots(xin, xfin, nvertices)\n\nReturns the knots of a BSpline with beginning xin, end xfin with nvertices components. The knots returned by this function are equispaced between xin and xfin.\n\n#Examples\n\njulia> buildKnots(1.,5.,10)\n10-element Vector{Float64}:\n 0.4285714285714286\n 1.0\n 1.5714285714285714\n 2.142857142857143\n 2.7142857142857144\n 3.2857142857142856\n 3.857142857142857\n 4.428571428571429\n 5.0\n 5.571428571428571\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.buildOm1-Tuple{Any}","page":"Home","title":"PSplines.buildOm1","text":"buildOm1(nvertices)\n\nReturns a matrix full of zeros with the appropiate dimensions to apply a selective penalization to the first differences vector. The default 0 value, implies no penalization on the first derivative.\n\n#Examples\n\njulia> buildOm1(5)\n4×4 Matrix{Float64}:\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.buildOm2-Tuple{Any}","page":"Home","title":"PSplines.buildOm2","text":"buildOm2(nvertices)\n\nReturns a matrix full of zeros with the appropiate dimensions to apply a selective penalization to the second differences vector. The default 0 value, implies no penalization on the second derivative.\n\n#Examples\n\njulia> buildOm2(5)\n3×3 Matrix{Float64}:\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n 0.0  0.0  0.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.buildOm3-Tuple{Any}","page":"Home","title":"PSplines.buildOm3","text":"buildOm3(nvertices)\n\nReturns a matrix full of zeros with the appropiate dimensions to apply a selective penalization to the third differences vector. The default 0 value, implies no penalization on the third derivative.\n\n#Examples\n\njulia> buildOm3(5)\n2×2 Matrix{Float64}:\n 0.0  0.0\n 0.0  0.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.buildW-Tuple{Any}","page":"Home","title":"PSplines.buildW","text":"buildW(nvertices)\n\nGiven the total number of experimental points, build the weight matrix for that number of points. Note that by default all the weights are equal to 1. #Examples\n\njulia> buildOm3(5)\n2×2 Matrix{Float64}:\n 0.0  0.0\n 0.0  0.0\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.equiSpacedPosition-Tuple{Float64, Vector{Float64}}","page":"Home","title":"PSplines.equiSpacedPosition","text":"equiSpacedPosition(x, knotsIn)\n\nReturns the bin in which x is contained for the bins defined by the  limits in knotsIn. By default if the value is in the limit it is pPosition at the right bin except for the las limit which is positioned to the left. #Examples\n\njulia> equiSpacedPosition(2, [1,2,3,4,5,6,7,8])\n2\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.evalBSpline-Tuple{Any, Any, Any}","page":"Home","title":"PSplines.evalBSpline","text":"evalBSpline(x, knots, vertices) Evaluation the cubic BSpline at x for vertices values vertices and the BSpline basis defined by evalBSplineBasis(x, knots). #Examples\n\njulia> evalBSpline(9.,buildKnots(1.,10., 12), ones(12))\n0.9999999999999999\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.evalBSplineBasis-Tuple{Float64, Vector{Float64}}","page":"Home","title":"PSplines.evalBSplineBasis","text":"evalBSplineBasis(x, knots)\n\nBuilds the evaluation of the cubic BSpline basis corresponding to the knots at x. #Examples\n\njulia> evalBSplineBasis(10.,buildKnots(1.,10., 12))\n12-element Vector{Float64}:\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.0\n 0.16666666666666663\n 0.6666666666666666\n 0.16666666666666666\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.knots-Tuple{BSpline}","page":"Home","title":"PSplines.knots","text":"knots(bs::BSpline)\n\nReturns the field knots for a certain BSpline.\n\n\n\n\n\n","category":"method"},{"location":"#PSplines.vertices-Tuple{BSpline}","page":"Home","title":"PSplines.vertices","text":"vertices(bs::BSpline)\n\nReturns the field vertices for a certain BSpline.\n\n\n\n\n\n","category":"method"}]
}
