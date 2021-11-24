"""
    D1matrix(nvertices)
Returns the first difference of dimension `nvercies-1 x nvertices`.
It is important to note that the numerical first derivative of a 
vector produces a vector with 1 component less than the input vector.

#Examples
```julia julia-repl
julia> D1matrix(3)
2×3 Matrix{Float64}:
 -1.0   1.0  0.0
  0.0  -1.0  1.0
```
"""
function D1matrix(nvertices::Int)
    D1matrix = zeros(Float64,nvertices-1,nvertices)
    for i=1:nvertices-1
        D1matrix[i,i] = -1.
        D1matrix[i,i+1] = 1.
    end
    return D1matrix::Array{Float64,2}
end

"""
    D2matrix(nvertices)
Returns the second difference of dimension `nvercies-2 x nvertices`.
It is important to note that the numerical second derivative of a 
vector produces a vector with 2 component less than the input vector.

#Examples
```julia julia-repl
julia> D2matrix(4)
2×4 Matrix{Float64}:
 1.0  -2.0   1.0  0.0
 0.0   1.0  -2.0  1.0
```
"""
function D2matrix(nvertices::Int)

    D2matrix = zeros(Float64,nvertices-2,nvertices)

    for i=1:nvertices-2
        D2matrix[i,i] = 1.
        D2matrix[i,i+1] = -2.
        D2matrix[i,i+2] = 1.
    end

    return D2matrix::Array{Float64,2}
end


"""
    D3matrix(nvertices)
Returns the second difference of dimension `nvercies-3 x nvertices`.
It is important to note that the numerical second derivative of a 
vector produces a vector with 2 component less than the input vector.

#Examples
```julia julia-repl
julia> D3matrix(5)
2×5 Matrix{Float64}:
 -1.0   3.0  -3.0   1.0  0.0
  0.0  -1.0   3.0  -3.0  1.0
```
"""
function D3matrix(nvertices::Int)

    D3matrix = zeros(Float64,nvertices-3,nvertices)

    for i=1:nvertices-3
        D3matrix[i,i] = -1.
        D3matrix[i,i+1] = 3.
        D3matrix[i,i+2] = -3.
        D3matrix[i,i+3] = 1.
    end

    return D3matrix::Array{Float64,2}
end


"""
    buildOm1(nvertices)
Returns a matrix full of zeros with the appropiate dimensions to apply
a selective penalization to the first differences vector. The default
0 value, implies no penalization on the first derivative.

#Examples
```julia julia-repl
julia> buildOm1(5)
4×4 Matrix{Float64}:
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
```
"""
function buildOm1(nvertices)
    return zeros(Float64, nvertices - 1, nvertices - 1)
end

"""
    buildOm2(nvertices)
Returns a matrix full of zeros with the appropiate dimensions to apply
a selective penalization to the second differences vector. The default
0 value, implies no penalization on the second derivative.
#Examples
```julia julia-repl
julia> buildOm2(5)
3×3 Matrix{Float64}:
 0.0  0.0  0.0
 0.0  0.0  0.0
 0.0  0.0  0.0
```
"""
function buildOm2(nvertices)
    return zeros(Float64, nvertices - 2, nvertices - 2)
end

"""
    buildOm3(nvertices)
Returns a matrix full of zeros with the appropiate dimensions to apply
a selective penalization to the third differences vector. The default
0 value, implies no penalization on the third derivative.
#Examples
```julia julia-repl
julia> buildOm3(5)
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0
```
"""
function buildOm3(nvertices)
    return zeros(Float64, nvertices - 3, nvertices - 3)
end