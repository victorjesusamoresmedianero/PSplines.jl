

function Base.promote_rule(::Type{ForwardDiff.Dual{T}}, ::Type{ForwardDiff.Dual{R}}) where {T,R} 
    ForwardDiff.Dual{promote_type(T,R)}
  end
  
function Base.promote_rule(::Type{ForwardDiff.Dual{T}}, ::Type{R}) where {T,R}
    ForwardDiff.Dual{promote_type(T,R)}
end

function Base.convert(::Type{ForwardDiff.Dual{T}}, x::ForwardDiff.Dual) where T
    ForwardDiff.Dual(convert(T, x.x), convert(T, x.ϵ))
  end
  
function Base.convert(::Type{ForwardDiff.Dual{T}}, x::Real) where T
    ForwardDiff.Dual(convert(T, x), zero(T))
end
