using PSplines
using Test

@testset "PSplines.jl" begin
    include("test_MatricesBuilding.jl")
    include("test_types.jl")
    include("test_utils.jl")
end

