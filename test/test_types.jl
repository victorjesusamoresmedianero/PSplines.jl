@testset "testing length" begin
    @test length(buildKnots(1.0, 4., 10)) == 10
end

@testset "testing bounds" begin
    @test buildKnots(1.0, 4., 10)[2] == 1.0 && buildKnots(1.0, 4., 10)[end-1] == 4.0 
end