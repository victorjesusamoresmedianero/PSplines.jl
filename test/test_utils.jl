@testset "at the beggining of equiSpacedPosition" begin
    @test equiSpacedPosition(1, [1,2,3,4,5,6,7,8]) == 1
end

@testset "inside of equiSpacedPosition 1" begin
    @test equiSpacedPosition(3.5, [1,2,3,4,5,6,7,8]) == 3
end

@testset "inside of equiSpacedPosition 2" begin
    @test equiSpacedPosition(2, [1,2,3,4,5,6,7,8]) == 2
end

@testset "at the end of equiSpacedPosition" begin
    @test equiSpacedPosition(8, [1,2,3,4,5,6,7,8]) == 7
end


@testset "the beginning of evalBSplineBasis" begin
    @test all(evalBSplineBasis(1.,buildKnots(1.,10., 12)) .≈ vcat([1.0/6.0, 4.0/6.0, 1.0/6.0, 0.], zeros(Float64,8)))
end


@testset "the end of evalBSplineBasis" begin
    @test all(evalBSplineBasis(10.,buildKnots(1.,10., 12)) .≈ vcat(zeros(Float64, 8), [0., 1.0/6.0, 4.0/6.0, 1.0/6.0]))
end