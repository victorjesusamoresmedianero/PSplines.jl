# Testing the function D1matrix. The first derivative of a constant is zero.
@testset "D1matrix building" begin
    @test all(D1matrix(5)*ones(5).==0)
end

@testset "D2matrix building" begin
    @test all(D2matrix(4).==D1matrix(3)*D1matrix(4))
end

@testset "D3matrix building" begin
    @test all(D3matrix(5).==D2matrix(4)*D1matrix(5))
end