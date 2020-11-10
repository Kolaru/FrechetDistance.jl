using FrechetDistance
using Test

@testset "Input validity" begin
    P = Float64[
        0 0
        0 1
    ] |> transpose

    Q = Float64[
        0 0 0
        0 0 1
    ] |> transpose

    @test_throws DimensionMismatch frechet(P, Q)
end

@testset "Correctness" begin
    P = Float64[
        0 0 0
        0 0 1
        0 0 2
    ] |> transpose

    Q = Float64[
        0 0 0
        0 0 1
        0 1 1
        0 0 1
        0 0 2
    ] |> transpose

    @test frechet(P, P) == 0
    @test frechet(Q, Q) == 0
    @test frechet(P, Q) == 1
end