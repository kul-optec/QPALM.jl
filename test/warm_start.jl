using Test

@testset "Warm start" begin

    using QPALM
    using LinearAlgebra
    using SparseArrays

    n, m = 2, 3
    act = 3

    Q = sparse(I, n, n)
    A = sparse(ones(m, n))
    q = ones(n)
    bmin = -ones(m)
    bmax = ones(m)


    x0 = ones(n)
    y0 = ones(m)

    @testset "Not setup" begin
        model = QPALM.Model()

        @test_throws ErrorException QPALM.warm_start!(model, x_warm_start=x0)
    end

    @testset "Normal usage (1)" begin
        model = QPALM.Model()

        QPALM.setup!(model, Q=Q, q=q, A=A, bmin=bmin, bmax=bmax)
        QPALM.warm_start!(model, x_warm_start=x0, y_warm_start=y0)

        results = QPALM.solve!(model)
        @test results.info.status == :Solved
    end

    @testset "Normal usage (2)" begin
        model = QPALM.Model()

        QPALM.setup!(model, Q=Q, q=q, A=A, bmin=bmin, bmax=bmax)
        QPALM.warm_start!(model, y_warm_start=y0)
        QPALM.warm_start!(model, x_warm_start=x0)

        results = QPALM.solve!(model)
        @test results.info.status == :Solved
    end

end
