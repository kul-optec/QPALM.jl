using Test

@testset "Update" begin

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

    model = QPALM.Model()

    QPALM.setup!(model, Q=Q, q=q, A=A, bmin=bmin, bmax=bmax)

    @testset "Settings" begin

        QPALM.update!(model; max_iter=42)

        workspace = unsafe_load(model.workspace)
        settings = unsafe_load(workspace.settings)

        @test settings.max_iter == 42

    end

    @testset "Linear term" begin

        q_new = similar(q)

        workspace = unsafe_load(model.workspace)
        data = unsafe_load(workspace.data)

        QPALM.update!(model; q=2.0*ones(n))

        unsafe_copyto!(pointer(q_new), data.q, n)
        @test q_new ≈ 2.0*ones(n)

    end

    @testset "Bounds" begin

        bmin_new = similar(bmin)
        bmax_new = similar(bmax)

        workspace = unsafe_load(model.workspace)
        data = unsafe_load(workspace.data)

        QPALM.update!(model; bmin=-2.0*ones(m))

        unsafe_copyto!(pointer(bmin_new), data.bmin, m)
        @test bmin_new ≈ -2.0*ones(m)
        unsafe_copyto!(pointer(bmax_new), data.bmax, m)
        @test bmax_new ≈ 1.0*ones(m)

        QPALM.update!(model; bmax=3.0*ones(m))

        unsafe_copyto!(pointer(bmin_new), data.bmin, m)
        @test bmin_new ≈ -2.0*ones(m)
        unsafe_copyto!(pointer(bmax_new), data.bmax, m)
        @test bmax_new ≈ 3.0*ones(m)

        QPALM.update!(model; bmin=-4.0*ones(m), bmax=4.0*ones(m))

        unsafe_copyto!(pointer(bmin_new), data.bmin, m)
        @test bmin_new ≈ -4.0*ones(m)
        unsafe_copyto!(pointer(bmax_new), data.bmax, m)
        @test bmax_new ≈ 4.0*ones(m)

    end

end
