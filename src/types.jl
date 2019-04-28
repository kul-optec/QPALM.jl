using SparseArrays
using SuiteSparse: CHOLMOD

# Integer type from C
if Sys.WORD_SIZE == 64   # 64bit system
    const Cc_int = Clonglong
else  # 32bit system
    const Cc_int = Cint
end

struct Carray_element
    x::Cdouble
    i::Cc_int
end

# struct Ccsc
#     nzmax::Cc_int
#     m::Cc_int
#     n::Cc_int
#     p::Ptr{Cc_int}
#     i::Ptr{Cc_int}
#     x::Ptr{Cdouble}
#     nz::Cc_int
# end
#
#
# struct ManagedCcsc
#     nzmax::Cc_int
#     m::Cc_int
#     n::Cc_int
#     p::Vector{Cc_int}
#     i::Vector{Cc_int}
#     x::Vector{Cdouble}
#     nz::Cc_int
#
# end
#
# # Construct ManagedCcsc matrix from SparseMatrixCSC
# function ManagedCcsc(M::SparseMatrixCSC)
#
#     # Get dimensions
#     m = M.m
#     n = M.n
#
#     # Get vectors of data, rows indices and column pointers
#     x = convert(Array{Float64,1}, M.nzval)
#     # C is 0 indexed
#     i = convert(Array{Cc_int,1}, M.rowval .- 1)
#     # C is 0 indexed
#     p = convert(Array{Cc_int,1}, M.colptr .- 1)
#
#     # Create new ManagedCcsc matrix
#     ManagedCcsc(length(M.nzval), m, n, p, i, x, -1)
# end
#
# # Returns an Ccsc matrix. The vectors are *not* GC tracked in the struct.
# # Use this only when you know that the managed matrix will outlive the Ccsc
# # matrix.
# Ccsc(m::ManagedCcsc) =
#     Ccsc(m.nzmax, m.m, m.n, pointer(m.p), pointer(m.i), pointer(m.x), m.nz)

struct Data
    n::Cc_int
    m::Cc_int
    Q::Ptr{CHOLMOD.C_Sparse}
    A::Ptr{CHOLMOD.C_Sparse}
    q::Ptr{Cdouble}
    bmin::Ptr{Cdouble}
    bmax::Ptr{Cdouble}
end

struct Solution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
end

struct CInfo
    iter::Cc_int
    iter_out::Cc_int
    status::NTuple{32,Cchar}
    status_val::Cc_int
    pri_res_norm::Cdouble
    dua_res_norm::Cdouble
    dua2_res_norm::Cdouble
    setup_time::Cdouble
    solve_time::Cdouble
    run_time::Cdouble
end

struct Settings
    max_iter::Cc_int
    eps_abs::Cdouble
    eps_rel::Cdouble
    eps_abs_in::Cdouble
    eps_rel_in::Cdouble
    rho::Cdouble
    eps_prim_inf::Cdouble
    eps_dual_inf::Cdouble
    theta::Cdouble
    delta::Cdouble
    tau_init::Cdouble
    proximal::Cc_int
    gamma::Cdouble
    gamma_upd::Cdouble
    gamma_max::Cdouble
    scaling::Cc_int
    verbose::Cc_int
    warm_start::Cc_int
end

struct LBFGS
    curridx::Cc_int
    currmem::Cc_int
    reset_lbfgs::Cc_int
    s::Ptr{Cdouble}
    y::Ptr{Cdouble}
    ys::Cdouble
    Sbuffer::Ptr{Cdouble}
    Ybuffer::Ptr{Cdouble}
    YSbuffer::Ptr{Cdouble}
    H::Cdouble
    alpha::Ptr{Cdouble}
    q::Ptr{Cdouble}
end

struct Workspace
    data::Ptr{QPALM.Data}

    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    Ax::Ptr{Cdouble}
    Qx::Ptr{Cdouble}
    Aty::Ptr{Cdouble}
    x_prev::Ptr{Cdouble}
    initialized::Cc_int

    temp_m::Ptr{Cdouble}
    temp_n::Ptr{Cdouble}
    sigma::Ptr{Cdouble}
    Axys::Ptr{Cdouble}
    z::Ptr{Cdouble}
    pri_res::Ptr{Cdouble}
    pri_res_in::Ptr{Cdouble}
    yh::Ptr{Cdouble}
    Atyh::Ptr{Cdouble}
    df::Ptr{Cdouble}
    x0::Ptr{Cdouble}
    xx0::Ptr{Cdouble}
    dphi::Ptr{Cdouble}
    neg_dphi::Ptr{Cdouble}
    dphi_prev::Ptr{Cdouble}
    d::Ptr{Cdouble}

    tau::Cdouble
    Qd::Ptr{Cdouble}
    Ad::Ptr{Cdouble}
    sqrt_sigma::Ptr{Cdouble}
    sqrt_delta::Cdouble
    eta::Cdouble
    beta::Cdouble
    delta::Ptr{Cdouble}
    alpha::Ptr{Cdouble}
    temp_2m::Ptr{Cdouble}
    delta2::Ptr{Cdouble}
    delta_alpha::Ptr{Cdouble}
    s::Ptr{Carray_element}
    index_L::Cc_int
    index_P::Cc_int
    index_J::Cc_int

    eps_pri::Cdouble
    eps_dua::Cdouble
    eps_dua_in::Cdouble

    delta_y::Ptr{Cdouble}
    Atdelta_y::Ptr{Cdouble}

    delta_x::Ptr{Cdouble}
    Qdelta_x::Ptr{Cdouble}
    Adelta_x::Ptr{Cdouble}

    D_temp::Ptr{Cdouble}
    E_temp::Ptr{Cdouble}

    chol::Ptr{Nothing}
    settings::Ptr{QPALM.Settings}
    scaling::Ptr{Nothing}
    solution::Ptr{QPALM.Solution}
    info::Ptr{QPALM.CInfo}

    timer::Ptr{Nothing}
end
