using SparseArrays
linsys = "ladel" #"ladel" or "cholmod"

if linsys == "cholmod"
    using SuiteSparse: CHOLMOD
end

const Cc_int = Int64
const Cc_uint = UInt64
const Cc_float = Cdouble


struct Carray_element
    x::Cc_float
    i::Csize_t
end

struct ladel_sparse_matrix
    nzmax::Cc_int
    nrow::Cc_int
    ncol::Cc_int
    p::Ptr{Cc_int}
    i::Ptr{Cc_int}
    x::Ptr{Cc_float}
    nz::Ptr{Cc_int}
    values::Cc_int
    symmetry::Cc_int
end

if linsys == "ladel"
    sparse_matrix_ptr = Ref{QPALM.ladel_sparse_matrix}
elseif linsys == "cholmod"
    sparse_matrix_ptr = Ptr{CHOLMOD.C_Sparse}
else
    error("Unrecognised linsys")
end

struct Data
    n::Csize_t
    m::Csize_t
    Q::sparse_matrix_ptr
    A::sparse_matrix_ptr
    q::Ptr{Cc_float}
    c::Cc_float
    bmin::Ptr{Cc_float}
    bmax::Ptr{Cc_float}
end

struct Solution
    x::Ptr{Cc_float}
    y::Ptr{Cc_float}
end

struct CInfo
    iter::Cc_int
    iter_out::Cc_int
    status::NTuple{32,Cchar}
    status_val::Cc_int
    pri_res_norm::Cc_float
    dua_res_norm::Cc_float
    dua2_res_norm::Cc_float
    objective::Cc_float
    dual_objective::Cc_float
    setup_time::Cc_float
    solve_time::Cc_float
    run_time::Cc_float
end

struct Settings
    max_iter::Cc_int
    inner_max_iter::Cc_int
    eps_abs::Cc_float
    eps_rel::Cc_float
    eps_abs_in::Cc_float
    eps_rel_in::Cc_float
    rho::Cc_float
    eps_prim_inf::Cc_float
    eps_dual_inf::Cc_float
    theta::Cc_float
    delta::Cc_float
    sigma_max::Cc_float
    sigma_init::Cc_float
    proximal::Cc_int
    gamma_init::Cc_float
    gamma_upd::Cc_float
    gamma_max::Cc_float
    scaling::Cc_int
    nonconvex::Cc_int
    verbose::Cc_int
    print_iter::Cc_int
    warm_start::Cc_int
    reset_newton_iter::Cc_int
    enable_dual_termination::Cc_int
    dual_objective_limit::Cc_float
    time_limit::Cc_float
    ordering::Cc_int
    factorization_method::Cc_int
    max_rank_update::Cc_int
    max_rank_update_fraction::Cc_float
end

struct Workspace
    data::Ptr{QPALM.Data}

    x::Ptr{Cc_float}
    y::Ptr{Cc_float}
    Ax::Ptr{Cc_float}
    Qx::Ptr{Cc_float}
    Aty::Ptr{Cc_float}
    x_prev::Ptr{Cc_float}
    initialized::Cc_int

    temp_m::Ptr{Cc_float}
    temp_n::Ptr{Cc_float}
    sigma::Ptr{Cc_float}
    sigma_inv::Ptr{Cc_float}
    sqrt_sigma_max::Cc_float
    nb_sigma_changed::Cc_int
    gamma::Cc_float
    gamma_maxed::Cc_int
    Axys::Ptr{Cc_float}
    z::Ptr{Cc_float}
    pri_res::Ptr{Cc_float}
    pri_res_in::Ptr{Cc_float}
    yh::Ptr{Cc_float}
    Atyh::Ptr{Cc_float}
    df::Ptr{Cc_float}
    x0::Ptr{Cc_float}
    xx0::Ptr{Cc_float}
    dphi::Ptr{Cc_float}
    neg_dphi::Ptr{Cc_float}
    dphi_prev::Ptr{Cc_float}
    d::Ptr{Cc_float}

    tau::Cc_float
    Qd::Ptr{Cc_float}
    Ad::Ptr{Cc_float}
    sqrt_sigma::Ptr{Cc_float}
    sqrt_delta::Cc_float
    eta::Cc_float
    beta::Cc_float
    delta::Ptr{Cc_float}
    alpha::Ptr{Cc_float}
    temp_2m::Ptr{Cc_float}
    delta2::Ptr{Cc_float}
    delta_alpha::Ptr{Cc_float}
    s::Ptr{Carray_element}
    index_L::Ptr{Cc_int}
    index_P::Ptr{Cc_int}
    index_J::Ptr{Cc_int}

    eps_pri::Cc_float
    eps_dua::Cc_float
    eps_dua_in::Cc_float
    eps_abs_in::Cc_float
    eps_rel_in::Cc_float

    delta_y::Ptr{Cc_float}
    Atdelta_y::Ptr{Cc_float}

    delta_x::Ptr{Cc_float}
    Qdelta_x::Ptr{Cc_float}
    Adelta_x::Ptr{Cc_float}

    D_temp::Ptr{Cc_float}
    E_temp::Ptr{Cc_float}

    solver::Ptr{Nothing}
    settings::Ptr{QPALM.Settings}
    scaling::Ptr{Nothing}
    solution::Ptr{QPALM.Solution}
    info::Ptr{QPALM.CInfo}

    timer::Ptr{Nothing}
end
