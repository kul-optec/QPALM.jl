const LIBQPALM_PATH = "../deps/libqpalm.dylib"

const QPALM_INFTY = 1e20

const SOLUTION_PRESENT = [:Solved, :Max_iter_reached]

const status_map = Dict{Int,Symbol}(
     1  => :Solved,
    -2  => :Max_iter_reached,
    -3  => :Primal_infeasible,
    -4  => :Dual_infeasible,
    -10 => :Unsolved
)
