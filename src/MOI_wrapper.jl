import MathOptInterface as MOI

# Inspired from `Clp.jl/src/MOI_wrapper.jl`
MOI.Utilities.@product_of_sets(
    _LPProductOfSets,
    MOI.EqualTo{T},
    MOI.GreaterThan{T},
    MOI.LessThan{T},
    MOI.Interval{T},
)

const OptimizerCache = MOI.Utilities.GenericModel{
    Float64,
    MOI.Utilities.ObjectiveContainer{Float64},
    MOI.Utilities.VariablesContainer{Float64},
    MOI.Utilities.MatrixOfConstraints{
        Float64,
        MOI.Utilities.MutableSparseMatrixCSC{
            Float64,
            Int64,
            MOI.Utilities.OneBasedIndexing,
        },
        MOI.Utilities.Hyperrectangle{Float64},
        _LPProductOfSets{Float64},
    },
}

Base.show(io::IO, ::Type{OptimizerCache}) = print(io, "QPALM.OptimizerCache")

const BOUND_SETS = Union{
    MOI.EqualTo{Float64},
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},
    MOI.Interval{Float64},
}

"""
    Optimizer()

Create a new QPALM optimizer.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    model::Union{Nothing,Model}
    start_x::Union{Nothing,Vector{Float64}}
    start_y::Union{Nothing,Vector{Float64}}
    m::Int
    n::Int
    result::Union{Nothing,Results}
    has_dual::Bool
    objective_constant::Float64
    max_sense::Bool
    silent::Bool
    options::Dict{Symbol,Any}

    function Optimizer()
        return new(
            nothing,
            nothing,
            nothing,
            0,
            0,
            nothing,
            false,
            0,
            false,
            false,
            Dict{Symbol,Any}(),
        )
    end
end

function MOI.default_cache(::Optimizer, ::Type)
    return MOI.Utilities.UniversalFallback(OptimizerCache())
end

# ====================
#   empty functions
# ====================

function MOI.is_empty(optimizer::Optimizer)
    return isnothing(optimizer.model)
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.model = nothing
    optimizer.start_x = nothing
    optimizer.start_y = nothing
    optimizer.result = nothing
    return
end

MOI.get(::Optimizer, ::MOI.SolverName) = "QPALM"

# MOI.RawOptimizerAttribute

function MOI.supports(::Optimizer, param::MOI.RawOptimizerAttribute)
    return hasfield(Settings, Symbol(param.name))
end

function MOI.set(optimizer::Optimizer, param::MOI.RawOptimizerAttribute, value)
    if !MOI.supports(optimizer, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    optimizer.options[Symbol(param.name)] = value
    return
end

function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
    if !MOI.supports(optimizer, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    return optimizer.parameters[Symbol(param.name)]
end

# MOI.Silent

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

# ========================================
#   Supported constraints and objectives
# ========================================

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.ScalarAffineFunction{Float64}},
    ::Type{<:BOUND_SETS},
)
    return true
end

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{F},
) where {F<:Union{MOI.ScalarAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64}}}
    return true
end

# ===============================
#   Optimize and post-optimize
# ===============================

function _flip_sense(optimizer::Optimizer, obj)
    return optimizer.max_sense ? -obj : obj
end

function terms_to_vector(model, terms::Vector{MOI.ScalarAffineTerm{Float64}}, n)
    c = zeros(n)
    for term in terms
        c[term.variable.value] += _flip_sense(model, term.coefficient)
    end
    return c
end

function terms_to_matrix(model, terms::Vector{MOI.ScalarQuadraticTerm{Float64}}, n)
    I = Int[]
    J = Int[]
    V = Float64[]
    for term in terms
        coef = _flip_sense(model, term.coefficient)
        push!(I, term.variable_1.value)
        push!(J, term.variable_2.value)
        push!(V, coef)
        if term.variable_1.value != term.variable_2.value
            push!(I, term.variable_2.value)
            push!(J, term.variable_1.value)
            push!(V, coef)
        end
    end
    return sparse(I, J, V, n, n)
end

function MOI.copy_to(dest::Optimizer, src::OptimizerCache)
    dest.model = Model()
    Ab = src.constraints
    A = convert(SparseMatrixCSC{Float64,Int64}, Ab.coefficients)
    row_bounds = src.constraints.constants
    dest.max_sense = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    F = MOI.get(src, MOI.ObjectiveFunctionType())
    obj = MOI.get(src, MOI.ObjectiveFunction{F}())
    if F == MOI.ScalarAffineFunction{Float64}
        q = terms_to_vector(dest, obj.terms, A.n)
        Q = nothing
    else
        @assert F == MOI.ScalarQuadraticFunction{Float64}
        q = terms_to_vector(dest, obj.affine_terms, A.n)
        Q = terms_to_matrix(dest, obj.quadratic_terms, A.n)
    end
    dest.objective_constant = MOI.constant(obj)
    options = dest.options
    if dest.silent
        options = copy(options)
        options[:verbose] = 0
    end
    @show size(A)
#    if size(A, 1) == 0 && haskey(options, :enable_dual_termination) && options[:enable_dual_termination]
#        # Otherwise throws: LADEL ERROR: MATRIX (POSSIBLY) NOT FULL RANK (diagonal element of 0.000000e+00)
#        @warn("Disabling `enable_dual_termination`, not supported for a problem with no constraint")
#        options[:enable_dual_termination] = false
#    end
    setup!(
        dest.model;
        Q,
        q,
        A,
        bmin=row_bounds.lower,
        bmax=row_bounds.upper,
        options...,
    )
    dest.has_dual = false
    dest.m, dest.n = size(A)
    dest.start_x = nothing
    dest.start_y = nothing
    dest.result = nothing
    return MOI.Utilities.identity_index_map(src)
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value,
)
    if !isnothing(value) && isnothing(optimizer.start_x)
        optimizer.start_x = zeros(optimizer.n)
    end
    optimizer.start_x[vi.value] = something(value, 0)
    return
end

function MOI.set(
    optimizer::Optimizer,
    ::MOI.ConstraintDualStart,
    vi::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}},
    value,
)
    if !isnothing(value) && isnothing(optimizer.start_y)
        optimizer.start_x = zeros(optimizer.m)
    end
    optimizer.start_y[vi.value] = something(value, 0)
    return
end

function MOI.copy_to(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache},
)
    index_map = MOI.copy_to(dest, src.model)
    MOI.Utilities.pass_attributes(
        dest,
        src,
        index_map,
        MOI.get(src, MOI.ListOfVariableIndices()),
    )
    # The `ObjectiveSense` and `ObjectiveFunction`
    # have already been handled, we need to gracefully
    # error if there are other ones
    MOI.Utilities.pass_attributes(
        dest,
        MOI.Utilities.ModelFilter(src) do attr
            return !(attr isa MOI.ObjectiveSense) &&
                !(attr isa MOI.ObjectiveFunction)
        end,
        index_map,
    )
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        MOI.Utilities.pass_attributes(
            dest,
            src,
            index_map,
            MOI.get(src, MOI.ListOfConstraintIndices{F,S}()),
        )
    end
    return index_map
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    cache = MOI.default_cache(dest, Float64)
    index_map = MOI.copy_to(cache, src)
    MOI.copy_to(dest, cache)
    return index_map, false
end

function MOI.optimize!(dest::Optimizer)
    dest.result = solve!(dest.model)
    return
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec)
    return optimizer.result.info.run_time
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    if isnothing(optimizer.result)
        return "Optimize not called"
    else
        return string(optimizer.result.info.status)
    end
end

const _TERMINATION_STATUS_MAP = Dict{Int,MOI.TerminationStatusCode}(
     0  => MOI.OTHER_ERROR,
     1  => MOI.OPTIMAL,
     2  => MOI.OTHER_ERROR, # ???
    -2  => MOI.ITERATION_LIMIT,
    -3  => MOI.INFEASIBLE,
    -4  => MOI.DUAL_INFEASIBLE,
    -5  => MOI.TIME_LIMIT,
    -10 => MOI.OPTIMIZE_NOT_CALLED,
)

# Implements getter for result value and statuses
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    return isnothing(optimizer.result) ? MOI.OPTIMIZE_NOT_CALLED :
           _TERMINATION_STATUS_MAP[optimizer.result.info.status_val]
end

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    return _flip_sense(optimizer, optimizer.result.info.objective) + optimizer.objective_constant
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    if !optimizer.has_dual
        return NaN
    end
    return _flip_sense(optimizer, optimizer.result.info.dual_objective) + optimizer.objective_constant
end

const _PRIMAL_STATUS_MAP = Dict{Int,MOI.ResultStatusCode}(
     0  => MOI.UNKNOWN_RESULT_STATUS,
     1  => MOI.FEASIBLE_POINT,
     2  => MOI.UNKNOWN_RESULT_STATUS, # ???
    -2  => MOI.UNKNOWN_RESULT_STATUS,
    -3  => MOI.INFEASIBLE_POINT,
    -4  => MOI.INFEASIBILITY_CERTIFICATE,
    -5  => MOI.UNKNOWN_RESULT_STATUS,
    -10 => MOI.NO_SOLUTION,
)

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return _PRIMAL_STATUS_MAP[optimizer.result.info.status_val]
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.result.x[vi.value]
end

const _DUAL_STATUS_MAP = Dict{Int,MOI.ResultStatusCode}(
     0  => MOI.UNKNOWN_RESULT_STATUS,
     1  => MOI.FEASIBLE_POINT,
     2  => MOI.UNKNOWN_RESULT_STATUS, # ???
    -2  => MOI.UNKNOWN_RESULT_STATUS,
    -3  => MOI.INFEASIBILITY_CERTIFICATE,
    -4  => MOI.INFEASIBLE_POINT,
    -5  => MOI.UNKNOWN_RESULT_STATUS,
    -10 => MOI.NO_SOLUTION,
)

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if !optimizer_has_dual || attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    return _DUAL_STATUS_MAP[optimizer.result.info.status_val]
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},S},
) where {S}
    MOI.check_result_index_bounds(optimizer, attr)
    if !optimizer.has_dual
        return NaN
    end
    return _flip_sense(optimizer, optimizer.result.y[ci.value])
end

function MOI.get(optimizer::Optimizer, ::MOI.ResultCount)
    if isnothing(optimizer.result)
        return 0
    else
        return 1
    end
end
