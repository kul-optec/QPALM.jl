import MathOptInterface as MOI

# Inspired from `Zaphod.jl/src/MOI_wrapper.jl`
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
    objective_constant::Float64
    max_sense::Bool
    silent::Bool
    options::Dict{Symbol,Any}

    function Optimizer()
        return new(
            nothing,
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
    return isempty(optimizer.model)
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.model = nothing
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
    setfield!(optimizer.options, Symbol(param.name), value)
    return
end

function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
    if !MOI.supports(optimizer, param)
        throw(MOI.UnsupportedAttribute(param))
    end
    return getfield(optimizer.parameters, Symbol(param.name))
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
    ::Type{MOI.VariableIndex},
    ::Type{<:BOUND_SETS},
)
    return true
end

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
    ::MOI.ObjectiveFunction{<:Union{MOI.ScalarAffineFunction{Float64},MOI.ScalarQuadraticFunction{Float64}}},
)
    return true
end

# ===============================
#   Optimize and post-optimize
# ===============================

function _flip_sense(optimizer::Optimizer, obj)
    return optimizer.max_sense ? -obj : obj
end

function terms_to_vector(terms::Vector{MOI.ScalarAffineTerm}, n)
    c = zeros(n)
    for term in terms
        c[term.variable.value] += term.coefficient
    end
end

function terms_to_matrix(terms::Vector{MOI.ScalarAffineTerm}, n)
    I = Int[]
    J = Int[]
    V = Float64[]
    for term in terms
        push!(I, term.variable_1)
        push!(J, term.variable_2)
        push!(V, term.coefficient)
        if term.variable_1 != term.variable_2
            push!(I, term.variable_2)
            push!(J, term.variable_1)
            push!(V, term.coefficient)
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
        q = terms_to_vector(obj.terms, A.n)
        Q = nothing
    else
        @assert F == MOI.ScalarQuadraticFunction{Float64}
        q = terms_to_vector(obj.affine_terms, A.n)
        Q = terms_to_matrix(obj.quadratic_terms, A.n)
    end
    dest.objective_constant = MOI.constant(obj)
    setup!(
        model;
        Q,
        q,
        A,
        bmin=row_bounds.lower,
        bmax=row_bounds.upper,
        dest.options...,
    )
    return MOI.identity_index_map(src)
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    cache = MOI.default_cache(dest, Float64)
    index_map = MOI.copy_to(cache, src)
    MOI.copy_to(dest, cache)
    return index_map, false
end

function MOI.optimize!(dest::Optimizer)
    solve!(dest.model)
end
