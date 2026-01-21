module TestMOIWrapper

using Test

import QPALM
import MathOptInterface as MOI

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_SolverName()
    @test MOI.get(QPALM.Optimizer(), MOI.SolverName()) == "QPALM"
    return
end

function test_supports_default_copy_to()
    @test !MOI.supports_incremental_interface(QPALM.Optimizer())
    return
end

function test_runtests()
    # This is what JuMP would construct
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(QPALM.Optimizer; with_bridge_type = Float64),
    )
    @test model.optimizer.model.model_cache isa
          MOI.Utilities.UniversalFallback{QPALM.OptimizerCache}
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("enable_dual_termination"), true)
    config = MOI.Test.Config(;
        atol = 1e-3,
        rtol = 1e-3,
        exclude = Any[
            MOI.ObjectiveBound,
            MOI.SolverVersion,
            MOI.VariableBasisStatus,
            MOI.ConstraintBasisStatus,
            MOI.DualObjectiveValue, # FIXME
            MOI.ConstraintDual, # FIXME
        ],
    )
    MOI.Test.runtests(
        model,
        config,
        exclude = [
            "INFEASIBILITY",
            "INFEASIBLE",
            "infeasible",
            "unbounded",
            "test_model_copy_to_UnsupportedAttribute", # FIXME DivideError: integer division error
        ],
    )
    return
end

end  # module TestMOIWrapper

TestMOIWrapper.runtests()
