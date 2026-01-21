# QPALM.jl

This repository provides a Julia interface to the [QPALM](https://github.com/kul-optec/QPALM)
QP solver.

## Installation

In the Julia console, press `]` to enter the Pkg REPL and install QPALM using:
```sh
add QPALM
```

To test the correct installation of the Julia interface, you can run the unit tests:
```julia
include("test/runtests.jl")
```

## Usage

Given the QP

```
minimize        ½ xᵀQx + qᵀx

subject to      l ≤ Ax ≤ u
```

this is solved by

```julia
using QPALM
model = QPALM.Model()
QPALM.setup!(model, Q=Q, q=q, A=A, bmin=bmin, bmax=bmax, settings...)
results = QPALM.solve!(model)
```

where `settings...` are keyword arguments specifying the solver options
to use. They have the same name and type as the underlying C API,
so please refer to [QPALM's documentation](https://kul-optec.github.io/QPALM/Doxygen/structQPALMSettings.html)
on how to set these and their semantics. Leaving the settings unspecified
will run the solver with default options.

## Use with JuMP

> [!INFO]
> For the dual values to be computed and be available with `JuMP.dual` and `JuMP.objective_dual`,
> you need set the attribute `enable_dual_termination` before calling `JuMP.optimize!` as shown below.

To use QPALM with JuMP, use `QPALM.Optimizer`:

```julia
using JuMP, QPALM
model = Model(QPALM.Optimizer)
# If duals are needed
set_attribute(model, "enable_dual_termination", true)
```
