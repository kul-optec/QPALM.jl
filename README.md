# QPALM.jl

This repository provides a Julia interface to the [QPALM](https://github.com/Benny44/QPALM)
QP solver.

## Installation

From the Julia (>= 1.0) REPL, do

```julia
] add https://github.com/kul-forbes/QPALM.jl
```

In order to use QPALM.jl, you will need the binaries, which you can find [here](https://bintray.com/benny44/generic/QPALM/1.0#files/). You will end up having a library file named `libqpalm.so`
(on Linux), `libqpalm.dylib` (on macOS), or `libqpalm.dll` (on Windows), which
you will have to copy as follows:

```bash
# on Linux
cp libqpalm.so <PATH-TO-QPALM.jl>/deps/lib/

# on macOS
cp libqpalm.dylib <PATH-TO-QPALM.jl>/deps/lib/

# on Windows
cp libqpalm.dll <PATH-TO-QPALM.jl>/deps/lib/
```

## Usage

Given the QP

```
minimize    (1/2) x' Q x + q' x
subject to  bmin <= A x <= bmax
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
so please refer [QPALM's documentation](https://benny44.github.io/QPALM/structQPALMSettings.html)
on how to set these and their semantics. Leaving the settings unspecified
will run the solver with default options.
