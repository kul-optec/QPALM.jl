using QPALM
using Random
using LinearAlgebra
using SparseArrays
using SuiteSparse: CHOLMOD

using Test

Random.seed!(0)

F = randn(5, 3)
Q = F*F'
q = randn(5)
A = randn(10, 5)
bmin = -rand(10)
bmax = +rand(10)

model = QPALM.Model()
QPALM.setup!(model, Q=Q, q=q, A=A, bmin=bmin, bmax=bmax)
result = QPALM.solve!(model)

@test result.info.status == :Solved
