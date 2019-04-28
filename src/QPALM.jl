module QPALM

const Maybe{T} = Union{T, Nothing} where T

include("const.jl")
include("types.jl")
include("wrappers.jl")

end # module
