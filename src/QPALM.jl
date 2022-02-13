module QPALM

using QPALM_jll
const LIBQPALM_PATH = QPALM_jll.libqpalm_jll_path

const Maybe{T} = Union{T, Nothing} where T

include("const.jl")
include("types.jl")
include("wrappers.jl")

end # module
