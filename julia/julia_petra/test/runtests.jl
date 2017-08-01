using julia_petra
using Base.Test

include("InheritanceUtilTests.jl")
include("SerialCommTests.jl")

include("BlockMapTests.jl")
include("MultiVectorTests.jl")

# do MPI tests at the end so that other errors are found faster since the MPI tests take the longest
include("MPICommTests.jl")
