using julia_petra
using Base.Test

include("InheritanceUtilTests.jl")
include("SerialCommTests.jl")

include("Import-Export Tests.jl")

include("BlockMapTests.jl")
include("MultiVectorTests.jl")

# do MPI tests at the end so that other errors are found faster since the MPI tests take the longest
include("MPITestsStarter.jl")
