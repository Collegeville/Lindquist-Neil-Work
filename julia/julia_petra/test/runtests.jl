using julia_petra
using Base.Test

@. ARGS = lowercase(ARGS)
# check for command line arguments requesting parts to not be tested
noMPI = in("--mpi", ARGS) #don't run multi-process tests
noComm = in("--comm", ARGS) #don't run comm framework tests
noDataStructs = in("--data", ARGS) #don't run tests on data structures
noUtil = in("--util", ARGS) #don't run tests on Misc Utils


if !noUtil
end

if !noComm
    include("SerialCommTests.jl")
    include("Import-Export Tests.jl")
    include("BlockMapTests.jl")
end

if !noDataStructs
    include("MultiVectorTests.jl")
end

# do MPI tests at the end so that other errors are found faster since the MPI tests take the longest
if !noMPI
    include("MPITestsStarter.jl")
end
