using julia_petra
using Base.Test

@. ARGS = lowercase(ARGS)
# check for command line arguments requesting parts to not be tested
const noMPI = in("--mpi", ARGS) #don't run multi-process tests
const noComm = in("--comm", ARGS) #don't run comm framework tests
const noDataStructs = in("--data", ARGS) #don't run tests on data structures
const noUtil = in("--util", ARGS) #don't run tests on Misc Utils


@testset "Serial tests" begin

    if !noUtil
        @testset "Util Tests" begin
        end
    end

    if !noComm
        @testset "Comm Tests" begin
            include("SerialCommTests.jl")
            include("Import-Export Tests.jl")
            include("BlockMapTests.jl")
        end
    end

    if !noDataStructs
        @testset "Data Structure Tests" begin
            include("MultiVectorTests.jl")
            multiVectorTests(SerialComm{UInt64, UInt16, UInt32}())
            
            include("LocalCRSGraphTests.jl")
            include("CRSGraphTests.jl")
            include("CSRMatrixTests.jl")
        end
    end
end

# do MPI tests at the end so that other errors are found faster since the MPI tests take the longest
if !noMPI
    include("MPITestsStarter.jl")
end
