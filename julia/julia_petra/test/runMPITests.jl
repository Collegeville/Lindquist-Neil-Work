using julia_petra
using Base.Test

#use distinct types
const comm = MPIComm(UInt64, UInt16, UInt32)

const pid = myPid(comm)

#only print errors from one process
if pid != 1
    #redirect_stdout()
    #redirect_stderr()
end

#try

@testset "Comm MPI Tests" begin
    include("MPICommTests.jl")
    include("MPIBlockMapTests.jl")
    include("MPIimport-export Tests.jl")
end

#catch err
#    sleep(10)
#    throw(err)
#end