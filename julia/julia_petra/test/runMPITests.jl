using julia_petra
using Base.Test

#use distinct types
comm = MPIComm(UInt64, UInt16, UInt32)

pid = myPid(comm)

#only print errors from one process
if pid != 1
    #redirect_stdout()
    #redirect_stderr()
end

#try

include("MPICommTests.jl")
include("MPIBlockMapTests.jl")

#catch err
#    sleep(10)
#    throw(err)
#end