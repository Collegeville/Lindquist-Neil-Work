using julia_petra
using Base.Test

#use distinct types
comm = MPIComm(UInt64, UInt16, UInt32)

#only print errors from one process
if myPid(comm) != 1
    #redirect_stdout()
    #redirect_stderr()
end

try
    @test 4 == numProc(comm)
    @test isa(numProc(comm), UInt16)
    
    @test 1 <= myPid(comm) <= 4
    @test isa(myPid(comm), UInt16)
    
    @test [1, 8] == broadcastAll(comm, [myPid(comm), 8], 1)
    @test [2, 5] == broadcastAll(comm, [myPid(comm), 5], 2)
    @test [3, 7] == broadcastAll(comm, [myPid(comm), 7], 3)
    @test [4, 6] == broadcastAll(comm, [myPid(comm), 6], 4)
    
    @test [1, 2, 3, 4] == gatherAll(comm, [myPid(comm)])
    @test ([1, 2, 3, 2, 4, 6, 3, 6, 9, 4, 8, 12] 
            == gatherAll(comm, [myPid(comm), myPid(comm)*2, myPid(comm)*3]))
    
    #check for hangs and such, hard to test if all processes are at the same spot
    barrier(comm)
    
    @test [10] == sumAll(comm, [myPid(comm)])
    @test [32, 12, 10, 8] == sumAll(comm, [8, 3, myPid(comm), 2])
    
    @test [4] == maxAll(comm, [myPid(comm)])
    @test [4, -1, 8] == maxAll(comm, [myPid(comm), -Int(myPid(comm)), 8])
    
    @test [1] == minAll(comm, [myPid(comm)])
    @test [1, -4, 6] == minAll(comm, [myPid(comm), -Int(myPid(comm)), 6])
    
    @test [sum(1:myPid(comm))] == scanSum(comm, [myPid(comm)])
    @test ([myPid(comm)*5, sum(-2:-2:(-2*Int(myPid(comm)))), myPid(comm)*3] 
            == scanSum(comm, [5, -2*Int(myPid(comm)), 3]))
    
    #TODO test distributor
    dist = createDistributor(comm)
    
    @test isa(dist, Distributor{UInt64, UInt16, UInt32})
    
catch err
    throw(err)
end