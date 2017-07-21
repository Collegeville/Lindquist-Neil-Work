
type SerialComm <: Comm
end

const serialCommInstance = SerialComm()

function SerialComm()
    serialCommInstance
end


# most of these functions are no-ops or identify functions since there is only
# one processor

function barrier(comm::SerialComm)
end


function broadcast(comm::SerialComm, myvals::Array{T}, root::Integer)::Array{T} where T
    if root != 1
        error("SerialComm can only accept PID of 1")
    end
    myVals
end

function gatherAll(comm::SerialComm, myVals::Array{T})::Array{T} where T
    myVals
end

function sumAll(comm::SerialComm, partialsums::Array{T})::Array{T} where T
    partialsums
end

function maxAll(comm::SerialComm, partialmaxes::Array{T})::Array{T} where T
    partialmaxes
end

function minAll(comm::SerialComm, partialmins::Array{T})::Array{T} where T
    partialmins
end

function scanSum(comm::SerialComm, myvals::Array{T})::Array{T} where T
    myvals
end

function myPid(comm::SerialComm)::UInt8
    1
end

function numProc(comm::SerialComm)::UInt8
    1
end

function createDistributor(comm::SerialComm)::SerialDistributor
    SerialDistributor()
end

function createDirectory(comm::SerialComm, map::Map)::SerialDirectory
    BasicDirectory(map)
end