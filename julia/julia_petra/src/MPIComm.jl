
import MPI

export MPIComm

#TODO make MPIComm immutable then use atexit(f) to finalize MPI
"""
    MPIComm()
    MPIComm(comm::MPI.Comm)

An implementation of Comm using MPI
The no argument constructor uses MPI.COMM_WORLD
"""
struct MPIComm{GID <: Integer, PID <:Integer, LID <: Integer} <: Comm{GID, PID, LID}
    mpiComm::MPI.Comm
end

function MPIComm(GID::Type, PID::Type, LID::Type)
    MPI.Init()
    comm = MPIComm{GID, PID, LID}(MPI.COMM_WORLD)
    
    atexit(() -> MPI.Finalize())
    
    comm
end


function barrier(comm::MPIComm)
    MPI.Barrier(comm.mpiComm)
end

function broadcastAll(comm::MPIComm, myvals::Array{T}, root::Integer)::Array{T} where T
    vals = copy(myvals)
    result = MPI.Bcast!(vals, root-1, comm.mpiComm)
    result
end

function gatherAll(comm::MPIComm, myVals::Array{T})::Array{T} where T
    MPI.Allgather(myVals, comm.mpiComm)
end

function sumAll(comm::MPIComm, partialsums::Array{T})::Array{T} where T
    MPI.allreduce(partialsums, +, comm.mpiComm)
end

function maxAll(comm::MPIComm, partialmaxes::Array{T})::Array{T} where T
    MPI.allreduce(partialmaxes, max, comm.mpiComm)
end

function maxAll(comm::MPIComm, partialmaxes::Array{Bool})::Array{Bool}
    Array{Bool}(maxAll(comm, Array{UInt8}(partialmaxes)))
end

function minAll(comm::MPIComm, partialmins::Array{T})::Array{T} where T
    MPI.allreduce(partialmins, min, comm.mpiComm)
end

function minAll(comm::MPIComm, partialmins::Array{Bool})::Array{Bool}
    Array{Bool}(minAll(comm, Array{UInt8}(partialmins)))
end

function scanSum(comm::MPIComm, myvals::Array{T})::Array{T} where T
    MPI.Scan(myvals, length(myvals), MPI.SUM, comm.mpiComm)
end

function myPid(comm::MPIComm{GID, PID})::PID where {GID <: Integer, PID <: Integer}
    MPI.Comm_rank(comm.mpiComm) + 1
end

function numProc(comm::MPIComm{GID, PID})::PID where {GID <: Integer, PID <:Integer}
    MPI.Comm_size(comm.mpiComm)
end

function createDistributor(comm::MPIComm{GID, PID, LID})::MPIDistributor{GID, PID, LID}  where {GID <: Integer, PID <: Integer, LID <: Integer}
    MPIDistributor(comm)
end
