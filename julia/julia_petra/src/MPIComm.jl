
import MPI

"""
    MPIComm()
An implementation of Comm using MPI
The no argument constructor uses MPI.COMM_WORLD
"""
type MPIComm{GID <: Integer, PID <:Integer, LID <: Integer} <: Comm{GID, PID, LID}
    mpiComm::MPI.Comm
end

function MPIComm()
    MPIComm(MPI.COMM_WORLD)
end


function barrier(comm::MPIComm)
    MPI.Barrier(comm.mpiComm)
end

function broadcastAll(comm::MPIComm{GID, PID}, myvals::Array{T}, root::PID)::Array{T} where GID <: Integer where PID <: Integer
    vals = copy(myvals)
    
    MPI.Bcast!(vals, root, comm.mpiComm)
end

function gatherAll(comm::MPIComm, myVals::Array{T})::Array{T} where T
    MPI.Allgather(myVals, comm)
end

function sumAll(comm::MPIComm, partialsums::Array{T})::Array{T} where T
    MPI.allreduce(partialsums, +, comm.mpiComm)
end

function maxAll(comm::MPIComm, partialmaxes::Array{T})::Array{T} where T
    MPI.allreduce(partialsums, max, comm.mpiComm)
end

function minAll(comm::MPIComm, partialmins::Array{T})::Array{T} where T
    MPI.allreduce(partialsums, min, comm.mpiComm)
end

function scanSum(comm::MPIComm, myvals::Array{T})::Array{T} where T
    MPI.scan(myvals, length(myvals), +, comm.mpiComm)
end

function myPid(comm::MPIComm{GID, PID})::PID where GID <: Integer where PID <: Integer
    MPI.Comm_rank(comm.mpiComm) + 1
end

function numProc(comm::MPIComm)::Integer
    MPI.Comm_size(comm.mpiComm)
end

#TODO implement this
createDistributor(comm::MPIComm)::Distributor - Create a distributor object
