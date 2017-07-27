
export Comm
export barrier, broadcastAll, gatherAll, sumAll, maxAll, minAll, scanSum
export myPid, numProc, createDistributor

# methods (and docs) are currently based straight off Epetra_Comm
# tpetra's equivalent seemed to be a wrapper to other Trilinos packages

# following julia's convention, n processors are labled 1 through n
# count variables are removed, since that information is contained in the arrays

"""
The base type for types that represent communication in parallel computing.
All subtypes must have the following methods, with CommImpl standing in for the subtype:

barrier(comm::CommImpl) - Each processor must wait until all processors have arrived

broadcastAll(comm::CommImpl, myvals::Array{T}, Root::Integer)::Array{T} where T
    - Takes a list of input values from the root processor and sends to all
        other processors.  The values are returned (including on the root process)

gatherAll(comm::CommImpl, myVals::Array{T})::Array{T} where T
    - Takes a list of input values from all processors and returns an ordered
        contiguous list of those values on each processor

sumAll(comm::CommImpl, partialsums::Array{T})::Array{T} where T
    - Take a list of input values from all processors and returns the sum on each
        processor.  The method +(::T, ::T)::T must exist

maxAll(comm::CommImpl, partialmaxes::Array{T})::Array{T} where T
    - Takes a list of input values from all processors and returns the max to all
        processors.  The method <(::T, ::T)::Bool must exist

minAll(comm::CommImpl, partialmins::Array{T})::Array{T} where T
    - Takes a list of input values from all processors and returns the min to all
        processors.  The method <(::T, ::T)::Bool must exist

scanSum(comm::CommImpl, myvals::Array{T})::Array{T} where T
    - Takes a list of input values from all processors, computes the scan sum and
        returns it to all processors such that processor i contains the sum of
        values from processor 1 up to and including processor i.  The method
        +(::T, ::T)::T must exist

myPid(comm::CommImpl)::Integer - Returns the process rank

numProc(comm::CommImpl)::Integer - Returns the total number of processes

createDistributor(comm::CommImpl)::Distributor - Create a distributor object

"""
abstract type Comm{GID <: Integer, PID <:Integer, LID <: Integer}
end

function Base.show(io::IO, comm::Comm)
    print(io, split(String(Symbol(typeof(comm))), ".")[2]," with PID ", myPid(comm),
                " and ", numProc(comm), " processes")
end

function broadcastAll(comm::Comm, myVal::T, root::Integer)::T where {T}
    broadcastAll(comm, [myVal], root)[1]
end

function gatherAll(comm::Comm, myVal::T)::Array{T} where {T}
    gatherAll(comm, [myVal])
end

function sumAll(comm::Comm, val::T)::T where {T}
    sumAll(comm, [val])[1]
end

function maxAll(comm::Comm, val::T)::T where {T}
    maxAll(comm, [val])[1]
end

function minAll(comm::Comm, val::T)::T where {T}
    minAll(comm, [val])[1]
end

function scanSum(comm::Comm, val::T)::T where {T}
    scanSum(comm, [val])[1]
end

