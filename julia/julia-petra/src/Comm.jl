
# methods (and docs) are currently based straight off Epetra_Comm
# tpetra's equivalent seemed to be a wrapper to other Trilinos packages

# following julia's convention, n processors are labled 1 through n
# count variables are removed, since that information is contained in the arrays

# DECISION: limit T to Array{S} where S <: Number? or continue to allow any value?
# DECISION: use specific datatype for process rank? UInt16 (UInt8?) is big enough
    # Apply decision to return value specs for getDirectoryEntries
"""
The base type for types that represent communication in parallel computing.
All subtypes must have the following methods, with CommImpl standing in for the subtype:

barrier(comm::CommImpl) - Each processor must wait until all processors have arrived

broadcast(comm::CommImpl, myvals::T, Root::Integer)::T where T
    - Takes a list of input values from the root processor and sends to all
        other processors.  The values are returned (including on the root process)

gatherAll(comm::CommImpl, myVals::T)::T where T
    - Takes a list of input values from all processors and returns an ordered
        contiguous list of those values on each processor

sumAll(comm::CommImpl, partialsums::T)::T where T
    - Take a list of input values from all processors and returns the sum on each
        processor.  The method +(::T, ::T)::T must exist

maxAll(comm::CommImpl, partialmaxes::T)::T where T
    - Takes a list of input values from all processors and returns the max to all
        processors.  The method <(::T, ::T)::Bool must exist

minAll(comm::CommImpl, partialmaxes::T)::T where T
    - Takes a list of input values from all processors and returns the min to all
        processors.  The method <(::T, ::T)::Bool must exist

scanSum(comm::CommImpl, myvals::T)::T where T
    - Takes a list of input values from all processors, computes the scan sum and
        returns it to all processors such that processor i contains the sum of
        values from processor 1 up to and including processor i.  The method
        +(::T, ::T)::T must exist

myPid(comm::CommImpl)::Integer - Returns the process rank

numProc(comm::CommImpl)::Integer - Returns the total number of processes

createDistributor(comm::CommImpl)::Distributor - Create a distributor object

createDirectory(comm::CommImpl, map::BlockMap)::Directory
    - Create a directory object for the given BlockMap

"""
abstract type Comm
end

function Base.show(io::IO, comm::Comm)
    println(typeof(comm)," with PID ", myPid(comm),
                " and ", numProc(comm), " processes")
end