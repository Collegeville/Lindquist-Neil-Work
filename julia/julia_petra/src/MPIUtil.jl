#Contains MPI related things that the library is missing
import MPI

function MPI_RSend{T}(buf::MPI.MPIBuffertype{T}, count::Integer,
                dest::Integer, tag::Integer, comm::MPI.Comm)
    
    ccall(MPI.MPI_RSEND, Void,
        (Ptr{T}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
           Ptr{Cint}),
        buf, &count, &mpitype(T), &dest, &tag, &comm.val, &0)
end

function MPI_RSend{T}(buf::Array{T}, dest::Integer, tag::Integer, comm::MPI.Comm)
    MPI_RSend(buf, length(buf), dest, tag, comm)
end

function MPI_RSend{T}(obj::T, dest::Integer, tag::Integer, comm::MPI.Comm)
    MPI_RSend([obj], dest, tag, comm)
end