
import MPI

"""
    MPIDistributor{GID, PID, LID}(comm::MPIComm{GID, PID, LID})
Creates an Distributor to work with MPIComm.  Created by
createDistributor(::MPIComm{GID, PID, LID})
"""
type MPIDistributor{GID <: Integer, PID <: Integer, LID <: Integer} <: Distributor{GID, PID, LID}
    comm::MPIComm
    
    lengths_to::Array{Integer}
    procs_to::Array{PID}
    indices_to::Array{Integer}
    
    lengths_from::Array{Integer}
    procs_from::Array{PID}
    indices_from::Array{Integer}
    
    resized::Bool
    sizes::Array{Integer}
    
    sizes_to::Array{Integer}
    starts_to::Array{Integer}
    #starts_to_ptr::Array{Integer}
    #indices_to_ptr::Array{Integer}
    
    sizes_from::Array{Integer}
    starts_from::Array{Integer}
    #sizes_from_ptr::Array{Integer}
    #starts_from_ptr::Array{Integer}
    
    numRecvs::Integer
    numSends::Integer
    numExports::Integer
    
    selfMsg::Integer
    
    maxSendLength::Integer
    totalRecvLength::Integer
    tag::Integer
    
    request::Array{MPI.Request}
    status::Array{MPI.Status}
    
    #sendArray::Array{UInt8}
    
    planReverse::Nullable{MPIDistributor}
    
    #never seem to be used
    #lastRoundBytesSend::Integer
    #lastRoundBytesRecv::Integer
    
    function MPIDistrbutor{GID <: Integer, PID <: Integer, LID <: Integer}(comm::MPIComm{GID, PID, LID})
        new{GID, PID, LID}(comm, [], [], [], [], [], [], false, [], [], [], [], [],
            0, 0, 0, 0, 0, 0,[], [],  Nullable{MPIDistributor}())
    end
end


#### internal methods ####
function createSendStructure(dist::MPIDistributor{GID, PID}, pid::PID, nProcs::PID, exportPIDs::Array{PID}) where GID <: Integer where PID <: Integer
    numExports = length(exportPIDs)
    dist.numExports = numExports
    
    #starts = Array{Integer}(nProcs + 1)
    #fill!(starts, 0)
    starts = zeros(nProcs+1)
    
    nactive = 0
    noSendBuff = true
    numDeadIndices = 0  #for GIDs not owned by any processors
    
    for i = 1:numExports
        if noSendBuff && i > 1 && exportPIDs[i] < exportPIDs[i-1]
            noSendBuff == false
        end
        if exportPID[i] >= 1
            starts[exportPIDs[i]] += 1
            nactive += 1
        else
            numDeadIndices += 1
        end
    end
    
    dist.selfMsg = starts[pid] != 1
    
    dist.nsends = 0
    
    if noSendBuff #grouped by processor, no send buffer or indicies_to needed
        for i = 1:nprocs
            if starts[i] > 0
                dist.nsends += 1
            end
        end
        
        dist.procs_to = Array{PID}(dist.nsends)
        dist.starts_to = Array{Integer}(dist.nsends)
        dist.lengths_to = Array{Integer}(dist.nsends)
        
        index = numDeadIndices
        for i = 1:dist.nsends
            dist.starts_to[i] = index
            proc = exportPIDs[index]
            dist.procs_to_[i] = proc
            index += dist.starts[proc]
        end
        
        perm = sortperm(dist.procs_to)
        dist.procs_to = dist.procs_to[perm]
        dist.starts_to = dist.starts_to[perm]
        
        # line 430
        
        dist.maxSendLength = 0
        
        for i = 1:dist.nsends
            proc = procs_to[i]
            lengths_to[i] = starts[proc]
            if (proc != pid) && (lengths_to[i] > dist.maxSendLength)
                maxSendLength = lengths_to[i]
            end
        end
    else #not grouped by processor, need send buffer and indices_to
        if starts[1] != 0
            dist.nsends = 1
        end
            
        for i = 2:nProcs
            if starts[i] != 0
                dist.nsends += 1
            end
            starts[i] += starts[i-1]
        end

        for i = nprocs:-1:1
            starts[i] = starts[i-1]
        end
        starts[1] = 0

        if nactive > 0
            dist.indices_to = Array{Integer}(nactive)
        end

        for i = 1:numExports
            if exportPIDs[i] >= 1
                dist.indices_to[starts[exportPIDs[i]]] = i
                starts[exportPIDs[i]] += 1
            end
        end

        #reconstruct starts array to index into indices_to

        for i = nprocs:-1:1
            starts[i] = starts[i-1]
        end
        starts[1] = 0
        starts[nProcs] = nactive

        if dist.nsends > 0
            dist.length_to = Array{Integer}(dist.nsends)
            dist.procs_to = Array{PID}(dist.nsends)
            dist.starts_to = Array{Integer}(dist.nsends)
        end

        j = 0
        dist.maxSendLength = 0
        for i = 1:nProcs
            if starts[i+1] != starts[i]
                dist.lengths_to[j] = starts[i+1] - starts[i]
                dist.starts_to[j] = starts[i]
                if (i != pid) && (dist.lengths_to[j] > dist.maxSendLength)
                    dist.maxSendLength = dist.lengths_to[j]
                end
                dist.procs_to[j] = i
                j += 1
            end
        end
    end
    
    dist.nsends -= dist.selfMsg
    
    dist
end

    
function computeRecvs(dist::MPIDistributor{GID, PID, LID}, myProc::PID, nProcs::PID) where GID <: Integer where PID <: Integer where LID <: Integer
    
    status = MPI.Status()
    
    msgCount = zeros(Int, nProcs)
    counts = ones(Int, nProcs)
    
    for i = 1:(dist.nsends + dist.selfMsg)
        msgCount[dist.procs_to[i]] = 1
    end
    
    
    #bug fix for reduce-scatter bug applied since no reduce_scatter is present in julia's MPI
    counts = MPI.reduce(msgCounts, MPI.SUM, 0, dist.comm.mpiComm)
    dist.nRecvs = MPI.scatter(counts, 1, 0, dist.comm.mpiComm)
    
    if dist.nRecvs > 0
        dist.lengths_from = zeros(Int, nProcs)
        dist.procs_from = zeros(PID, nProcs)
    end
    
    #using NEW_COMM_PATTERN (see line 590)
    
    if dist.nRecvs > 0
        if dist.request != []
            dist.request = Array{MPI.Request}(dist.nRecvs - dist.selfMsg)
            dist.status = Array{MPI.Status}(dist.nRecvs - dist.selfMsg)
        end
    end
    
    #line 616
    
    for i = 1:(dist.nRecvs - dist.selfMsg)
        dist.request[i] = MPI.Irecv!([dist.lengths_from[i]], MPI.ANY_SOURCE, dist.tag, dist.comm.MPIComm)
    end
    
    barrier(dist.comm)
    
    for i = 1:(dist.nsends+dist.selfMsg)
        if dist.procs_to[i] != myProc
            #have to use RSend in MPIUtil
            MPI_RSend(dist.lengths_to[i], dist.procs_to[i], dis.tag, dist.comm.MPIComm)
        else
            dist.lengths_from[dist.nRecvs] = dist.length_to[i]
            dist.procs_from[dist.nRecvs] = myProc
        end
    end
    
    if dist.nRecvs > dist.selfMsg
        dist.status = MPI.Waitall!(dist.request)
    end
    
    for i = 1:(dist.nRecvs - dist.selfMsg)
        dist.procs_from[i] = MPI.Get_source(dist.status[i])
    end
    
    perm = sortPerm(dist.procs_from)
    dist.procs_from = dist.procs_from[perm]
    dist.lengths_from = dist.lengths_from[perm]
    
    dist.starts_from = Array{Integer}(dist.nRecvs)
    j = 1
    for i = 1:dist.nRecvs
        dist.starts_from[i] = j
        j += dist.lengths_from[i]
    end
    
    dist.totalRecvLength = 0
    for i = 1:dist.nRecvs
        dist.totalRecvLength += dist.lengths_from[i]
    end
    
    dist.nRecvs -= dist.selfMsg
    
    barrier(dist.comm)
    
    dist        
end

#### Distributor interface ####

function createFromSends(dist::MPIDistributor{GID, PID, LID}, exportPIDs::Array{PID})::Integer where GID <:Integer where PID <:Integer where LID <:Integer
    pid = myPid(dist.comm)
    nProcs = numProcs(dist.comm)
    createSendStructure(dist, pid, nProcs, exportPIDs)
    computeRecvs(dist, pid, nprocs)
    
    if dist.numRecvs > 0
        if dist.request == []
            dist.request = Array{MPI.request}(numRecvs)
            dist.status = Array{MPI.status}(runRecvs)
        end
    end
    dist.totalRecvLength
end

"""
createFromRecvs(dist::MPIDistributor, remoteGIDs::Array{GID},
        removePIDs::Array{PID})::Tuple{Array{GID}, Array{PID}}
        where GID <: Integer where PID <: Integer
    - sets up the Distributor object using a list of remote global IDs and
        corresponding PIDs.  Returns a tuple with the global IDs and their
        respective processor IDs being sent to me.

resolvePosts(dist::MPIDistributor, exportObjs::Array)
    - Post buffer of export objects (can do other local work before executing
        Waits).  Otherwise, as do(::DistributorImpl, ::Array{T})::Array{T}

resolveWaits(dist::MPIDistributor)::Array - wait on a set of posts

resolveReversePosts(dist::MPIDistributor, exportObjs::Array)
    - Do reverse post of buffer of export objects (can do other local work
        before executing Waits).  Otherwise, as
        doReverse(::DistributorImpl, ::Array{T})::Array{T}

resolveReverseWaits(dist::MPIDistributor)::Array - wait on a reverse set of posts
"""
function asdf()#dummy function to allow me to "block comment" the other stubs
end