
import MPI

export MPIDistributor

"""
    MPIDistributor{GID, PID, LID}(comm::MPIComm{GID, PID, LID})
Creates an Distributor to work with MPIComm.  Created by
createDistributor(::MPIComm{GID, PID, LID})
"""
mutable struct MPIDistributor{GID <: Integer, PID <: Integer, LID <: Integer} <: Distributor{GID, PID, LID}
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
    
    importObjs::Nullable{Array{Array{UInt8}}}
    
    #never seem to be used
    #lastRoundBytesSend::Integer
    #lastRoundBytesRecv::Integer
    
    function MPIDistributor(comm::MPIComm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
        new{GID, PID, LID}(comm, [], [], [], [], [], [], false, [], [], [], [], [],
            0, 0, 0, 0, 0, 0, 0, [], [],  Nullable{MPIDistributor}(),
            Nullable{Array{UInt8}}())
    end
end


#### internal methods ####
function createSendStructure(dist::MPIDistributor{GID, PID}, pid::PID, nProcs::PID, exportPIDs::Array{PID}) where GID <: Integer where PID <: Integer
    
    numExports = length(exportPIDs)
    dist.numExports = numExports
    
    #starts = Array{Integer}(nProcs + 1)
    #fill!(starts, 0)
    starts = zeros(GID, nProcs+1)
    
    nactive = 0
    noSendBuff = true
    numDeadIndices = 0  #for GIDs not owned by any processors
    
    for i = 1:numExports
        if noSendBuff && i > 1 && exportPIDs[i] < exportPIDs[i-1]
            noSendBuff = false
        end
        if exportPIDs[i] >= 1
            starts[exportPIDs[i]] += 1
            nactive += 1
        else
            numDeadIndices += 1
        end
    end
    
    dist.selfMsg = starts[pid] != 0
    dist.numSends = 0
    
    if noSendBuff #grouped by processor, no send buffer or indices_to needed
        for i = 1:nProcs
            if starts[i] > 0
                dist.numSends += 1
            end
        end
        
        dist.procs_to = Array{PID}(dist.numSends)
        dist.starts_to = Array{Integer}(dist.numSends)
        dist.lengths_to = Array{Integer}(dist.numSends)
        
        index = numDeadIndices+1
        for i = 1:dist.numSends
            dist.starts_to[i] = index
            proc = exportPIDs[index]
            dist.procs_to[i] = proc
            index += starts[proc]
        end
        
        perm = sortperm(dist.procs_to)
        dist.procs_to = dist.procs_to[perm]
        dist.starts_to = dist.starts_to[perm]
        
        # line 430
        
        dist.maxSendLength = 0
        
        for i = 1:dist.numSends
            proc = dist.procs_to[i]
            dist.lengths_to[i] = starts[proc]
            if (proc != pid) && (dist.lengths_to[i] > dist.maxSendLength)
                maxSendLength = dist.lengths_to[i]
            end
        end
    else #not grouped by processor, need send buffer and indices_to
        if starts[1] != 0
            dist.numSends = 1
        end
            
        for i = 2:nProcs
            if starts[i] != 0
                dist.numSends += 1
            end
            starts[i] += starts[i-1]
        end

        for i = nProcs:-1:2
            starts[i] = starts[i-1] + 1
        end
        starts[1] = 1

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

        for i = nProcs:-1:2
            starts[i] = starts[i-1]
        end
        starts[1] = 1
        starts[nProcs+1] = nactive+1

        
        if dist.numSends > 0
            dist.lengths_to = Array{Integer}(dist.numSends)
            dist.procs_to = Array{PID}(dist.numSends)
            dist.starts_to = Array{Integer}(dist.numSends)
        end

        j = 1
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
    
    dist.numSends -= dist.selfMsg
    
    dist
end

    
function computeRecvs(dist::MPIDistributor{GID, PID, LID}, myProc::PID, nProcs::PID) where GID <: Integer where PID <: Integer where LID <: Integer
    
    status = MPI.Status()
    
    msgCount = zeros(Int, nProcs)
    counts = ones(Int, nProcs)
    
    for i = 1:(dist.numSends + dist.selfMsg)
        msgCount[dist.procs_to[i]] += 1
    end
    
    #bug fix for reduce-scatter bug applied since no reduce_scatter is present in julia's MPI
    counts = MPI.Reduce(msgCount, MPI.SUM, 0, dist.comm.mpiComm)
    if counts == nothing
        counts = Int[]
    end
    dist.numRecvs = MPI.Scatter(counts, 1, 0, dist.comm.mpiComm)[1]
    
    dist.lengths_from = zeros(Int, nProcs)
    dist.procs_from = zeros(PID, nProcs)
    
    #using NEW_COMM_PATTERN (see line 590)
    
    if dist.request == []
        dist.request = Array{MPI.Request}(dist.numRecvs - dist.selfMsg)
    end
    
    #line 616
    
    lengthWrappers = [Array{Int, 1}(1) for i in 1:(dist.numRecvs - dist.selfMsg)]
    for i = 1:(dist.numRecvs - dist.selfMsg)
        dist.request[i] = MPI.Irecv!(lengthWrappers[i], MPI.ANY_SOURCE, dist.tag, dist.comm.mpiComm)
    end
    
    for i = 1:(dist.numSends+dist.selfMsg)
        if dist.procs_to[i] != myProc
            #have to use Rsend in MPIUtil
            MPI_Rsend(dist.lengths_to[i], dist.procs_to[i]-1, dist.tag, dist.comm.mpiComm)
        else
            dist.lengths_from[dist.numRecvs] = dist.lengths_to[i]
            dist.procs_from[dist.numRecvs] = myProc
        end
    end
    
    if dist.numRecvs > dist.selfMsg
        dist.status = MPI.Waitall!(dist.request)
    end
    
    for i = 1:(dist.numRecvs - dist.selfMsg)
        dist.lengths_from[i] = lengthWrappers[i][1]
    end
    
    
    for i = 1:(dist.numRecvs - dist.selfMsg)
        dist.procs_from[i] = MPI.Get_source(dist.status[i])+1
    end
    
    perm = sortperm(dist.procs_from)
    dist.procs_from = dist.procs_from[perm]
    dist.lengths_from = dist.lengths_from[perm]
    
    dist.starts_from = Array{Integer}(dist.numRecvs)
    j = 1
    for i = 1:dist.numRecvs
        dist.starts_from[i] = j
        j += dist.lengths_from[i]
    end
    
    barrier(dist.comm) #testing
    dist.totalRecvLength = 0
    for i = 1:dist.numRecvs
        dist.totalRecvLength += dist.lengths_from[i]
    end
    
    dist.numRecvs -= dist.selfMsg
    
    barrier(dist.comm)
    
    dist        
end

function computeSends(dist::MPIDistributor{GID, PID, LID}, remoteGIDs::Array{GID}, remotePIDs::Array{PID})::Tuple{Array{GID}, Array{PID}} where GID <:Integer where PID <:Integer where LID <:Integer
    numImports = length(remoteGIDs)

    tmpPlan = MPIDistributor(dist.comm)
    
    procList = copy(remotePIDs)
    importObjs = Array{Tuple{GID, PID}}(numImports)
    for i = 1:numImports
        importObjs[i] = (remoteGIDs[i], myPid(dist.comm))#remotePIDs[i])
    end
    
    numExports = createFromSends(tmpPlan, procList)
    
    exportIDs = Array{GID}(numExports)
    exportProcs = Array{PID}(numExports)
    
    exportObjs = resolve(tmpPlan, importObjs)
    for i = 1:numExports
        exportIDs[i] = exportObjs[i][1]
        exportProcs[i] = exportObjs[i][2]
    end
    (exportIDs, exportProcs)
end

"""
Creates a reverse distributor for the given MPIDistributor
"""
function createReverseDistributor(dist::MPIDistributor)
    myProc = myPid(dist.comm)
    
    if isnull(dist.planReverse)
        totalSendLength = reduce(+, dist.lengths_to)
        
        maxRecvLength = 0
        for i = 1:dist.numRecvs
            if dist.procs_from[i] != myProc
                maxRecvLength = max(maxRecvLength, dist.lengths_from[i])
            end
        end
        
        reverse = MPIDistributor(dist.comm)
        dist.planReverse = Nullable(reverse)
        
        reverse.lengths_to = dist.lengths_from
        reverse.procs_to = dist.procs_from
        reverse.indices_to = dist.indices_from
        reverse.starts_to = dist.starts_from
        
        reverse.lengths_from = dist.lengths_to
        reverse.procs_from = dist.procs_to
        reverse.indices_from = dist.indices_to
        reverse.starts_from = dist.starts_to
        
        reverse.numSends = dist.numRecvs
        reverse.numRecvs = dist.numSends
        reverse.selfMsg = dist.selfMsg
        
        reverse.maxSendLength = maxRecvLength
        reverse.totalRecvLength = totalSendLength
        
        reverse.request = Array{MPI.Request}(reverse.numRecvs)
        reverse.status  = Array{MPI.Status}(reverse.numRecvs)
    end
end


#### Distributor interface ####

function createFromSends(dist::MPIDistributor{GID, PID, LID}, exportPIDs::Array{PID})::Integer where GID <:Integer where PID <:Integer where LID <:Integer
    pid = myPid(dist.comm)
    nProcs = numProc(dist.comm)
    createSendStructure(dist, pid, nProcs, exportPIDs)
    computeRecvs(dist, pid, nProcs)
    if dist.numRecvs > 0
        if dist.request == []
            dist.request = Array{MPI.request}(numRecvs)
            dist.status = Array{MPI.status}(runRecvs)
        end
    end
    dist.totalRecvLength
end


function createFromRecvs(dist::MPIDistributor{GID, PID}, remoteGIDs::Array{GID}, remotePIDs::Array{PID})::Tuple{Array{GID}, Array{PID}} where GID <: Integer where PID <: Integer
    if length(remoteGIDs) != length(remotePIDs)
        throw(InvalidArgumentError("remote lists must be the same length"))
    end
    (exportGIDs, exportPIDs) = computeSends(dist, remoteGIDs, remotePIDs)
    
    createFromSends(dist, exportPIDs)
    
    exportGIDs, exportPIDs
end

function resolvePosts(dist::MPIDistributor, exportObjs::Array{T}) where T
    myProc = myPid(dist.comm) 
    
    selfRecvAddress = 0
    
    nBlocks = dist.numSends + dist.selfMsg
    procIndex = 1
    while procIndex <= nBlocks && dist.procs_to[procIndex] < myProc
        procIndex += 1
    end
    if procIndex == nBlocks
        procIndex = 1
    end
    
    
    exportBytes = Array{Array{UInt8}}(dist.numRecvs + dist.selfMsg)
    
    buffer = IOBuffer()
    if dist.indices_to == [] #data already grouped by processor
        j = 1
        for i = 1:nBlocks
            Serializer.serialize(buffer, exportObjs[j:j+dist.lengths_to[i]-1])
            j += dist.lengths_to[i]
            exportBytes[i] = take!(buffer)
        end
    else #data not grouped by proc, must be grouped first
        for i = 1:nBlocks
            j = dist.starts_to[i]
            sendArray = Array{T}(dist.lengths_to[i])
            for k = 1:dist.lengths_to[i]
                sendArray[k] = exportObjs[dist.indices_to[j+k-1]]
            end
            Serializer.serialize(buffer, sendArray)
            exportBytes[i] = take!(buffer)
        end
    end
    
    
    ## get sizes of data begin received ##
    lengthRequests = Array{MPI.Request}(dist.numRecvs)
    lengths = Array{Array{Int, 1}}(dist.numRecvs + dist.selfMsg)
    for i = 1:(dist.numRecvs + dist.selfMsg)
        lengths[i] = Array{Int}(1)
    end
    
    j = 1
    for i = 1:(dist.numRecvs + dist.selfMsg)
        if dist.procs_from[i] != myProc
            lengthRequests[j] = MPI.Irecv!(lengths[i], dist.procs_from[i]-1, dist.tag, dist.comm.mpiComm)
            j += 1
        end
    end
    
    barrier(dist.comm)
           
    for i = 1:(dist.numSends + dist.selfMsg)
        p = procIndex + 1
        if p > nBlocks
            p -= nBlocks
        end
        if dist.procs_to[i] != myProc
            MPI_Rsend(length(exportBytes[i]), dist.procs_to[i]-1, dist.tag, dist.comm.mpiComm)
        else
            lengths[i][1] = length(exportBytes[i])
        end
    end
    
    MPI.Waitall!(lengthRequests)
    #at this point `lengths` should contain the sizes of incoming data
    
    importObjs = Array{Array{UInt8}}(dist.numRecvs+dist.selfMsg)
    for i = 1:length(importObjs)
        importObjs[i] = Array{UInt8}(lengths[i][1])
    end
    
    dist.importObjs = Nullable(importObjs)
    
    ## back to the regularly scheduled program ##
    
    k = 0
    j = 0
    for i = 1:(dist.numRecvs + dist.selfMsg)
        if dist.procs_from[i] != myProc
            MPI.Irecv!(importObjs[i], dist.procs_from[i]-1, dist.tag, dist.comm.mpiComm)
        else
            selfRecvAddress = i
        end
    end
    
    barrier(dist.comm)
    
    #line 844
    
    selfNum = 0
    
#    if dist.indices_to == [] #data already grouped by processor
    for i = 1:nBlocks
        p = i + procIndex - 1
        if p > nBlocks
            p -= nBlocks
        end
        if dist.procs_to[p] != myProc
            MPI_Rsend(exportBytes[p], dist.procs_to[p]-1, dist.tag, dist.comm.mpiComm)
        else
            selfNum = p
        end
    end

    if dist.selfMsg
        importObjs[selfRecvAddress] = exportBytes[selfNum]
    end
end

            

function resolveWaits(dist::MPIDistributor)::Array
    barrier(dist.comm)#run into issues deserializing otherwise
    if dist.numRecvs> 0
        dist.status = MPI.Waitall!(dist.request)
    end
    
    if isnull(dist.importObjs)
        throw(InvalidStateError("Cannot resolve waits when no posts have been made"))
    end
    
    importObjs = get(dist.importObjs)
    deserializedObjs = Array{Array{Any}}(length(importObjs))
    for i = 1:length(importObjs)
        deserializedObjs[i] = Serializer.deserialize(IOBuffer(importObjs[i]))
    end
    
    dist.importObjs = Nullable{Array{Array{UInt8}}}()
    
    reduce(vcat, [], deserializedObjs)
end


function resolveReversePosts(dist::MPIDistributor, exportObjs::Array)
    if dist.indices_to != []
        throw(InvalidStateError("Cannot do reverse comm when data is not blocked by processor"))
    end
    
    if isnull(dist.planReverse)
        createReverseDistributor(dist)
    end
    
    resolvePosts(get(dist.planReverse), exportObjs)
end


function resolveReverseWaits(dist::MPIDistributor)::Array
    if isnull(dist.planReverse)
        throw(InvalidStateError("Cannot resolve reverse waits if there is no reverse plan"))
    end
    
    resolveWaits(get(dist.planReverse))
end
