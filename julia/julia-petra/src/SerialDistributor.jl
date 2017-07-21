
type SerialDistributor <: Distributor
    post::Array
    reversePost::Array
    
    function SerialDistributor()
        this(nothing, nothing)
    end
end


function createFromSends(dist::SerialDistributor,
        exportPIDs::Array{PID})::Integer where PID <:Integer
    for id in exportPIDs
        if id != 1
            error("SerialDistributor can only accept PID of 1")
        end
    end
    numExportIDs
end

function createFromRecvs(dist::SerialDistributor, remoteGIDs::Array{GID},
        remotePIDs::Array{PID})::Tuple{Array{GID}, Array{PID}}
        where GID <: Integer where PID <: Integer
    
    # DECISION: check if GIDs exist?
    for id in remotePIDs
        if id != 1
            error("SerialDistributor can only accept PID of 1")
        end
    end
    remoteGIDs,removePIDs
end

function resolve(dist::SerialDistributor, exportObjs::Array{T})::Array{T} where T
    exportObjs
end

function resolveReverse(dist::SerialDistributor, exportObjs::Array{T})::Array{T} where T
    exportObjs
end

function resolvePosts(dist::SerialDistributor, exportObjs::Array)
    dist.post = exportObjs
end

function resolveWaits(dist::SerialDistributor)::Array
    if dist.post == nothing
        error("Must post before waiting")
    end
    result = dist.post
    dist.post = nothing
    result
end

function resolveReversePosts(dist::SerialDistributor, exportObjs::Array) 
    dist.reversePost = exportObjs
end

function resolveReverseWaits(dist::SerialDistributor)::Array
    if dist.reversePost == nothing
        error("Must reverse post before reverse waiting")
    end
    result = dist.reversePost
    dist.reversePost = nothing
    result
end

