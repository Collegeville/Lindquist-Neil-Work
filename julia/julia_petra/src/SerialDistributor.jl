
export SerialDistributor

type SerialDistributor <: Distributor
    post::Nullable{Array}
    reversePost::Nullable{Array}
    
    function SerialDistributor()
        new(nothing, nothing)
    end
end


function createFromSends(dist::SerialDistributor,
        exportPIDs::Array{PID})::Integer where PID <:Integer
    for id in exportPIDs
        @assert id == 1 "SerialDistributor can only accept PID of 1"
    end
    length(exportPIDs)
end

function createFromRecvs(
        dist::SerialDistributor, remoteGIDs::Array{GID}, remotePIDs::Array{PID}
        )::Tuple{Array{GID}, Array{PID}} where GID <: Integer where PID <: Integer
    @assert(length(remoteGIDs) == length(remotePIDs),
        "Number of GIDs and PIDs must be the same")
    for id in remotePIDs
        @assert id == 1 "SerialDistributor can only accept PID of 1"
    end
    remoteGIDs,remotePIDs
end

function resolve(dist::SerialDistributor, exportObjs::Array{T})::Array{T} where T
    exportObjs
end

function resolveReverse(dist::SerialDistributor, exportObjs::Array{T})::Array{T} where T
    exportObjs
end

function resolvePosts(dist::SerialDistributor, exportObjs::Array)
    dist.post = Nullable(exportObjs)
end

function resolveWaits(dist::SerialDistributor)::Array
    @assert !isnull(dist.post) "Must post before waiting"
    
    result = get(dist.post)
    dist.post = Nullable{Array}()
    result
end

function resolveReversePosts(dist::SerialDistributor, exportObjs::Array) 
    dist.reversePost = Nullable(exportObjs)
end

function resolveReverseWaits(dist::SerialDistributor)::Array
    @assert !isnull(dist.reversePost) "Must reverse post before reverse waiting"
    
    result = get(dist.reversePost)
    dist.reversePost = Nullable{Array}()
    result
end

