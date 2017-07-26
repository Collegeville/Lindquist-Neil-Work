
export SerialDistributor

type SerialDistributor{GID <: Integer, PID <:Integer, LID <: Integer} <: Distributor{GID, PID, LID}
    post::Nullable{Array}
    reversePost::Nullable{Array}
    
    function SerialDistributor{GID, PID, LID}() where GID <: Integer where PID <: Integer where LID <: Integer
        new(nothing, nothing)
    end
end


function createFromSends(dist::SerialDistributor{GID, PID, LID},
        exportPIDs::Array{PID})::Integer where GID <: Integer where PID <: Integer where LID <: Integer
    for id in exportPIDs
        if id != 1
            throw(InvalidArgumentError("SerialDistributor can only accept PID of 1"))
        end
    end
    length(exportPIDs)
end

function createFromRecvs(
        dist::SerialDistributor{GID, PID, LID}, remoteGIDs::Array{GID}, remotePIDs::Array{PID}
        )::Tuple{Array{GID}, Array{PID}} where GID <: Integer where PID <: Integer where LID <: Integer
    if length(remoteGIDs) != length(remotePIDs)
        throw(InvalidArgumentError("Number of GIDs and PIDs must be the same"))
    end
    for id in remotePIDs
        if id != 1
            throw(InvalidArgumentError("SerialDistributor can only accept PID of 1"))
        end
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
    if isnull(dist.post)
        throw(InvalidStateException("Must post before waiting", :SerialDistributor))
    end
    
    result = get(dist.post)
    dist.post = Nullable{Array}()
    result
end

function resolveReversePosts(dist::SerialDistributor, exportObjs::Array) 
    dist.reversePost = Nullable(exportObjs)
end

function resolveReverseWaits(dist::SerialDistributor)::Array
     if isnull(dist.reversePost)
        throw(InvalidStateException("Must reverse post before reverse waiting",
                :SerialDistributor))
    end
    
    result = get(dist.reversePost)
    dist.reversePost = Nullable{Array}()
    result
end

