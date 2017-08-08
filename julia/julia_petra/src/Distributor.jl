

export Distributor
export createFromSends, createFromRecvs
export resolve, resolveReverse, resolveWaits
export resolvePosts, resolveReversePosts, resolveReverseWaits

# methods (and docs) are currently based straight off Epetra_Distributor to match Comm

"""
The base type for gather/scatter setup.
All subtypes must have the following methods, with DistributorImpl standing
in for the subtype:

createFromSends(dist::DistributorImpl,exportPIDs::Array{PID})
        ::Integer where PID <:Integer
    - sets up the Distributor object using a list of process IDs to which we
        export and the number of IDs being exported.  Returns the number of
        IDs this processor will be reciving

createFromRecvs(dist::DistributorImpl, remoteGIDs::Array{GID},
        remotePIDs::Array{PID})::Tuple{Array{GID}, Array{PID}}
        where GID <: Integer where PID <: Integer
    - sets up the Distributor object using a list of remote global IDs and
        corresponding PIDs.  Returns a tuple with the global IDs and their
        respective processor IDs being sent by me.

resolvePosts(dist::DistributorImpl, exportObjs::Array)
    - Post buffer of export objects (can do other local work before executing
        Waits).  Otherwise, as do(::DistributorImpl, ::Array{T})::Array{T}

resolveWaits(dist::DistributorImpl)::Array - wait on a set of posts

resolveReversePosts(dist::DistributorImpl, exportObjs::Array)
    - Do reverse post of buffer of export objects (can do other local work
        before executing Waits).  Otherwise, as
        doReverse(::DistributorImpl, ::Array{T})::Array{T}

resolveReverseWaits(dist::DistributorImpl)::Array - wait on a reverse set of posts

"""
abstract type Distributor{GID <: Integer, PID <: Integer, LID <: Integer}
end

"""
    createFromSends(dist::Distributor, exportPIDs::Array{<:Integer})::Integer

Sets up the Distributor object using a list of process IDs to which we
export and the number of IDs being exported.  Returns the number of
IDs this processor will be receiving
"""
function createFromSends(dist::Distributor{GID, PID, LID}, exportPIDs::Array{<:Integer}) where GID <: Integer where PID <: Integer where LID <: Integer
    createFromSends(dist, Array{PID}(exportPIDs))
end

"""
    createFromRecvs(dist::Distributor{GID, PID, LID}, remoteGIDs::Array{<:Integer}, remotePIDs::Array{<:Integer})::Tuple{Array{GID}, Array{PID}}

Sets up the Distributor object using a list of remote global IDs and
corresponding PIDs.  Returns a tuple with the global IDs and their
respective processor IDs being sent by me.
"""
function createFromRecvs(dist::Distributor{GID, PID, LID}, remoteGIDs::Array{<:Integer},
        remotePIDs::Array{<:Integer}) where GID <: Integer where PID <: Integer where LID <: Integer
    createFromRecvs(dist, Array{GID}(remoteGIDs), Array{PID}(remotePIDs))
end

"""
    resolve(dist::Distributor, exportObjs::Array{T})::Array{T}

Execute the current plan on buffer of export objects and return the
objects set to this processor
"""
function resolve(dist::Distributor, exportObjs::Array{T})::Array{T} where T
    resolvePosts(dist, exportObjs)
    resolveWaits(dist)
end

"""
    resolveReverse(dist::Distributor, exportObjs::Array{T})::Array{T}

Execute the reverse of the current plan on buffer of export objects and
return the objects set to this processor
"""
function resolveReverse(dist::Distributor, exportObjs::Array{T})::Array{T} where T
    resolveReversePosts(dist, exportObjs)
    resolveReverseWaits(dist)
end


#### required method documentation stubs ####

"""
    resolvePosts(dist::Distributor, exportObjs::Array)

Post buffer of export objects (can do other local work before executing
Waits.  Otherwise, as resolve(::DistributorImpl, ::Array{T})::Array{T}
"""
function resolvePosts end

"""
    resolveWaits(dist::Distributor)::Array

wait on a set of posts
"""
function resolveWaits end

"""
    resolveReversePosts(dist::Distributor, exportObjs::Array)
Do reverse post of buffer of export objects (can do other local work
before executing Waits).  Otherwise, as
doReverse(::DistributorImpl, ::Array{T})::Array{T}
"""
function resolveReversePosts end

"""
    resolveReverseWaits(dist::Distributor)::Array

wait on a reverse set of posts
"""
function resolveReverseWaits end
