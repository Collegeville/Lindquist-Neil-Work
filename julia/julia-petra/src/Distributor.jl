

# methods (and docs) are currently based straight off Epetra_Distributor to match Comm

# TODO resolve any limitations on PID/GID/LID as decided
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
        removePIDs::Array{PID})::Tuple{Array{GID}, Array{PID}}
        where GID <: Integer where PID <: Integer
    - sets up the Distributor object using a list of remote global IDs and
        corresponding PIDs.  Returns a tuple with the global IDs and their
        respective processor IDs being sent to me.

do(dist::DistributorImpl, exportObjs::Array{T})::Array{T} where T
    - Execute the current plan on buffer of export objects and return the
        objects set to this processor

doReverse(dist::DistributorImpl, exportObjs::Array{T})::Array{T} where T
    - Execute the reverse of the current plan on buffer of export objects and
        return the objects set to this processor

doPosts(dist::DistributorImpl, exportObjs::Array)
    - Post buffer of export objects (can do other local work before executing
        Waits).  Otherwise, as do(::DistributorImpl, ::Array{T})::Array{T}

doWaits(dist::DistributorImpl)::Array - wait on a set of posts

doReversePosts(dist::DistributorImpl, exportObjs::Array)
    - Do reverse post of buffer of export objects (can do other local work
        before executing Waits).  Otherwise, as
        doReverse(::DistributorImpl, ::Array{T})::Array{T}

doReverseWaits(dist::DistributorImpl)::Array - wait on a reverse set of posts

"""
abstract type Distributor
end