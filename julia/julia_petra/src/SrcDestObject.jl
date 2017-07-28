
"""
A base type for supporting flexible source distributed objects for import/export operations.

Subtypes must implement a map(::Impl{GID, PID, LID})::BlockMap{GID, PID, LID} method,
where Impl is the subtype
"""
abstract type SrcDestObject{GID <: Integer, PID <: Integer, LID <: Integer}
end

"""
Returns true if this object is a distributed global
"""
function distributedGlobal(obj::SrcDestObject)
    distributedGlobal(map(obj))
end


"""
Get's the Comm instance being used by this object
"""
function comm(obj::SrcDestObject)
    comm(map(obj))
end