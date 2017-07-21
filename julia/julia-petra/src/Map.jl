
# methods and docs based straight off Epetra_Map to match Comm

# ignoring indexBase methods
# ignoring elementSize methods since, type information is carried anyways

# This class is the abstract parent for julia implmentations of Epetra_BlockMap,
    # Epetra_Map and Epetra_LocalMap

# TODO figure out firstPointInElement and related
# TODO figure out expert users and developers only functions
"""
The base class for partitioning block element vectors and matrices.
All subtypes must have the following methods, with MapImpl standing in for the subtype:

--local/global ID accessor methods--

removeIDList(map::MapImpl, gidList::Array{GID})::Tuple{Array{PID}, Array{LID}}
        where GID <: Integer where PID <: Integer where LID <: Integer
    - Return the processor ID and local index value for a given list of global indices

lid(map::MapImpl, gid::GID)::LID where GID <: Integer where LID <: Integer
    - Return local ID of global ID, or 0 if not found on this processor

gid(map::MapImpl, lid::LID)::GID where GID <: Integer where LID <: Integer
    - Return global ID of local ID, or 0 if not found on this processor

findLocalElementID(map::MapImpl, pointID::Integer)::Tuple{LID, OFF}
        where LID <: Integer where OFF <: Integer
    - Return the LID of the element that contains the given local point ID
        and the offset of the point in that element

minAllGID(map::MapImpl)::GID where GID <: Integer
    - Return the minimum global ID across the entire map

maxAllGID(map::MapImpl)::GID where GID <: Integer
    - Return the maximum global ID across the entire map

minLID(map::MapImpl)::LID where LID <: Integer
    - Return the mimimum local index value on the calling processor

maxLID(map::MapImpl)::LID where LID <: Integer
    - Return the maximum local index value on the calling processor


--size/dimension accessor functions--

numGlobalElements(map::MapImpl)::GID where GID <: Integer
    - Return the number of elements across all processors

myGlobalElements(map{T}::MapImpl)::Array{T}
    - Return a list of global elements on this processor

numGlobalPoints(map::MapImpl)::GID where GID <: Integer
    - Return the number of global points for this map

numMyPoints(map::MapImpl)::LID where LID <: Integer
    - Return the number of local points for this map


--Miscellaneous boolean tests--

uniqueGIDs(map::MapImpl)::Bool
    - Return true if each map GID exists on at most 1 processor

globalIndicesType(map::MapImpl)::Type{GID} where GID <: Integer
    - Return the type used for global indices in the map

sameBlockMapDataAs(this::MapImpl, other::Map)::Bool
    - Return true if the maps have the same data 

sameAs(this::MapImpl, other::Map)::Bool
    - Return true if this and other are identical maps

pointSameAs(this::MapImpl, other::Map)::Bool
    - Return true if this and other have identical point-wise structure

linearMap(map::MapImpl)::Bool
    - Return true if the global ID space is contiguously divided (but
        not necessarily uniformly) across all processors


--Array accessor functions--

myGlobalElementIDs(map::MapImpl)::Array{GID} where GID <: Integer
    - Return list of global IDs assigned to the calling processor

pointToElementList(map::MapImpl)::Array{LID} where LID <: Integer
    - Return list where for each local point, indicates the local
        element ID that the point belongs to

Comm(map::MapImpl)::Comm - Return the Comm for the map
"""
abstract type Map{T} where {T}
end

"""
Return true if the GID passed in belongs to the calling processor in this
map, otherwise returns false.
"""
function myGID(map::Map, gid::GID) where GID <: Integer
    return GID(map, gid) != 0
end

"""
Return true if the LID passed in belongs to the calling processor in this
map, otherwise returns false.
"""
function myLID(map::Map, lid::LID) where LID <: Integer
    return LID(map, lid) != 0
end

"""
Return true if map is defined across more than one processor
"""
function distributedGlobal(map::Map)
    return numGlobalElements(map) != numLocalElements(map)
end

"""
Return the number of elements across the calling processor
"""
function numMyElements(map::Map)::Integer #::LID where LID <: Integer
    length(myGlobalElementIDs(map))
end

"""
Return the minimum global ID owned by this processor
"""
function minMyGID(map::Map)::Integer #::GID where GID <: Integer
    minimum(myGlobalElementsIDs(map))
end
    
"""
Return the maximum global ID owned by this processor
"""
function maxMyGID(map::Map)::Integer #::GID where GID <: Integer
    maximum(myGlobalElementsIDs(map))
end