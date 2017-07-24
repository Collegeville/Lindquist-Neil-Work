
export BlockMap
export remoteIDList, lid, gid, findLocalElementID
export minAllGID, maxAllGID, minMyGID, maxMyGID, minLID, maxLID
export numGlobalElements, myGlobalElements, numGlobalPoints, numMyPoints
export uniqueGIDs, globalIndicesType, sameBlockMapDataAs, sameAs, pointSameAs
export linearMap, myGlobalElementIDs, pointToElementList, Comm
export myGID, myLID, distributedGlobal, numMyElements

# TODO implement BlockMap


# methods and docs based straight off Epetra_BlockMap to match Comm

# ignoring indexBase methods and sticking with 1-based indexing
# ignoring elementSize methods since, type information is carried anyways


# TODO figure out firstPointInElement and related
# TODO figure out expert users and developers only functions


"""
A type for partitioning block element vectors and matrices
"""
type BlockMap{T}
    data::BlockMapData
end

"""
Return true if the GID passed in belongs to the calling processor in this
map, otherwise returns false.
"""
function myGID(map::BlockMap, gid::GID) where GID <: Integer
    GID(map, gid) != 0
end

"""
Return true if the LID passed in belongs to the calling processor in this
map, otherwise returns false.
"""
function myLID(map::BlockMap, lid::LID) where LID <: Integer
    LID(map, lid) != 0
end

"""
Return true if map is defined across more than one processor
"""
function distributedGlobal(map::BlockMap)
    #numGlobalElements(map) != numLocalElements(map)
    map.data.distributedGlobal
end

"""
Return the number of elements across the calling processor
"""
function numMyElements(map::BlockMap)::Integer #::LID where LID <: Integer
    #length(myGlobalElementIDs(map))
    map.data.numMyElements
end

"""
Return the minimum global ID owned by this processor
"""
function minMyGID(map::BlockMap)::Integer #::GID where GID <: Integer
    #minimum(myGlobalElementsIDs(map))
    map.data.minMyGID
end
    
"""
Return the maximum global ID owned by this processor
"""
function maxMyGID(map::BlockMap)::Integer #::GID where GID <: Integer
    #maximum(myGlobalElementsIDs(map))
    map.data.maxMyGID
end


##local/global ID accessor methods##

"""
Return the processor ID and local index value for a given list of global indices.
The returned value is a tuple containing
    1 - an Array of processors owning the global ID's in question
    2 - an Array of local IDs of the global on the owning processor
"""
function remoteIDList(map::BlockMap, gidList::Array{GID}
        )::Tuple{Array{Integer}, Array{Integer}} where GID <: Integer
        #where PID <: Integer where LID <: Integer
    data = map.data
    if isnull(data.directory)
        data.directory = createDirectory(data.comm, map)
    end
    
    getDirectoryEntries(map, gidList)
end

"""
Return local ID of global ID, or 0 if not found on this processor
"""
function lid(map::BlockMap, gid::GID)::Integer where GID <: Integer
        #where LID <: Integer
    data = map.data
    if (gid < data.minMyGID) || (gid > data.maxMyGID)
        return 0
    end
    if data.linearMap
        return gid-data.minMyGID
    end
    if gid >= data.myGlobalElements[1] && gid <= data.lastContiguousGID
        return gid - data.myGlobalElements[1]
    end
    
    return data.lidHash[gid]
end

"""
Return global ID of local ID, or 0 if not found on this processor
"""
function gid(map::BlockMap, lid::LID)::Integer where LID <: Integer 
        #where GID <: Integer 
    data = map.data
    if (data.numMyElements == 0) || (lid < data.minLID) || (lid > data.maxLID)
        return 0
    end
    
    if data.linearMap
        return lid + data.minMyGID
    end
    return data.myGlobalElements[lid]
end


#TODO figure this out, it has to do with element size
#"""
#Return the LID of the element that contains the given local point ID
#and the offset of the point in that element
#"""
#function findLocalElementID(map::BlockMap, pointID::Integer)::Tuple{LID, OFF}
#        where LID <: Integer where OFF <: Integer
#end

"""
Return the minimum global ID across the entire map
"""
function minAllGID(map::BlockMap)::Integer #::GID where GID <: Integer
    map.data.minAllGID
end

"""
Return the maximum global ID across the entire map
"""
function maxAllGID(map::BlockMap)::Integer
        #::GID where GID <: Integer
    map.data.maxAllGID
end

"""
Return the mimimum local index value on the calling processor
"""
function minLID(map::BlockMap)::Integer #::LID where LID <: Integer
    map.data.minLID
end

"""
Return the maximum local index value on the calling processor
"""
function maxLID(map::BlockMap)::Integer#::LID where LID <: Integer
    map.data.maxLID
end

##size/dimension accessor functions##

"""
Return the number of elements across all processors
"""
function numGlobalElements(map::BlockMap)::Integer
        #::GID where GID <: Integer
    map.data.numGlobalElements
end

"""
Return a list of global elements on this processor
"""
function myGlobalElements(map::BlockMap)::Array{Integer}
    map.data.myGlobalElements
end

"""
Return the number of global points for this map
"""
function numGlobalPoints(map::BlockMap)::Integer
        #::GID where GID <: Integer
    map.data.numGlobalPoints
end

"""
Return the number of local points for this map
"""
function numMyPoints(map::BlockMap)::Integer
        #::LID where LID <: Integer
    map.data.numMyPoints
end


##Miscellaneous boolean tests##

"""
Return true if each map GID exists on at most 1 processor
"""
function uniqueGIDs(map::BlockMap)::Bool
    isOneToOne(map)
end

#DECISION figure out if this is nessacery, since I've been using generic Integer
#function globalIndicesType(map::BlockMap)::Type{GID} where GID <: Integer
#    - Return the type used for global indices in the map

"""
Return true if the maps have the same data
"""
function sameBlockMapDataAs(this::BlockMap, other::BlockMap)::Bool
    this.data == other.data
end

"""
Return true if this and other are identical maps
"""
function sameAs(this::BlockMap, other::BlockMap)::Bool
    tData = this.data
    oData = other.data
    if tData == oData
        return true
    end
    
    if ((tData.minAllGID != oData.minAllGID)
        || (tData.maxAllGID != oData.maxAllGID)
        || (tData.numGlobalElements != oData.numGlobalElements))
        return false
    end
        
    mySameMap = 1
    
    #TODO add checks for element size if added
    
    if tData.numMyElements != oData.numMyElements
        mySameMap = 0
    end
    
    if tData.linearMap && oData.linearMap
        # For linear maps, just need to check whether lower bound is the same
        if tData.minMyGID != oData.minMyGID
            mySameMap = 0
        end
    else
        for i = 1:tData.numMyElements
            if gid(this, i) == gid(other, i)
                mySameMap = 0
                break
            end
        end
    end
    
    Bool(minAll(tData.comm, [mySameMap])[1])
end
    
"""
Return true if this and other have identical point-wise structure
"""
function pointSameAs(this::BlockMap, other::BlockMap)::Bool
    tData = this.data
    oData = this.data
    if tData == oData
        return true
    end
    
    if tData.numGlobalPoints != oData.numGlobalPoints
        return false
    end
    
    mySameMap = 1
    if tData.numMyPoints != oData.numMyPoints
        mySameMap = 0
    end
    
    Bool(minAll(tData.comm, [mySameMap])[1])
end


"""
Return true if the global ID space is contiguously divided (but
not necessarily uniformly) across all processors
"""
function linearMap(map::BlockMap)::Bool
    map.data.linearMap
end


##Array accessor functions##

"""
Return list of global IDs assigned to the calling processor
"""
function myGlobalElementIDs(map::BlockMap)::Array{Integer}
    #::Array{GID} where GID <: Integer
    data = map.data
    myGlobalElements = Array{Integer}(data.numMyElements)
    if length(data.myGlobalElements) == 0
        for i = 1:data.numMyElements
            myGlobalElements[i] = data.minMyGID + i
        end
    else
        for i = 1:data.numMyElements
            myGlobalElements[i] = data.myGlobalElements[i]
        end
    end
    
    myGlobalElements
end



# DECISION depends on element size
#"""
#Return list where for each local point, indicates the local
#element ID that the point belongs to
#"""
#function pointToElementList(map::BlockMap)::Array{Integer}
#        #::Array{LID} where LID <: Integer
#    data = map.data
#    
#    firstPointInElementList = Array{Integer}(data.numMyElements)
#    
#    if length(data.firstPointInElementList) == 0
#        firstPointInElementList[0] = 0


function isOneToOne(map::BlockMap)::Bool
    data = map.data
    if !(data.oneToOneIsDetermined)
        data.oneToOne = determineIsOneToOne()
        data.oneToOneIsDetermined = true
    end
    data.oneToOne
end

function determineIsOneToOne(map::BlockMap)::Bool
    data = map.data
    if numProc(data.comm) < 2
        true
    else
        if isnull(data.directory)
            data.directory = Nullable(createDirectory(data.comm, map))
        end
       gidsAllUniquelyOwned(data.directory)
    end 
end
                                                    
"""
Return the Comm for the map                                                    
"""
function Comm(map::BlockMap)::Comm
    map.data.comm
end