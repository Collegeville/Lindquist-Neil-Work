
export BlockMap
export remoteIDList, lid, gid, findLocalElementID
export minAllGID, maxAllGID, minMyGID, maxMyGID, minLID, maxLID
export numGlobalElements, myGlobalElements
export uniqueGIDs, globalIndicesType, sameBlockMapDataAs, sameAs
export linearMap, myGlobalElementIDs, comm
export myGID, myLID, distributedGlobal, numMyElements
#export pointToElementList, numGlobalPoints, numMyPoints, pointSameAs


# methods and docs based straight off Epetra_BlockMap to match Comm

# ignoring indexBase methods and sticking with 1-based indexing
# ignoring elementSize methods since, type information is carried anyways
# ignoring point-related code, since elementSize is ignored


# TODO figure out expert users and developers only functions


"""
A type for partitioning block element vectors and matrices
"""
type BlockMap{GID <: Integer, PID <:Integer, LID <: Integer}
    data::BlockMapData{GID, PID, LID}
    
    function BlockMap{GID, PID, LID}(data::BlockMapData) where GID <: Integer where PID <:Integer where LID <: Integer
        new(data)
    end
end


"""
constructor for Epetra-defined uniform linear distribution of elements
"""
function BlockMap(numGlobalElements::Integer, comm::Comm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    BlockMap(GID(numGlobalElements), comm)
end

function BlockMap(numGlobalElements::GID, comm::Comm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    if numGlobalElements < 0 
        throw(InvalidArgumentError("NumGlobalElements = $(numGlobalElements).  Should be >= 0"))
    end
    
    const data = BlockMapData(numGlobalElements, comm)
    const map = BlockMap{GID, PID, LID}(data)
    
    numProcVal = numProc(comm)
    data.linearMap = true
    
    myPIDVal = myPid(comm) - 1
    
    data.numMyElements = floor(typeof(data.numGlobalElements),
        data.numGlobalElements/numProcVal)
    remainder = data.numGlobalElements % numProcVal
    startIndex = myPIDVal * (data.numMyElements+1)
    
    if myPIDVal < remainder
        data.numMyElements += 1
    else
        startIndex -= (myPIDVal - remainder)
    end
    
    data.minAllGID = 1
    data.maxAllGID = data.minAllGID + data.numGlobalElements - 1
    data.minMyGID = startIndex + 1
    data.maxMyGID = data.minMyGID + data.numMyElements - 1
    data.distributedGlobal = isDistributedGlobal(map, data.numGlobalElements,
        data.numMyElements)
    
    EndOfConstructorOps(map)
    map
end

"""
constructor for user-defined linear distribution of elements
"""
function BlockMap(numGlobalElements::Integer, numMyElements::Integer, comm::Comm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    BlockMap(GID(numGlobalElements), LID(numMyElements), comm)
end

function BlockMap(numGlobalElements::GID, numMyElements::LID, comm::Comm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    if numGlobalElements < -1 
        throw(InvalidArgumentError("NumGlobalElements = $(numGlobalElements).  Should be >= -1"))
    end
    if numMyElements < 0
        throw(InvalidArgumentError("NumMyElements = $(numMyElements). Should be >= 0"))
    end
    
    const data = BlockMapData(numGlobalElements, comm)
    const map = BlockMap{GID, PID, LID}(data)
    
    data.numMyElements = numMyElements
    data.linearMap = true
    
    data.distributedGlobal = isDistributedGlobal(map, numGlobalElements, numMyElements)
    
    #Local Map and uniprocessor case: Each processor gets a complete copy of all elements
    if !data.distributedGlobal || numProc(comm) == 1
        data.numGlobalElements = data.numMyElements
        
        data.minAllGID = 1
        data.maxAllGID = data.minAllGID + data.numGlobalElements - 1
        data.minMyGID = 1
        data.maxMyGID = data.minMyGID + data.numMyElements - 1
    else
        tmp_numMyElements = data.numMyElements
        data.numGlobalElements = sumAll(data.comm, tmp_numMyElements)
        
        data.minAllGID = 1
        data.maxAllGID = data.minAllGID + data.numGlobalElements - 1
        
        tmp_numMyElements = data.numMyElements
        data.maxMyGID = scanSum(data.comm, tmp_numMyElements)
        
        startIndex = data.maxMyGID - data.numMyElements
        data.minMyGID = startIndex + 1
        data.maxMyGID = data.minMyGID + data.numMyElements - 1
    end
    checkValidNGE(map, numGlobalElements)
    
    EndOfConstructorOps(map)
    map
end


"""
constructor for user-defined arbitrary distribution of elements
"""
function BlockMap(numGlobalElements::Integer, numMyElements::Integer,
        myGlobalElements::Array{<:Integer}, comm::Comm{GID, PID,LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    BlockMap(GID(numGlobalElements), LID(numMyElements), Array{GID}(myGlobalElements), comm)
end

function BlockMap(numGlobalElements::GID, numMyElements::LID,
        myGlobalElements::Array{GID}, comm::Comm{GID, PID,LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    if numGlobalElements < -1 
        throw(InvalidArgumentError("NumGlobalElements = $(numGlobalElements).  Should be >= -1"))
    end
    if numMyElements < 0
        throw(InvalidArgumentError("NumMyElements = $(numMyElements). Should be >= 0"))
    end
    
    const data = BlockMapData(numGlobalElements, comm)
    const map = BlockMap{GID, PID, LID}(data)
    
    data.numMyElements = numMyElements
    
    linear = 1
    if numMyElements > 0
        data.myGlobalElements = Array{GID}(numMyElements)
        
        data.myGlobalElements[1] = myGlobalElements[1]
        data.minMyGID = myGlobalElements[1]
        data.maxMyGID = myGlobalElements[1]
        
        for i = 2:numMyElements
            data.myGlobalElements[i] = myGlobalElements[i]
            data.minMyGID = min(data.minMyGID, myGlobalElements[i])
            data.maxMyGID = max(data.maxMyGID, myGlobalElements[i])
            
            if myGlobalElements[i] != myGlobalElements[i-1] + 1
                linear = 0
            end
        end
    else
        data.minMyGID = 1
        data.maxMyGID = 0
    end
    
    data.linearMap = Bool(minAll(data.comm, linear))
    
    data.distributedGlobal = isDistributedGlobal(map, numGlobalElements, numMyElements)
    
    if !data.distributedGlobal || numProc(comm) == 1
        data.numGlobalElements = data.numMyElements
        checkValidNGE(map, numGlobalElements)
        data.minAllGID = data.minMyGID
        data.maxAllGID = data.maxMyGID
    else
        data.numGlobalElements = sumAll(data.comm, data.numMyElements)
        checkValidNGE(map, numGlobalElements)
        
        tmp_send = [
            -((data.numMyElements > 0
                    || data.numGlobalElements == 0)?
                data.minMyGID:Inf)
            , data.maxMyGID]
        
        tmp_recv = maxAll(data.comm, tmp_send)
        
        @assert typeof(tmp_recv[1]) <: Integer "Result type is $(typeof(tmp_recv[1])), should be subtype of Integer"
            
        data.minAllGID = -tmp_recv[1]
        data.maxAllGID =  tmp_recv[2]
            
    end
    
    EndOfConstructorOps(map)
    map
end

"""
constructor for user-defined arbitrary distribution of elements
will all information on globals provided by the user
"""
function BlockMap(numGlobalElements::Integer, numMyElements::Integer,
        myGlobalElements::Array{GID}, userIsDistributedGlobal::Bool,
        userMinAllGID::Integer, userMaxAllGID::Integer, comm::Comm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    BlockMap(GID(numGlobalElements), LID(numMyElements), Array{GID}(myGlobalElements), userIsDistributedGlobal,
        GID(userMinAllGID), GID(userMaxAllGID), comm)
end

function BlockMap(numGlobalElements::GID, numMyElements::LID,
        myGlobalElements::Array{GID}, userIsDistributedGlobal::Bool,
        userMinAllGID::GID, userMaxAllGID::GID, comm::Comm{GID, PID, LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    if numGlobalElements < -1 
        throw(InvalidArgumentError("NumGlobalElements = $(numGlobalElements).  Should be >= -1"))
    end
    if numMyElements < 0
        throw(InvalidArgumentError("NumMyElements = $(numMyElements). Should be >= 0"))
    end
    if userMinAllGID < 1 
        throw(InvalidArgumentError("Minimum global element index = $(data.minAllGID).  Should be >= 1"))
    end
    
    const data = BlockMapData(numGlobalElements, comm)
    const map = BlockMap{GID, PID, LID}(data)
    
    data.numMyElements = numMyElements
    
    linear = 1
    if numMyElements > 0
        data.myGlobalElements = Array{GID}(numMyElements)
        
        data.myGlobalElements[1] = myGlobalElements[1]
        data.minMyGID = myGlobalElements[1]
        data.maxMyGID = myGlobalElements[1]
        
        for i = 2:numMyElements
            data.myGlobalElements[i] = myGlobalElements[i]
            data.minMyGID = min(data.minMyGID, myGlobalElements[i])
            data.maxMyGID = max(data.maxMyGID, myGlobalElements[i])
            
            if myGlobalElements[i] != myGlobalElements[i-1] + 1
                linear = 0
            end
        end
        
    else
        data.minMyGID = 1
        data.maxMyGID = 0
    end
    
    data.linearMap = Bool(minAll(comm, linear))
    
    data.distributedGlobal = userIsDistributedGlobal
    
     if !data.distributedGlobal || numProc(comm) == 1
        data.numGlobalElements = data.numMyElements
        checkValidNGE(map, numGlobalElements)
        data.minAllGID = data.minMyGID
        data.maxAllGID = data.maxMyGID
    else
        if numGlobalElements == -1
            data.numGlobalElements = sumAll(data.comm, data.numMyElements)
        else
            data.numGlobalElements = numGlobalElements
        end
        checkValidNGE(data.numGlobalELements)
        
        data.minAllGID = userMinAllGID
        data.maxAllGID = userMaxAllGID
    end
    EndOfConstructorOps(map)
    map
end

"""
Copy constructor
"""
function BlockMap(map::BlockMap)
    BlockMap(map.data)
end



##### internal construction methods #####
function isDistributedGlobal(map::BlockMap{GID, PID, LID}, numGlobalElements::GID,
        numMyElements::LID) where GID <: Integer where PID <: Integer where LID <: Integer
    data = map.data
    if numProc(data.comm) > 1 
        localReplicated = numGlobalElements == numMyElements
        !Bool(minAll(data.comm, localReplicated))
    else
        false
    end
end

function EndOfConstructorOps(map::BlockMap)
    map.data.minLID = 1
    map.data.maxLID = max(map.data.numMyElements, 1)
    
    GlobalToLocalSetup(map);
end

function GlobalToLocalSetup(map::BlockMap)
    data = map.data
    numMyElements = data.numMyElements
    myGlobalElements = data.myGlobalElements
    
    if data.linearMap || numMyElements == 0
        return map
    end
    if length(data.numGlobalElements) == 0
        return map
    end
    
    
    val = myGlobalElements[1]
    i = 1
    for i = 1:numMyElements
        if val != myGlobalElements[i]
            break
        end
        val += 1
    end
    
    data.lastContiguousGIDLoc = i
    if data.lastContiguousGIDLoc <= 1
        data.lastContiguousGID = myGlobalElements[1]
    else
        data.lastContiguousGID = myGlobalElements[data.lastContiguousGIDLoc]
    end
    
    if i < numMyElements
        data.lidHash = empty!(data.lidHash)
        
        sizehint!(data.lidHash, numMyElements - i + 2)
        
        for i = i:numMyElements
            data.lidHash[myGlobalElement[i]] = i
        end
    end
    map
end

function checkValidNGE(map::BlockMap{GID, PID, LID}, numGlobalElements::GID) where GID <: Integer where PID <: Integer where LID <: Integer
    if (numGlobalElements != -1) && (numGlobalElements != map.data.numGlobalElements)
        throw(InvalidArgumentError("Invalid NumGlobalElements.  "
              * "NumGlobalElements = $(numGlobalElements)"
              * ".  Should equal $(map.data.numGlobalElements)"
              * ", or be set to -1 to compute automatically"))
    end
end

##### external methods #####

"""
Return true if the GID passed in belongs to the calling processor in this
map, otherwise returns false.
"""
function myGID(map::BlockMap{GID, PID, LID}, gidVal::Integer) where GID <: Integer where PID <: Integer where LID <: Integer
    lid(map, gidVal) != 0
end

"""
Return true if the LID passed in belongs to the calling processor in this
map, otherwise returns false.
"""
function myLID(map::BlockMap{GID, PID, LID}, lidVal::Integer) where GID <: Integer where PID <: Integer where LID <: Integer
    gid(map, lidVal) != 0
end

"""
Return true if map is defined across more than one processor
"""
function distributedGlobal(map::BlockMap)
    map.data.distributedGlobal
end

"""
Return the number of elements across the calling processor
"""
function numMyElements(map::BlockMap{GID, PID, LID})::LID where GID <: Integer where PID <: Integer where LID <: Integer
    map.data.numMyElements
end

"""
Return the minimum global ID owned by this processor
"""
function minMyGID(map::BlockMap{GID, PID, LID})::GID where GID <: Integer where PID <: Integer where LID <: Integer
    map.data.minMyGID
end
    
"""
Return the maximum global ID owned by this processor
"""
function maxMyGID(map::BlockMap{GID, PID, LID})::GID where GID <: Integer where PID <: Integer where LID <: Integer
    map.data.maxMyGID
end


##local/global ID accessor methods##

"""
Return the processor ID and local index value for a given list of global indices.
The returned value is a tuple containing
    1 - an Array of processors owning the global ID's in question
    2 - an Array of local IDs of the global on the owning processor
"""
function remoteIDList(map::BlockMap{GID, PID, LID}, gidList::Array{<:Integer}
        )::Tuple{Array{PID}, Array{LID}} where GID <: Integer where PID <: Integer where LID <: Integer
    remoteIDList(map, Array{GID}(gidList))
end

function remoteIDList(map::BlockMap{GID, PID, LID}, gidList::Array{GID}
        )::Tuple{Array{PID}, Array{LID}} where GID <: Integer where PID <: Integer where LID <: Integer
    data = map.data
    if isnull(data.directory)
        data.directory = createDirectory(data.comm, map)
    end

    getDirectoryEntries(get(data.directory), map, gidList)
end

"""
Return local ID of global ID, or 0 if not found on this processor
"""
function lid(map::BlockMap{GID, PID, LID}, gid::Integer)::LID where GID <: Integer where PID <: Integer where LID <: Integer
    data = map.data
    if (gid < data.minMyGID) || (gid > data.maxMyGID)
        return 0
    end
    if data.linearMap
        return gid - data.minMyGID + 1
    end
    if gid >= data.myGlobalElements[1] && gid <= data.lastContiguousGID
        return gid - data.myGlobalElements[1] + 1
    end

    return data.lidHash[gid]
end

"""
Return global ID of local ID, or 0 if not found on this processor
"""
function gid(map::BlockMap{GID, PID, LID}, lid::Integer)::GID where GID <: Integer where PID <: Integer where LID <: Integer 
    data = map.data
    if (data.numMyElements == 0) || (lid < data.minLID) || (lid > data.maxLID)
        return 0
    end
    
    if data.linearMap
        return lid + data.minMyGID - 1
    end
    return data.myGlobalElements[lid]
end


"""
Return the minimum global ID across the entire map
"""
function minAllGID(map::BlockMap{GID})::GID where GID <: Integer
    map.data.minAllGID
end

"""
Return the maximum global ID across the entire map
"""
function maxAllGID(map::BlockMap{GID})::GID where GID <: Integer
    map.data.maxAllGID
end

"""
Return the mimimum local index value on the calling processor
"""
function minLID(map::BlockMap{GID, PID, LID})::LID where GID <: Integer where PID <: Integer where LID <: Integer
    map.data.minLID
end

"""
Return the maximum local index value on the calling processor
"""
function maxLID(map::BlockMap{GID, PID, LID})::LID where GID <: Integer where PID <: Integer where LID <: Integer
    map.data.maxLID
end

##size/dimension accessor functions##

"""
Return the number of elements across all processors
"""
function numGlobalElements(map::BlockMap{GID})::GID where GID <: Integer
    map.data.numGlobalElements
end

"""
Return a list of global elements on this processor
"""
function myGlobalElements(map::BlockMap{GID})::Array{GID} where GID <: Integer
    data = map.data
    
    if length(data.myGlobalElements) == 0
        myGlobalElements = Array{Integer}(data.numMyElements)
        for i = 1:data.numMyElements
            myGlobalElements[i] = data.minMyGID + i - 1
        end
        data.myGlobalElements = myGlobalElements
    else
        data.myGlobalElements
    end
end

#"""
#Return the number of global points for this map
#"""
#function numGlobalPoints(map::BlockMap)::Integer
#        #::GID where GID <: Integer
#    map.data.numGlobalPoints
#end
#
#"""
#Return the number of local points for this map
#"""
#function numMyPoints(map::BlockMap)::Integer
#        #::LID where LID <: Integer
#    map.data.numMyPoints
#end


##Miscellaneous boolean tests##

"""
Return true if each map GID exists on at most 1 processor
"""
function uniqueGIDs(map::BlockMap)::Bool
    isOneToOne(map)
end


"""
	globalIndicesType(map::BlockMap{GID})::Type{GID}
Return the type used for global indices in the map
"""
function globalIndicesType(map::BlockMap{GID})::Type{GID} where GID <: Integer
	GID
end

"""
Return true if the maps have the same data
"""
function sameBlockMapDataAs(this::BlockMap, other::BlockMap)::Bool
    this.data == other.data
end

"""
Return true if this and other are identical maps
"""
function sameAs(this::BlockMap{GID, PID, LID}, other::BlockMap{GID, PID, LID})::Bool where GID <: Integer where PID <: Integer where LID <: Integer
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
            if gid(this, i) != gid(other, i)
                mySameMap = 0
                break
            end
        end
    end
    
    Bool(minAll(tData.comm, mySameMap))
end
    
#"""
#Return true if this and other have identical point-wise structure
#"""
#function pointSameAs(this::BlockMap, other::BlockMap)::Bool
#    tData = this.data
#    oData = this.data
#    if tData == oData
#        return true
#    end
#    
#    if tData.numGlobalPoints != oData.numGlobalPoints
#        return false
#    end
#    
#    mySameMap = 1
#    if tData.numMyPoints != oData.numMyPoints
#        mySameMap = 0
#    end
#    
#    Bool(minAll(tData.comm, [mySameMap])[1])
#end


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
function myGlobalElementIDs(map::BlockMap{GID})::Array{GID} where GID <: Integer
    data = map.data
    if length(data.myGlobalElements) == 0
        base = 0:data.numMyElements-1
        rng = data.minMyGID + base
        myGlobalElements = collect(rng)
    else
        myGlobalElements = copy(data.myGlobalElements)
    end
    
    myGlobalElements
end



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
        data.oneToOne = determineIsOneToOne(map)
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
       gidsAllUniquelyOwned(get(data.directory))
    end 
end
                                                    
"""
Return the Comm for the map                                                    
"""
function comm(map::BlockMap{GID, PID, LID})::Comm{GID, PID, LID} where GID <: Integer where PID <: Integer where LID <: Integer
    map.data.comm
end
