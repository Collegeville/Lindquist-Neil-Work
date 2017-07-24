"""
Contains the data for a BlockMap
"""
type BlockMapData
    comm::Comm
    directory::Nullable{Directory}
    lid::Array{Integer}
    myGlobalElements::Array{Integer}
    firstPointInElementList::Array{Integer}
#    elementSizeList::Array{Integer}
    pointToElementList::Array{Integer}
    
    numGlobalElements::Integer
    numMyElements::Integer
    indexBase::Integer
#    elementSize::Integer
#    minMyElementSize::Integer
#    maxMyElementSize::Integer
#    minElementSize::Integer
#    maxElementSize::Integer
    minAllGID::Integer
    maxAllGID::Integer
    minMyGID::Integer
    maxMyGID::Integer
    minLID::Integer
    maxLID::Integer
    numGlobalPoints::Integer
    numMyPoints::Integer
    
#    constantElementSize::Bool
    linearMap::Bool
    distributedGlobal::Bool
    oneToOneIsDetermined::Bool
    oneToOne::Bool
    lastContiguousGID::Integer
    lastContiguousGIDLoc::Integer
    LIDHash::Dict{Integer, Integer}
end