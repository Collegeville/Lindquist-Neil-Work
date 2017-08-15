export getProfileType, getColMap
export resumeFill, fillComplete
export insertLocalIndices, insertGlobalIndices

#TODO document these methods

getRowMap(graph::CRSGraph) = graph.rowMap
getColMap(graph::CRSGraph) = get(graph.colMap)
getDomainMap(graph::CRSGraph) = get(graph.domainMap)
getRangeMap(graph::CRSGraph) = get(graph.rangeMap)
getImporter(graph::CRSGraph) = get(graph.importer)
getExporter(graph::CRSGraph) = get(graph.exporter)


function insertLocalIndices(graph::CRSGraph{GID, PID, LID}, localRow::Integer,
        numEntries::Integer, inds::Array{<:Integer, 1}) where {GID, PID, LID <: Integer}
    insertLocalIndices(graph, LID(localRow), LID(numEntries), Array{LID, 1}(inds))
end
function insertLocalIndices(graph::CRSGraph{GID, PID, LID}, localRow::LID,
        numEntries::LID, inds::Array{LID, 1}) where {
        GID, PID, LID <: Integer}
    indicesView = view(inds, 1:numEntries)
    insertLocalIndices(graph, localRow, indsT)
end

function insertLocalIndices(graph::CRSGraph{GID, PID, LID}, localRow::Integer,
        inds::Array{<:Integer, 1}) where {GID, PID, LID <: Integer}
    insertLocalIndices(graph, LID(localRow), Array{LID, 1}(inds))
end
function insertLocalIndices(graph::CRSGraph{GID, PID, LID},
        localRow::LID, indices::Array{LID, 1}) where{
        GID, PID, LID <: Integer}
    if !isFillActive(graph)
        throw(InvalidStateError("insertLocalIndices requires that fill is active"))
    end
    if isGloballyIndexed(graph)
        throw(InvalidStateError("graph indices are global, use insertGlobalIndices(...) instead"))
    end
    if !hasColMap(graph)
        throw(InvalidStateError("Cannot insert local indices without a column map"))
    end
    if !myLID(map(graph), localRow)
        throw(InvalidArgumentError("Row does not belong to this process"))
    end
    if !hasRowInfo(graph)
        throw(InvalidStateError("Row information was deleted"))
    end
    
    debug = @debug graph
    if debug
        colMap = getColMap(graph)
        badColIndices = [ind for ind in indices if myLID(colMap, ind)]
        
        if length(badColIndices) != 0
            throw(InvalidArgumentError(
                "Attempting to insert entries in owned row $localRow, "
                * "at the following column indices: $indices.\n"
                                
                * "Of those, the following indices are not in "
                * "the column map on this process: $badColIndices.\n"
        
                * "Since the graph has a column map already, it is "
                * "invalid to insert entries at those locations"))
        end
    end
       
    insertLocalIndicesImpl(graph, localRow, indices)
    
    if debug
        @assert isLocallyIndexed(graph) "Post condtion violated"
    end
end    
    

function insertGlobalIndices(graph::CRSGraph{GID, PID, LID}, globalRow::Integer,
        numEntries::Integer, inds::AbstractArray{<: Integer, 1}) where {
        GID <: Integer, PID, LID <: Integer}
    insertGlobalIndices(graph, GID(globalRow), LID(numEntries), Array{GID, 1}(inds))
end

function insertGlobalIndices(graph::CRSGraph{GID, PID, LID}, globalRow::GID,
        numEntries::LID, inds::AbstractArray{GID, 1}) where {
        GID <: Integer, PID, LID <: Integer}
    indicesView = view(inds, 1:numEntries)
    insertGlobalIndices(graph, globalRow, indsT)
end

function insertGlobalIndices(graph::CRSGraph{GID, PID, LID}, globalRow::Integer,
        inds::AbstractArray{<: Integer, 1}) where {GID <: Integer, PID, LID <: Integer}
    insertGlobalIndices(graph, GID(globalRow), Array{GID, 1}(inds))
end

function insertGlobalIndices(graph::CRSGraph{GID, PID, LID}, globalRow::GID,
        indices::AbstractArray{GID, 1}) where {GID <: Integer, PID, LID <: Integer}
    if isLocallyIndexed(graph)
        throw(InvalidStateError("Graph indices are local, use insertLocalIndices()"))
    end
    if !hasRowInfo(graph)
        throw(InvalidStateError("Graph row information was deleted"))
    end
    if isFillComplete(graph)
        throw(InvalidStateError("Cannot call this method if the fill is not active"))
    end
    
    myRow = lid(graph.rowMap, globalRow)
    if myRow != 0
        if @debug graph
            if hasColMap(graph)
                colMap = get(graph.colMap)

                #appearently jupyter can't render the generator if properly
                badColInds = [index for index in indices 
                                        if myGid(colMap, index)==0]
                if length(badColInds) != 0
                    throw(InvalidArgumentError("$(myPid(comm(graph))): "
                        * "Attempted to insert entries in owned row $globalRow, "
                        * "at the following column indices: $indices.\n"

                        * "Of those, the following indices are not in the "
                        * "column Map on this process: $badColInds.\n"

                        * "Since the matrix has a column map already, it "
                        * "is invalid to insert entries at those locations"))
               end
           end
        end
        insertGlobalIndicesImpl(graph, myRow, indices)
    else
        append!(graph.nonlocalRow, indices)
    end
end


function insertGlobalIndicesFiltered(graph::CRSGraph{GID, PID, LID}, globalRow::GID,
        indices::AbstractArray{GID, 1}) where{GID, PID, LID}
    if isLocallyIndexed(graph)
        throw(InvalidStateError(
                "graph indices are local, use insertLocalIndices(...)"))
    end
    if !hasRowInfo(graph)
        throw(InvalidStateError("Graph row information was deleted"))
    end
    if isFillComplete(graph)
        throw(InvalidStateError("Cannot insert into fill complete graph"))
    end
    
    myRow = lid(graph.rowMap, globalRow)
    if myRow != 0
        #if column map present, use it to filter the entries
        if hasColMap(graph)
            colMap = getColMap(graph)
            indices = [index for index in indices if myLid(colMap, index)]
        end
        insertGlobalIndicesImpl(myRow, indices)
    else
        #nonlocal row
        append!(graph.nonlocals[globalRow], indices)
    end
end
                        
                        

function getNumEntriesInGlobalRow(graph::CRSGraph{GID}, globalRow::Integer)::Integer where {GID <: Integer}
    localRow = lid(graph.rowMap, GID(globalRow))
    if hasRowInfo(graph) && localRow != 0
        getRowInfo(graph, localRow).numEntries
    else
        -1
    end
end

function getNumEntriesInLocalRow(graph::CRSGraph{GID, PID, LID}, localRow::Integer)::Integer where {GID, PID, LID <: Integer}
    if hasRowInfo(graph) && myLID(graph.rowMap, LID(localRow))
        getRowInfo(graph, LID(localRow)).numEntries
    else
        -1
    end
end

function getGlobalView(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}) where {GID <: Integer, PID, LID <: Integer}
    if rowInfo.allocSize > 0
        if length(graph.globalIndices1D) != 0
            range = rowInfo.offset1D : rowInfo.offset1D + rowInfo.allocSize
            view(graph.globalIndices1D, range)
        elseif length(graph.globalIndices2D[rowInfo.localRow]) == 0
            globalIndices2D[rowInfo.localRow]
        else
            GID[]
        end
    else
        GID[]
    end
end
        
function getLocalView(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}) where {GID <: Integer, PID, LID <: Integer}
    if rowInfo.allocSize > 0
        if length(graph.localIndices1D) != 0
            range = rowInfo.offset1D : rowInfo.offset1D + rowInfo.allocSize
            return view(graph.localIndices1D, range)
        elseif length(graph.localIndices2D[rowInfo.localRow]) == 0
            return localIndices2D[rowInfo.localRow]
        end
    end
    return LID[]
end

function getGlobalRowCopy(graph::CRSGraph{GID}, globalRow::GID)::Array{GID, 1} where {GID <: Integer}
    Array{GID, 1}(getGlobalRowView(graph, globalRow))
end
        
function getLocalRowCopy(graph::CRSGraph{GID, PID, LID}, globalRow::LID)::Array{LID, 1} where {GID, PID, LID <: Integer}
    Array{LID, 1}(getLocalRowView(graph, globalRow))
end

function getGlobalRowView(graph::CRSGraph{GID}, globalRow::GID)::AbstractArray{GID, 1} where {GID <: Integer}
    debug = @debug graph
    if isLocallyIndexed(graph)
        throw(InvalidArgumentError("The graph's indices are currently stored as local indices, so a view with global column indices cannot be returned.  Use getGlobalRowCopy(::CRSGraph) instead"))
    end

    if debug
        @assert hasRowInfo(graph) "Graph row information was deleted"
    end
    rowInfo = getRowInfoFromGlobalRow(graph, globalRow)
    
    if rowInfo.localRow != 0 && rowInfo.numEntries > 0
        indices = view(getGlobalView(graph, rowInfo), 1:rowInfo.numEntries)
        if debug
            @assert(length(indices) == getNumEntriesInGlobalRow(graph, globalRow),
                "length(indices) = $(length(indices)) "
                * "!= getNumEntriesInGlobalRow(graph, $globalRow) "
                * "= $(getNumEntriesInGlobalRow(graph, globalRow))")
        end
        indices
    else
        GID[]
    end
end
        
function getLocalRowView(graph::CRSGraph{GID}, localRow::GID)::AbstractArray{GID, 1} where {GID <: Integer}
    debug = @debug graph
    if isGloballyIndexed(graph)
        throw(InvalidArgumentError("The graph's indices are currently stored as global indices, so a view with local column indices cannot be returned.  Use getLocalRowCopy(::CRSGraph) instead"))
    end

    if debug
        @assert hasRowInfo() "Graph row information was deleted"
    end
    rowInfo = getRowInfoFromLocalRowIndex(localRow)
    
    if rowInfo.localRow != 0 && rowInfo.numEntries > 0
        indices = view(getLocalView, 1:rowInfo.numEntries)
        if debug
            @assert(length(indices) == getNumEntriesInLocalRow(graph, localRow),
                "length(indices) = $(length(indices)) "
                * "!= getNumEntriesInLocalRow(graph, $localRow) "
                * "= $(getNumEntriesInLocalRow(graph, localRow))")
        end
        indices
    else
        LID[]
    end
end
    


resumeFill(graph::CRSGraph; plist...) = resumeFill(graph, Dict(Array{Tuple{Symbol, Any}, 1}(plist)))

function resumeFill(graph::CRSGraph, plist::Dict{Symbol})
    if !hasRowInfo(graph)
        throw(InvalidStateError("Cannot resume fill of the CRSGraph, "
                * "since the graph's row information was deleted."))
    end
    
    clearGlobalConstants(graph)
    graph.plist = plist
    graph.lowerTriangle = false
    graph.upperTriangle = false
    graph.indicesAreSorted = true
    graph.noRedundancies = true
    graph.fillComplete = false
    
    if get(plist, :debug, false)
        @assert isFillActive(graph) && !isFillComplete(graph) "Post condition violated"
    end
end


fillComplete(graph::CRSGraph; plist...) = fillComplete(graph, Dict(Array{Tuple{Symbol, Any}, 1}(plist)))

function fillComplete(graph::CRSGraph, plist::Dict{Symbol})
    if isnull(graph.domainMap)
        domMap = graph.rowMap
    else
        domMap = get(graph.domainMap)
    end
    
    if isnull(graph.rangeMap)
        ranMap = graph.colMap
    else
        ranMap = get(graph.rangeMap)
    end
    
    fillComplete(graph, ranMap, domMap, plist)
end

function fillComplete(graph::CRSGraph{GID, PID, LID},
        domainMap::BlockMap{GID, PID, LID}, rangeMap::BlockMap{GID, PID, LID};
        plist...) where {GID, PID, LID}
    fillComplete(graph, Dict(Array{Tuple{Symbol, Any}, 1}(plist)))
end

function fillComplete(graph::CRSGraph{GID, PID, LID}, domainMap::BlockMap{GID, PID, LID}, rangeMap::BlockMap{GID, PID, LID}, plist::Dict{Symbol}) where {GID, PID, LID}
    if !isFillActive(graph) || isFillComplete(graph)
        throw(InvalidStateError("Graph fill state must be active to call fillComplete(...)"))
    end
       
    const numProcs = numProc(comm(graph))
    
    if haskey(plist, :sortColumnMapGhostGIDs)
        graph.sortGhostsAssociatedWithEachProcessor = get(plist, :sortColumnMapGhostGIDs, :VERY_BAD)
    end
    
    const assertNoNonlocalInserts = get(plist, :noNonlocalChanges, false)
    
    const mayNeedGlobalAssemble = !assertNoNonlocalInserts && numProcs > 1
    if mayNeedGlobalAssemble
        globalAssemble(graph)
    else
        if numProcs == 1 && length(graph.nonLocals) > 0
            throw(InvalidStateError("Only one process, but nonlocal entries are present"))
        end
    end
    
    setDomainRangeMaps(graph, domainMap, rangeMap)

    if !hasColMap(graph)
        makeColMap(graph)
    end

    makeIndicesLocal(graph)

    sortAndMergeAllIndices(graph, isSorted(graph), isMerged(graph))
    
    makeImportExport(graph)
    computeGlobalConstants(graph)
    fillLocalGraph(graph, plist)
    graph.fillComplete(true)
    
    if get(plist, :debug, false)
        @assert !isFillActive(graph) && isFillComplete(graph) "post conditions violated"
    end
    
    checkInternalState(graph)
end

function makeColMap(graph::CRSGraph{GID, PID, LID}) where {GID, PID, LID}
    debug = @debug graph
    const localNumRows = getNodeNumElements(graph)
    
    #TODO get rid of this order retention stuff, it has to do with epetra interop
    const sortEachProcsGIDs = graph.sortGhostsAssociatedWithEachProcessr
    
    #TODO look at FIXME on line 4898
    
    errCode, colMap = __makeColMap(graph, graph.domainMap, sortEachProcsGIDs)
    if debug
        comm = julia_petra.comm(graph)
        localSuccess = (errCode == 0)? 1 : 0
        globalSuccess = minAll(comm, localSuccess)
        
        if globalSuccess != 1
            error("makeColMap reports an error on at least one process")
        end
    end
    
    graph.colMap = colMap
    
    checkInternalState(graph)
end

"""
hasColMap(::CRSGraph)

Whether the graph has a column map
"""
hasColMap(graph::CRSGraph) = !isnull(graph.colMap)

"""
    isSorted(::CRSGraph)

Whether the indices are sorted
"""
isSorted(graph::CRSGraph) = graph.indicesAreSorted

"""
    isMerged(::CRSGraph)

Whether duplicate column indices in each row have been merged
"""
isMerged(graph::CRSGraph) = graph.noRedundancies

"""
    setAllIndices(graph::CRSGraph{GID, PID, LID}, rowPointers::Array{LID, 1}, columnIndices::Array{LID, 1})

Sets the graph's data directly, using 1D storage
"""
function setAllIndices(graph::CRSGraph{GID, PID, LID},
        rowPointers::Array{LID, 1},columnIndices::Array{LID, 1}) where {
        GID, PID, LID <: Integer}
    
    localNumRows = getNodeNumRows(graph)
    
    if isnull(graph.colMap)
        throw(InvalidStateError("The graph must have a "
                * "column map before calling setAllIndices"))
    end
    if length(rowPointers) != localNumRows + 1
        throw(InvalidArgumentError("length(rowPointers) = $(length(rowPointers)) "
                * "!= localNumRows+1 = $(localNumRows+1)"))
    end
    
    localNumEntries = rowPointers[localNumRows+1]
    
    graph.indicesType    = LOCAL_INDICES
    graph.pftype         = STATIC_PROFILE
    graph.localIndices1D = columnIndices
    graph.rowOffsets     = rowPointers
    graph.nodeNumEntries = localNumEntries
    graph.storageStatus  = STORAGE_1D_UNPACKED
    
    graph.localGraph     = LocalCRSGraph(columnIndices, rowPointers)
    
    checkInternalState(graph)
end


"""
    isStorageOptimized(::CRSGraph)

Whether the graph's storage is optimized
"""
function isStorageOptimized(graph::CRSGraph)
    const isOpt = length(graph.numRowEntries) == 0 && getNodeNumRows(graph) > 0
    if isOpt && @debug graph
        @assert(getProfileType(graph) == STATIC_PROFILE,
            "Matrix claims optimized storage by profile type "
            * "is dynamic.  This shouldn't happend.")
    end
    isOpt
end

"""
    isFillComplete(graph)::Bool

Whether the graph is fill complete
"""
isFillComplete(g::CRSGraph) = g.fillComplete

"""
    getNodeNumRows(graph)

Gets the number of rows on this processor
"""
function getNodeNumRows(graph::CRSGraph{<:Integer, <:Integer, LID}) where LID <: Integer
    numMyElements(graph.rowMap)
end
    

"""
    isGloballyIndexed(::CRSGraph)

Whether the graph uses global indexes
"""
function isGloballyIndexed(graph::CRSGraph)
    graph.indicesType == GLOBAL_INDICES
end
   
"""
    isLocallyIndexed(::CRSGraph)

Whether the graph uses local indexes
"""
function isLocallyIndexed(graph::CRSGraph)
    graph.indicesType == LOCAL_INDICES
end

"""
    getProfileType(::CRSGraph)

Gets the profile type of the graph
"""
function getProfileType(graph::CRSGraph)
    graph.pftype
end


function map(graph::CRSGraph)
    graph.rowMap
end


#TODO implement RowGraph methods

#### DistObject methods ####
function checkSizes(source::RowGraph{GID, PID, LID},
        target::CRSGraph{GID, PID, LID}) where {GID, PID, LID}
    #T and E petra's don't do any checks
    true
end

function copyAndPermute(source::RowGraph{GID, PID, LID},
        target::CRSGraph{GID, PID, LID}, numSameIDs::LID,
        permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1}) where {
        GID, PID, LID}
    copyAndPermuteNoViewMode(source, target,
        numSameIDs, permuteToLIDs, permuteFromLIDs)
end

function copyAndPermute(source::CRSGraph{GID, PID, LID},
        target::CRSGraph{GID, PID, LID}, numSameIDs::LID,
        permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1}) where {
        GID, PID, LID}
    if isFillComplete(target)
        throw(InvalidStateError("Target cannot be fill complete"))
    end
    if isFillComplete(source)
            copyAndPermuteNoViewMode(source, target,
        numSameIDs, permuteToLIDs, permuteFromLIDs)
    else
        if length(permuteToLIDs) != length(permuteFromLIDs)
            throw(InvalidArgumentError("permuteToLIDs and "
                    * "permuteFromLIDs must have the same size"))
        end
        
        srcRowMap = getRowMap(source)
        tgtRowMap = getRowMap(target)
        
        #copy part
        for myid = 1:numSameIDs
            myGID = gid(srcRowMap, myid)
            row = getGlobalRowView(source, myGID)
            insertGlobalIndices(target, myGID, row)
        end
        
        #permute part
        for i = 1:length(permuteToLIDs)
            srcGID = gid(srcRowMap, permuteFromLIDs[i])
            tgtGID = gid(tgtRowMap, permuteToLIDs[i])
            row = getGlobalRowView(source, srcGID)
            insertGlobalIndices(target, tgtGID, row)
        end
    end
end
        
            
function copyAndPermuteNoViewMode(source::RowGraph{GID, PID, LID},
        target::CRSGraph{GID, PID, LID}, numSameIDs::LID,
        permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1}) where {
        GID, PID, LID}
    if length(permuteToLIDs) != length(premuteFromLIDs)
        throw(InvalidArgumentError("permuteToLIDs and "
                * "permuteFromLIDs must have the same size"))
    end
    
    srcRowMap = getRowMap(source)
    tgtRowMap = getRowMap(target)
    
    #copy part
    for myid = 1:numSameIDs
        myGID = gid(srcRowMap, myid)
        rowCopy = getGlobalRowCopy(sourceRowGraph, myGID)
        insertGlobalIndices(target, myGID, rowCopy)
    end
    
    #permute part
    for i = 1:length(permuteToLIDs)
        tgtGID = gid(tgtRowMap, permuteToLIDs[i])
        srcGID = gid(srcRowMap, permuteFromLIDs[i])
        rowCopy = getGlobalRowCopy(source, srcGID)
        insertGlobalIndices(target, tgtGID, rowCopy)
    end
end
        
    

function packAndPrepare(source::RowGraph{GID, PID, LID},
        target::CRSGraph{GID, PID, LID}, exportLIDs::Array{LID, 1},
        distor::Distributor{GID, PID, LID})::Array{GID, 1} where {
        GID, PID, LID}

    pack(source, exportLIDs, distor)
end

function unpackAndCombine(target::CRSGraph{GID, PID, LID},
        importLIDs::Array{LID, 1}, imports::Array,
        distor::Distributor{GID, PID, LID}, cm::CombineMode) where {
        GID, PID, LID}
    #should be caught else where
    @assert(isFillActive(target),
        "Import and Export operations require a fill active graph")
    
    tgtMap = map(target)
    
    for i = 1:length(importLIDs)
        row = imports[i]
        insertGlobalIndicesFiltered(target, gid(tgtMap, importLIDs[i]), row)
    end
end


function pack(source::CRSGraph{GID, PID, LID}, exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID})::Array{Array{LID, 1}, 1} where {GID, PID, LID}
    srcMap = map(source)
    [getGlobalRowCopy(source, gid(srcMap, lid)) for lid in exportLIDs]
end
    
    
