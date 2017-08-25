#DECISION put this somewhere else?  Its only an internal grouping
struct RowInfo{LID <: Integer}
    graph::CRSGraph{<:Integer, <:Integer, LID}
    localRow::LID
    allocSize::LID
    numEntries::LID
    offset1D::LID
end

#TODO implement getLocalDiagOffsets(::CRSGraph)
getLocalGraph(graph::CRSGraph) = graph.localGraph


function updateGlobalAllocAndValues(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}, newAllocSize::Integer, rowValues::Array{Data, 1})::RowInfo{LID} where {Data, GID, PID, LID}
    
    resize!(graph.globalIndices2D[rowInfo.localRow], newAllocSize)
    resize!(rowVals, newAllocSize)
    
    getRowInfo(graph, rowInfo.localRow)
end

function insertIndicesAndValues(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}, newInds::Union{AbstractArray{GID, 1}, AbstractArray{LID, 1}}, oldRowVals::AbstractArray{Data, 1}, newRowVals::AbstractArray{Data, 1}, lg::IndexType) where {Data, GID, PID, LID}
    numNewInds = insertIndices(graph, rowInfo, newInds, lg)
    oldInd = rowInfo.numEntries+1
    
    oldRowVals[range(oldInd, 1, numNewInds)] = newRowVals[1:numNewInds]
end
    
function insertIndices(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}, newInds::Union{AbstractArray{GID, 1}, AbstractArray{LID, 1}}, lg::IndexType) where {GID, PID, LID}
    numNewInds = 0
    if lg == GLOBAL_INDICES
        newGInds = AbstractArray{GID, 1}(newInds)
        numNewInds = length(newGInds)
        if isGloballyIndexed(graph)
            gIndView = getGlobalView(graph, rowInfo)
            gIndView[range(rowInfo.numEntries+1, 1, length(newGInds))] = newGInds[:]
        else
            lIndView = getLocalView(graph, rowInfo)
            colMap = graph.colMap
            
            dest = range(rowInfo.numEntries, 1, length(newGInds))
            lIndView[dest] = [lid(colMap, gid) for gid in newGInds]
        end
    elseif lg == LOCAL_INDICES
        newLInds = AbstractArray{LID, 1}(newInds)
        numNewInds = length(newLInds)
        if isLocallyIndexed(graph)
            lIndView = getLocalView(graph, rowInfo)
            lIndView[range(rowInfo.numEntries, 1, length(newLInds))] = newLInds[:]
        else
            @assert(false,"lg=LOCAL_INDICES, isGloballyIndexed(g) not implemented, "
            * "because it doesn't make sense")
        end
    end
    
    graph.numRowEntries[rowInfo.localRow] += numNewInds
    graph.nodeNumEntries += numNewInds
    setLocallyModified(graph)
    
    numNewInds
end
            

function computeGlobalConstants(graph::CRSGraph{GID, PID, LID}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    
    #short circuit if already computed
    graph.haveGlobalConstants && return
    
    if @debug graph
        @assert !null(graph.colMap) "The graph must have a column map at this point"
    end
    
    computeLocalConstants(graph)
    
    commObj = comm(map(graph))

    #if graph.haveGlobalConstants == false  #short circuited above
    graph.globalNumEntries, graph.globalNumDiags = sumAll(commObj,
        [GID(graph.nodeNumEntries), GID(graph.nodeNumDiags)])
        
    graph.globalMaxNumRowEntries = maxAll(commObj, GID(graph.nodeMaxNumRowEntries))
    graph.haveGlobalConstants = true
end

function clearGlobalConstants(graph::CRSGraph)
    graph.globalNumEntries = 0
    graph.globalNumDiags = 0
    graph.globalMaxNumRowEntries = 0
    graph.haveGlobalConstants = false
end


function computeLocalConstants(graph::CRSGraph{GID, PID, LID}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    
    #short circuit if already computed
    graph.haveLocalConstants && return
    
    if @debug graph
        @assert !null(graph.colMap) "The graph must have a column map at this point"
    end
    
    #if graph.haveLocalConstants == false  #short circuited above
    
    graph.upperTriangle = true
    graph.lowerTriangle = true
    graph.nodeMaxNumRowEntries = 0
    graph.nodeNumDiags = 0
    
    rowMap = graph.rowMap
    colMap = get(graph.colMap)
    
    #indicesAreAllocated => true
    if  hasRowInfo(graph)
        const numLocalRows = numMyElements(rowMap)
        for localRow = LID(1):numLocalRows
            const globalRow = gid(rowMap, localRow)
            const rowLID = lid(colMap, globalRow)
            
            const rowInfo = getRowInfo(graph, localRow)
            rowView = getLocalView(rowInfo)

            if rowLID in rowView #appears to be broken IDK how or why
                graph.nodeNumDiags += 1
            end

            const smallestCol = rowView[1]
            const largestCol  = rowView[end]
            if smallestCol < localRow
                graph.upperTriangle = false
            end
            if localRow < largestCol
                graph.lowerTriangle = false
            end
            graph.nodeMaxNumRowEntries = max(graph.nodeMaxNumRowEntries, rowInfo.numEntries)
        end
    end
    graph.haveLocalConstants = true
end
              

hasRowInfo(graph::CRSGraph) = (getProfileType(graph) != STATIC_PROFILE 
                                || length(graph.rowOffsets) != 0)

function getRowInfoFromGlobalRow(graph::CRSGraph{GID, PID, LID},
        row::Integer)::RowInfo{LID} where {GID, PID, LID <: Integer}
    getRowInfo(graph, lid(graph.rowMap, row))
end
    
function getRowInfo(graph::CRSGraph{GID, PID, LID}, row::LID)::RowInfo{LID} where {GID, PID, LID <: Integer}
    if @debug graph
        @assert hasRowInfo(graph) "Graph does not have row info anymore.  Should have been caught earlier"
    end
    
    if !hasRowInfo(graph) || !myLID(graph.rowMap, row)
        return RowInfo{LID}(graph, 0, 0, 0, 0)
    end
    
    offset1D = 0
    allocSize = 0
    
    if getProfileType(graph) == STATIC_PROFILE
        if length(graph.rowOffsets) != 0
            offset1D  = graph.rowOffsets[row]
            allocSize = graph.rowOffsets[row+1] - graph.rowOffsets[row]
        end
        numEntries = (length(graph.numRowEntries) == 0 ?
            allocSize : graph.numRowEntries[row])
    else #dynamic profile
        if isLocallyIndexed(graph) && length(graph.localIndices2D) == 0
            allocSize = length(graph.localIndices2D[row])
            
        elseif isGloballyIndexed(graph) && length(graph.globalIndices2D) == 0
            allocSize = length(graph.globalIndices2D[row])
        end
        numEntries = (length(graph.numRowEntries) == 0 ?
            0 : graph.numRowEntries[row])
    end
    RowInfo{LID}(graph, row, allocSize, numEntries, offset1D)
end

function getLocalView(rowInfo::RowInfo{LID})::AbstractArray{LID, 1} where LID <: Integer
    graph = rowInfo.graph
    if rowInfo.allocSize == 0
        LID[]
    elseif length(graph.localIndices1D) != 0
        start = rowInfo.offset1D
        len = rowInfo.allocSize
        
        view(graph.localIndices1D, range(start, 1, len))
    elseif length(graph.localIndices2D[rowInfo.localRow]) != 0
        graph.localIndices2D[rowInfo.localRow]
    else
        LID[]
    end
end

function allocateIndices(graph::CRSGraph{GID, <:Integer, LID},
        lg::IndexType, numAllocPerRow::Array{<:Integer, 1}) where {
        GID <: Integer, LID <: Integer}
    numRows = getLocalNumRows(graph)
    @assert(length(numAllocPerRow) == numRows,
        "numAllocRows has length = $(length(numAllocPerRow)) "
        * "!= numRows = $numRows")
    allocateIndices(graph, lg, numAllocPerRow, i -> numAllocPerRow[i])
end

function allocateIndices(graph::CRSGraph{GID, <:Integer, LID},
        lg::IndexType, numAllocPerRow::Integer) where {
        GID <: Integer, LID <: Integer}
    allocateIndices(graph, lg, numAllocPerRow, i-> numAllocPerRow)
end
    
function allocateIndices(graph::CRSGraph{GID, <:Integer, LID},
        lg::IndexType, numAlloc, numAllocPerRow::Function) where {
        GID <: Integer, LID <: Integer}
    
    @assert(isLocallyIndexed(graph) == (lg == LOCAL_INDICES),
        "Graph is $(isLocallyIndexed(graph)?"":"not ")locally indexed, but lg=$lg")
    @assert(isGloballyIndexed(graph) == (lg == GLOBAL_INDICES),
        "Graph is $(isGloballyIndexed(graph)?"":"not ")globally indexed but lg=$lg")
    
    numRows = getLocalNumRows(graph)

    if getProfileType(graph) == STATIC_PROFILE
        rowPtrs = Array{LID, 1}(numRows + 1)
            
        computeOffsets(rowPtrs, numAlloc)
        
        graph.rowOffsets = rowPtrs
        numInds = rowPtrs[numRows+1]
        
        if lg == LOCAL_INDICES
            graph.localIndices1D = Array{LID, 1}(numInds)
        else
            graph.globalIndices1D = Array{GID, 1}(numInds)
        end
        graph.storageStatus = STORAGE_1D_UNPACKED
    else
        if lg == LOCAL_INDICES
            graph.localIndices2D = Array{Array{LID, 1}, 1}(numRows)
            for row = 1:numRows
                graph.localIndices2D[row] = Array{LID, 1}(numAllocPerRow(row))
            end
        else #lg == GLOBAL_INDICES
            graph.globalIndices2D = Array{Array{GID, 1}, 1}(numRows)
            for row = 1:numRows
                graph.globalIndices2D[row] = Array{GID, 1}(numAllocPerRow(row))
            end
        end
        graph.storageStatus = STORAGE_2D
    end

    graph.indicesType = lg
    
    if numRows > 0
        numRowEntries = zeros(LID, numRows)
        graph.numRowEntries = numRowEntries
    end
    
    #let the calling constructor take care of this
    #checkInternalState(graph)
end
    
    
function makeImportExport(graph::CRSGraph{GID, PID, LID}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    @assert !isnull(graph.colMap) "Cannot make imports and exports without a column map"
    
    if isnull(graph.importer)
        if get(graph.domainMap) !== get(graph.colMap) && !sameAs(get(graph.domainMap), get(graph.colMap))
            graph.importer = Import(graph.domainMap, graph.colMap, graph.plist)
        end
    end
    
    if isnull(graph.exporter)
        if get(graph.rangeMap) !== graph.rowMap && !sameAs(get(graph.rangeMap), graph.rowMap)
            graph.exporter = Export(graph.rowMap, graph.rangeMap, graph.plist)
        end
    end
end
    
#TODO migrate this to testing
function checkInternalState(graph::CRSGraph)
    if @debug graph
        const localNumRows = getLocalNumRows(graph)
        
        @assert(isFillActive(graph) != isFillComplete(graph),
            "Graph must be either fill active or fill "
            * "complete$(isFillActive(graph)?"not both":"").")
        @assert(!isFillComplete(graph)
                || (!isnull(graph.colMap) 
                    && !isnull(graph.domainMap) 
                    && !isnull(graph.rangeMap)),
            "Graph is fill complete, but at least one of {column, range, domain} map is null")
        @assert((graph.storageStatus != STORAGE_1D_PACKED 
                    && graph.storageStatus != STORAGE_1D_UNPACKED) 
            || graph.pftype != DYNAMIC_PROFILE,
            "Graph claims 1D storage, but dynamic profile")
        if graph.storageStatus == STORAGE_2D
            @assert(graph.pftype != STATIC_PROFILE ,
                "Graph claims 2D storage, but static profile")
            @assert(!isLocallyIndexed(graph) 
                || length(graph.localIndices2D) == localNumRows,
                "Graph calims to be locally index and have 2D storage, "
                * "but length(graph.localIndices2D) = $(length(graph.localIndices2D)) "
                * "!= getLocalNumRows(graph) = $localNumRows")
            @assert(!isGloballyIndexed(graph)
                || length(graph.globalIndices2D) == localNumRows,
                "Graph calims to be globally index and have 2D storage, "
                * "but length(graph.globalIndices2D) = $(length(graph.globalIndices2D)) "
                * "!= getLocalNumRows(graph) = $localNumRows")
        end
        
        @assert(graph.haveGlobalConstants 
            || (graph.globalNumEntries == 0 
                && graph.globalNumDiags == 0 
                && graph.globalMaxNumRowEntries == 0),
            "Graph claims to not have global constants, "
            * "but some of the global constants are not 0")
        
        @assert(!graph.haveGlobalConstants 
            || (graph.globalNumEntries != 0 
                && graph.globalMaxNumRowEntries != 0), 
            "Graph claims to have global constants, but also says 0 global entries")
        
        @assert(!graph.haveGlobalConstants
            || (graph.globalNumEntries > graph.nodeNumEntries
                && graph.globalNumDiags > graph.nodeNumDiags
                && graph.globalMaxNumRowEntries > graph.nodeMaxNumRowEntries),
            "Graph claims to have global constants, but some of the local "
            * "constants are greater than their corresponding global constants")
        
        @assert(!isStorageOptimized(graph)
            || graph.pftype == STATIC_PROFILE,
            "Storage is optimized, but graph is not STATIC_PROFILE")
        
        @assert(!isGloballyIndexed(graph)
            || length(graph.rowOffsets) == 0
            || (length(graph.rowOffsets) == localNumRows +1
                && graph.rowOffsets[localNumRows+1] == length(graph.globalIndices1D)),
            "If rowOffsets has nonzero size and the graph is globally "
            * "indexed, then rowOffsets must have N+1 rows and rowOffsets[N+1] "
            * "must equal the length of globalIndices1D")
        
        @assert(!isLocallyIndexed(graph)
            || length(graph.rowOffsets) == 0
            || (length(graph.rowOffsets) == localNumRows +1
                && graph.rowOffsets[localNumRows+1] == length(graph.localIndices1D)),
            "If rowOffsets has nonzero size and the graph is globally "
            * "indexed, then rowOffsets must have N+1 rows and rowOffsets[N+1] "
            * "must equal the length of localIndices1D")
        
        if graph.pftype == DYNAMIC_PROFILE
            @assert(localNumRows == 0 
                || length(graph.localIndices2D) > 0 
                || length(graph.globalIndices2D) > 0,
                "Graph has dynamic profile, the calling process has nonzero "
                * "rows, but no 2-D column index storage is present.")
            @assert(localNumRows == 0
                || length(graph.numRowEntries) != 0,
                "Graph has dynamic profiles and the calling process has "
                * "nonzero rows, but numRowEntries is not present")
            
            @assert(length(graph.localIndices1D) == 0
                && length(graph.globalIndices1D) == 0,
                "Graph has dynamic profile, but 1D allocations are present")
            
            @assert(length(graph.rowOffsets) == 0,
                "Graph has dynamic profile, but row offsets are present")
            
        elseif graph.pftype == STATIC_PROFILE
            @assert(length(graph.localIndices1D) != 0 
                || length(graph.globalIndices1D) != 0,
                "Graph has static profile, but 1D allocations are not present")
            
            @assert(length(graph.localIndices2D) == 0
                && length(graph.globalIndices2D) == 0,
                "Graph has static profile, but 2D allocations are present")
        else
            error("Unknown profile type: $(graph.pftype)")
        end
        
        if graph.indicesType == LOCAL_INDICES
            @assert(length(graph.globalIndices1D) == 0
                && length(graph.globalIndices2D) == 0,
                "Indices are local, but global allocations are present")
            @assert(graph.nodeNumEntries == 0
                || length(graph.localIndices1D) > 0
                || length(graph.localIndices2D) > 0,
                "Indices are local and local entries exist, but there aren't local allocations present")
        elseif graph.indicesType == GLOBAL_INDICES
            @assert(length(graph.localIndices1D) == 0
                && length(graph.localIndices2D) == 0,
                "Indices are global, but local allocations are present")
            @assert(graph.nodeNumEntries == 0
                || length(graph.globalIndices1D) > 0
                || length(graph.globalIndices2D) > 0,
                "Indices are global and local entries exist, but there aren't global allocations present")
        else
            warn("Unknown indices type: $(graph.indicesType)")
        end

        #check actual allocations
        const lenRowOffsets = length(graph.rowOffsets)
        if graph.pftype == STATIC_PROFILE && lenRowOffsets != 0
            @assert(lenRowOffsets == localNumRows+1,
                "Graph has static profile, rowOffsets has a nonzero length "
                * "($lenRowOffsets), but is not equal to the "
                * "local number of rows plus one ($(localNumRows+1))")
            const actualNumAllocated = graph.rowOffsets[localNumRows+1]
            @assert(!isLocallyIndexed(graph)
                || length(graph.localIndices1D) == actualNumAllocated,
                "Graph has static profile, rowOffsets has a nonzero length, "
                * "but length(localIndices1D) = $(length(graph.localIndices1D)) "
                * "!= actualNumAllocated = $actualNumAllocated")
            @assert(!isGloballyIndexed(graph)
                || length(graph.globalIndices1D) == actualNumAllocated,
                "Graph has static profile, rowOffsets has a nonzero length, "
                * "but length(globalIndices1D) = $(length(graph.globalIndices1D)) "
                * "!= actualNumAllocated = $actualNumAllocated")
        end
    end
end

function setLocallyModified(graph::CRSGraph)
    graph.indicesAreSorted = false
    graph.noRedundancies = false
    graph.haveLocalConstants = false
end

function sortAndMergeAllIndices(graph::CRSGraph, sorted::Bool, merged::Bool)
    @assert(isLocallyIndexed(graph),
        "This method may only be called after makeIndicesLocal(graph)")
    @assert(merged || isStoragedOptimized(graph),
        "The graph is already storage optimized, "
        * "so we shouldn't be merging any indices.")
    
    if !sorted || !merged
        localNumRows = getLocalNumRows(graph)
        totalNumDups = 0
        for localRow = 1:localNumRows
            rowInfo = getRowInfo(graph, localRow)
            if !sorted
                sortRowIndices(graph, rowInfo)
            end
            if !merged
                numDups += mergeRowIndices(graph, rowInfo)
            end
        end
        graph.nodeNumEntries -= totalNumDups
        graph.indiciesAreSorted = true
        graph.noRedunancies = true
    end
end

function sortRowIndices(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}) where {GID, PID, LID <: Integer}
    if rowInfo.numEntries > 0
        localColumnIndices = getLocalView(graph, rowInfo)
        sort!(localColumnIndices)
    end
end

function mergeRowIndices(graph::CRSGraph{GID, PID, LID}, rowInfo::RowInfo{LID}) where {GID, PID, LID <: Integer}
    localColIndices = getLocalView(graph, rowInfo)
    localColIndices[:] = unique(localColIndices)
    mergedEntries = length(localColIndices)
    graph.numRowEntries[rowInfo.localRow] = mergedEntries
    
    rowInfo.numEntries - mergedEntries
end
        

function setDomainRangeMaps(graph::CRSGraph{GID, PID, LID}, domainMap::BlockMap{GID, PID, LID}, rangeMap::BlockMap{GID, PID, LID}) where {GID, PID, LID}
    if graph.domainMap != domainMap
        graph.domainMap = domainMap
        graph.importer = Nullable{Import{GID, PID, LID}}()
    end
    if graph.rangeMap != rangeMap
        graph.rangeMap = rangeMap
        graph.exporter = Nullable{Export{GID, PID, LID}}()
    end
end


function globalAssemble(graph::CRSGraph)
    @assert isFillActive(graph) "Fill must be active before calling globalAssemble(graph)"
    
    comm = julia_petra.comm(graph)
    myNumNonlocalRows = length(graph.nonlocals)
    
    maxNonlocalRows = maxAll(comm, myNumNonlocalRows)
    if maxNonlocalRows != 0
        return
    end
    
    #skipping: nonlocalRowMap = null
    
    numEntPerNonlocalRow = Array{LID, 1}(myNumNonlocalRows)
    myNonlocalGlobalRows = Array{GID, 1}(myNumNonlocalRows)
    
    for (i, (key, val)) = zip(1:length(graph.nonlocals), graph.nonlocals)
        myNonlocalGlobalRows[i] = key
        const globalCols = val #const b/c changing in place
        sort!(globalCols)
        globalCols[:] = unique(globalCols)
        numEntPerNonlocalRow[i] = length(globalCols)
    end
    
    myMinNonLocalGlobalRow = minimum(myNonLocalGlobalRows)
    
    globalMinNonlocalRow = minAll(comm, myMinNonlocalGlobalRow)
    
    nonlocalRowMap = BlockMap(-1, myNonlocalGlobalRows, comm)
    
    nonlocalGraph = CRSGraph(nonlocalRowMap, numEntPerNonlocalRow, STATIC_PROFILE)
    for (i, (key, val)) = zip(1:length(graph.nonlocals), graph.nonlocals)
        globalRow = key
        globalColumns = val
        numEnt = length(numEntPerNonlocalRow[i])
        insertGlobalIndices(nonLocalGraph, globalRow, numEnt, globalColumns)
    end
    
    const origRowMap = graph.rowMap
    const origRowMapIsOneToOne = isOneToOne(origRowMap)
    
    if origRowMapIsOneToOne
        exportToOrig = Export(nonlocalRowMap, origRowMap)
        doExport(nonLocalGraph, graph, exportToOrig, INSERT)
    else
        oneToOneRowMap = createOneToOne(origRowMap)
        exportToOneToOne = Export(nonlocalRowMap, oneToOneRowMap)
        
        oneToOneGraph = CRSGraph(oneToOneRowMap, 0)
        doExport(nonlocalGraph, oneToOneGraph, exportToOneToOne, INSERT)
        
        #keep memory highwater mark down
        #nonlocalGraph = null
        
        importToOrig(oneToOneRowMap, origRowMap)
        doImport(oneToOneGraph, graph, importToOrig, INSERT)
    end
    clear!(graph.nonLocals)
    
    checkInternalState(graph)
end
     
function makeIndicesLocal(graph::CRSGraph{GID, PID, LID}) where {GID, PID, LID}
    @assert hasColMap(graph) "The graph does not have a column map yet.  This method should never be called in that case"
    
    colMap = get(graph.colMap)
    localNumRows = getLocalNumRows(graph)
    
    if isGloballyIndexed(graph) && localNumRows != 0
        numRowEntries = graph.numRowEntries
        
        if getProfileType(graph) == STATIC_PROFILE
            if GID == LID
                graph.localIndices1D = graph.globalIndices1D

                println("Reusing global index array")
            else
                @assert(length(graph.rowOffsets) != 0,
                    "length(graph.rowOffsets) == 0.  "
                    * "This should never happen at this point")
                const numEnt = graph.rowOffsets[localNumRows+1]
                graph.localIndices1D = Array{LID, 1}(numEnt)
            end
            
            localColumnMap = getLocalMap(colMap)
            
            numBad = convertColumnIndicesFromGlobalToLocal(
                        graph.localIndices1D,
                        graph.globalIndices1D,
                        graph.rowOffsets,
                        localColumnMap,
                        numRowEntries)
            
            if numBad != 0
                throw(InvalidArgumentError("When converting column indices from "
                        * "global to local, we enchoundered $numBad indices that "
                        * "do not live in the column map on this process"))
            end
            
            graph.globalIndices1D = Array{LID, 1}(0)
        else #graph has dynamic profile
            graph.localIndices2D = Array{Array{LID, 1}, 1}(localNumRows)
            for localRow = 1:localNumRows
                if length(graph.globalIndices2D[localRow]) != 0
                    globalIndices = graph.globalIndices2D[localRow]
                    
                    graph.localIndices2D[localRow] = [lid(colMap, gid) for gid in globalIndices]
                    if @debug graph
                        @assert(minimum(graph.localIndices2D[localRow]) > 0,
                            "Globalal indices were not found in the column Map")
                    end
                end
            end
            graph.globalIndices2D = Array{GID, 1}[]
        end
    end
    
    graph.localGraph = LocalCRSGraph(graph.localIndices1D, graph.rowOffsets)
    graph.indicesType = LOCAL_INDICES
    checkInternalState(graph)
end


function convertColumnIndicesFromGlobalToLocal(localColumnIndices::Array{LID, 1},
        globalColumnIndices::Array{GID, 1}, ptr::Array{LID, 1},
        localColumnMap::BlockMap{GID, PID, LID}, numRowEntries::Array{LID, 1})::LID where {
        GID, PID, LID}

    localNumRows = max(length(ptr)-1, 0)
    numBad = 0
    for localRow = 1:localNumRows
        offset = ptr[localRow]
        
        for j = 0:numRowEntries[localRow]-1
            gid = globalColumnIndices[offset+j]
            localColumnIndices[offset+j] = lid(localColumnMap, gid)
            if lid == 0
                numBad += 1
            end
        end
    end
    numBad
end

#covers the overlap between insert methods
macro insertIndicesImpl(indicesType, innards)
    indices1D = Symbol(indicesType*"Indices1D")
    indices2D = Symbol(indicesType*"Indices2D")

    esc(quote
            rowInfo = getRowInfo(graph, myRow)
            numNewIndices = length(indices)
            newNumEntries = rowInfo.numEntries + numNewIndices
            
            if newNumEntries > rowInfo.allocSize
                if getProfileType(graph) == STATIC_PROFILE
                    $innards
                else
                    newAllocSize = 2*rowInfo.allocSize
                    if newAllocSize < newNumEntries
                        newAllocSize = newNumEntries
                    end
                    resize!(graph.$indices2D[myRow], newAllocSize)
                end
            end

            if length(graph.$indices1D) != 0
                offset = rowInfo.offset1D + rowInfo.numEntries
                destRange = offset+1:offset+numNewIndices

                graph.$indices1D[destRange] = indices[:]
            else
                graph.$indices2D[myRow][rowInfo.numEntries+1:newNumEntries] = indices[:]
            end

            graph.numRowEntries[myRow] += numNewIndices
            graph.nodeNumEntries += numNewIndices
            setLocallyModified(graph)

            if @debug graph
                chkNewNumEntries = getNumEntriesInLocalRow(graph, myRow)
                @assert(chkNewNumEntries == newNumEntries,
                    "Internal Logic error: chkNewNumEntries = $chkNewNumEntries "
                    * "!= newNumEntries = $newNumEntries")
            end
    end)
end

function insertLocalIndicesImpl(graph::CRSGraph{GID, PID, LID},
        myRow::LID, indices::AbstractArray{LID, 1}) where {
        GID, PID, LID <: Integer}
    @insertIndicesImpl "local" begin  
        throw(InvalidArgumentError("new indices exceed statically allocated graph structure"))
    end
end

#TODO figure out if this all can be moved to @insertIndicesImpl
function insertGlobalIndicesImpl(graph::CRSGraph{GID, PID, LID},
        myRow::LID, indices::AbstractArray{GID, 1}) where {
        GID <: Integer, PID, LID <: Integer}
    @insertIndicesImpl "global" begin
        @assert(rowInfo.numEntries <= rowInfo.allocSize,
            "For local row $myRow, rowInfo.numEntries = $(rowInfo.numEntries) "
            * "> rowInfo.allocSize = $(rowInfo.allocSize).")

        dupCount = 0
        if length(graph.globalIndices1D) != 0
            curOffset = rowInfo.offset1D
            @assert(length(graph.globalIndices1D) >= curOffset,
                "length(graph.globalIndices1D) = $(length(graph.globalIndices1D)) "
                * ">= curOffset = $curOffset")
            @assert(length(graph.globalIndices1D) >= curOffset + rowInfo.offset1D,
                "length(graph.globalIndices1D) = $(length(graph.globalIndices1D)) "
                * ">= curOffset+rowInfo.offset1D = $(curOffset + rowInfo.offset1D)")

            range = curOffset:curOffset+rowInfo.numEntries
            globalIndicesCur = view(graph.globalIndices1D, range)
        else
            #line 1959

            globalIndices = graph.globalIndices2D[myRow]
            @assert(rowInfo.allocSize == length(globalIndices),
                "rowInfo.allocSize = $(rowInfo.allocSize) "
                * "== length(globalIndices) = $(length(globalIndices))")
            @assert(rowInfo.numEntries <= length(globalIndices), 
                "rowInfo.numEntries = $(rowInfo.numEntries) "
                * "== length(globalIndices) = $(length(globalIndices))")

            globalIndicesCur = view(globalIndices, 0, rowInfo.numEntries)
        end
        for newIndex = indices
            dupCount += count(old -> old==newIndex, globalIndicesCur)
        end

        numNewToInsert = numNewInds - dupCount
        @assert numNewToInsert >= 0 "More duplications than indices"

        if rowInfo.numEntries + numNewToInsert > rowInfo.allocSize
            throw(InvalidArgumentError("$(myPid(comm(graph))): "
                    * "For local row $myRow, even after excluding "
                    * "$dupCount duplicate(s) in input, the new number "
                    * "of entries $(rowInfo.numEntries + numNewToInsert) "
                    * "still exceeds this row's static allocation size "
                    * "$(rowInfo.allocSize).  You must either fix the upper "
                    * "bound on number of entries in this row, or switch "
                    * "to dynamic profile."))
        end

        if length(graph.globalIndices) != 0
            curOffset = rowInfo.offset1D
            globalIndicesCur = view(graph.globalIndices1D,
                range(curOffset, 1, rowInfo.numEntries))
            globalIndicesNew = view(graph.globalIndices1D,
                curOffset+rowInfo.numEntries+1 : currOffset+rowInfo.allocSize)
        else
            #line 2036

            globalIndices = graph.globalIndices2D[myRow]
            globalIndicesCur = view(globalIndices, 1:rowInfo.numEntries)
            globalIndicesNew = view(globalIndices,
                rowInfo.numEntries+1 : rowInfo.allocSize-rowInfo.numEntries)
        end

        curPos = 1
        for globalIndexToInsert = indices

            alreadyInOld = globalIndexToInsert in globalIndicesCur
            if !alreadyInOld
                @assert(curPos <= numNewToInsert,
                    "curPos = $curPos >= numNewToInsert = $newToInsert.")
                globalIndicesNew[curPos] = globalIndexToInsert
                curPos += 1
            end
        end

        graph.numRowEntries[myRow] = rowInfo.numEntries+numNewToInsert
        graph.nodeNumEntries += numNewToInsert
        setLocallyModified(graph)

        if @debug graph
            newNumEntries = rowInfo.numEntries + numNewToInsert
            chkNewNumEntries = getNumEntiresInLocalRow(graph, myRow)
            @assert(chkNewNumEntries == newNumEntries,
                "chkNewNumEntries = $chkNewNumEntries "
                * "!= newNumEntries = $newNumEntries")
        end
        return
    end
end



#internal implementation of makeColMap, needed to handle some return and debuging stuff
#returns Tuple(errCode, colMap)
function __makeColMap(graph::CRSGraph{GID, PID, LID}, wrappedDomMap::Nullable{BlockMap{GID, PID, LID}},
        sortEachProcsGIDs::Bool) where {GID, PID, LID}
    errCode = 0#TODO improve from int error code

    if isnull(wrappedDomMap)
        colMap = Nullable{BlockMap{GID, PID, LID}}()
    else
        myColumns = GID[]
        domMap = get(wrappedDomMap)

        if isLocallyIndexed(graph)
            wrappedColMap = graph.colMap

            if isnull(wrappedColMap)
                warn("$(myPid(comm(graph))): The graph is locally indexed, but does not have a column map")

                errCode = -1
            else
                colMap = get(wrappedColMap)
                if linearMap(colMap) #i think isContiguous(map) <=> linearMap(map)?
                    numCurGIDs = numMyElements(colMap)
                    myFirstGlobalIndex = minMyGIDs(colMap)
                    
                    myColumns = collect(range(myFirstGlobalIndex, 1, numCurGIDs))
                else
                    myColumns = copy(myGlobalElements(colMap))
                end
            end
        else #if graph.isGloballyIndexed
            const localNumRows = getLocalNumRows(graph)

            numLocalColGIDs = 0
            numRemoteColGIDs = 0

            gidIsLocal = zeros(Bool, localNumRows)
            remoteGIDSet = Set()
            remoteGIDUnorderedVector = GID[]

            #if rowMap != null
            const rowMap = graph.rowMap

            for localRow = 1:localNumRows
                globalRow = gid(rowMap, localRow)
                rowGIDs = getGlobalRowView(graph, globalRow)

                numEnt = length(rowGIDs)
                if numEnt != 0
                    for k = 1:numEnt
                        gid = rowGIDs[k]
                        lid = julia_petra.lid(domMap, gid)
                        if lid != 0
                            if !gidIsLocal[lid]
                                gidIsLocal[lid] = true
                                numLocalColGIDs += 1
                            end
                        else
                            if !in(remoteGIDSet, gid)
                                push!(remoteGIDSet, gid)
                                if !sortEachProcsGIDs
                                    #user wants order retained
                                    push!(remoteGIDUnorderedVector, gid)
                                end
                                numRemoteColGIDs += 1
                            end
                        end
                    end
                end
            end
        end

        #line 214, abunch of explanation of serial short circuit
        if numProc(comm(domMap)) == 1
            if numRemoteColGIDs != 0
                errCode = -2
            end
            if numLocalColGIDs == localNumRows
                return (errCode, domMap)
            end
            resize!(myColumns, numLocalColGIDs+numRemoteColGIDs)
            localColGIDs  = view(myColumns, 1:numLocalColGIDs)
            remoteColGIDs = view(myColumns, numLocalColGIDs+1:numRemoteColGIDs)

            if sortEachProcsGIDs
                remoteColGIDs[:] = [el for el in remoteGIDSet]
            else
                remoteColGIDs[:] = remoteGIDUnorderedVector
            end

            remotePIDs = Array{PID, 1}(numRemoteColGIDs)

            remotePIDs = remoteIDList(domMap, remoteColGIDs)
            if any(remotePIDs .== 0)
                if @debug graph
                    warn("Some column indices are not in the domain Map")
                end
                errCode = 3
            end

            order = sortperm(remotePIDs)
            permute!(remotePIDs, order)
            permute!(remoteColGIDs, order)

            #line 333

            numDomainELts = numMyElements(domMap)
            if numLocalColGIDs == numDomainElts
                if linearMap(domMap) #I think isContiguous() <=> linearMap()
                    localColGIDs[1:numLocalColGIDs] = minMyGIDs(domMap)
                else
                    domElts = myGlobalElements(domMap)
                    localColGIDs[1:length(domElts)] = domElts
                end
            else
                numLocalCount = 0
                if linearMap(domMap) #I think isContiguous() <=> linearMap()
                    curColMapGID = minMyGIDs(domMap)
                    for i = 1:numDomainElts
                        if gidIsLocal[i]
                            localColGIDs[numLocalCount] = curColMapGID
                            numLocalCount += 1
                        end
                        curColMapGID += 1
                    end
                else
                    domainElts = myGlobalElement(domMap)
                    for i = 1:numDomainElts
                        if gidIsLocal[i]
                            localColGIDs[numLocalCount] = domainElts[i]
                            numLocalCount += 1
                        end
                        curColMapGID += 1
                    end
                end

                if numLocalCount != numLocalColGIDs
                    if @debug graph
                        warn("$(myPid(comm(graph))): numLocalCount = $numLocalCount "
                            * "!= numLocalColGIDs = $numLocalColGIDs.  "
                            * "This should not happen.")
                    end
                    errCode = -4
                end
            end
        end
    end
    return(errCode, BlockMap(-1, -1, myColumns, comm(domMap)))
end
