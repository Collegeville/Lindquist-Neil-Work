
export CRSGraph, isLocallyIndexed, isGloballyIndexed, getProfileType

#TODO document this type and it's methods
#TODO ensure graphs that lack a colMap use global indexing

mutable struct CRSGraph{GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::Nullable{BlockMap{GID, PID, LID}}
    rangeMap::Nullable{BlockMap{GID, PID, LID}}
    domainMap::Nullable{BlockMap{GID, PID, LID}}
    
    #may be null if domainMap and colMap are the same
    importer::Nullable{Import{GID, PID, LID}}
    #may be null if rangeMap and rowMap are the same
    exporter::Nullable{Export{GID, PID, LID}}

    #TODO figure out if this really needs to be stored
    lclGraph::LocalCRSGraph

    #Local number of (populated) entries; must always be consistent
    nodeNumEntries::LID

    #Local number of (populated) diagonal entries.
    nodeNumDiags::LID

    #Local maximum of the number of entries in each row.
    nodeMaxNumRowEntries::LID

    #Global number of entries in the graph.
    globalNumEntries::GID

    #Global number of (populated) diagonal entries.
    globalNumDiags::GID

    #Global maximum of the number of entries in each row.
    globalMaxNumRowEntries::GID

    #Whether the graph was allocated with static or dynamic profile.
    pftype::ProfileType


    #TODO combine 1-D and 2-D storage

    ## 1-D storage (Static profile) data structures ##
    localIndices1D::Array{LID, 1}
    globalIndices1D::Array{GID, 1}
    rowOffsets::Array{LID, 1}  #Tpetra: k_rowPts_

    ## 2-D storage (Dynamic profile) data structures ##
    localIndices2D::Array{Array{LID, 1}, 1}
    globalIndices2D::Array{Array{GID, 1}, 1}
    #may exist in 1-D storage if not packed
    numRowEntries::Array{LID, 1}

    storageStatus::StorageStatus

    indicesAllowed::Bool
    indicesType::IndexType
    fillComplete::Bool

    lowerTriangle::Bool
    upperTriangle::Bool
    indicesAreSorted::Bool
    noRedundancies::Bool
    haveLocalConstants::Bool
    haveGlobalConstants::Bool
    sortGhostsAssociatedWithEachProcessor::Bool

    plist::Dict{Symbol}

    nonLocals::Dict{GID, Array{GID, 1}}

    #Large ammounts of duplication between the constructors, so group it in an inner constructor
    function CRSGraph(
        rowMap::BlockMap{GID, PID, LID},
        colMap::Nullable{BlockMap{GID, PID, LID}},
        rangeMap::Nullable{BlockMap{GID, PID, LID}},
        domainMap::Nullable{BlockMap{GID, PID, LID}},

        lclGraph::LocalCRSGraph,

        nodeNumEntries::LID,

        pftype::ProfileType,
        storageStatus::StorageStatus,

        indiciesType::IndexType,
        plist::Dict{Symbol}
    ) where {GID <: Integer, PID <: Integer, LID <: Integer}

        graph = new{GID, PID, LID}(
            rowMap,
            colMap,
            rangeMap,
            domainMap,

            Nullable{Import{GID, PID, LID}}(),
            Nullable{Export{GID, PID, LID}}(),

            lclGraph,

            #Local number of (populated) entries; must always be consistent
            nodeNumEntries,

            #using 0 to indicate uninitiallized, since -1 isn't gareenteed to work
            0, #nodeNumDiags
            0, #nodeMaxNumRowEntries
            0, #globalNumEntries
            0, #globalNumDiags
            0, #globalMaxNumRowEntries

            #Whether the graph was allocated with static or dynamic profile.
            pftype,


            ## 1-D storage (Static profile) data structures ##
            LID[],
            GID[],
            LID[],

            ## 2-D storage (Dynamic profile) data structures ##
            Array{Array{LID, 1}, 1}(0),
            Array{Array{GID, 1}, 1}(0),
            LID[],

            storageStatus,

            false,
            indiciesType,
            false,

            false,
            false,
            true,
            true,
            false,
            false,
            true,

            plist,
        
            Dict{GID, Array{GID, 1}}()
        )

        ## staticAssertions() 
        #skipping sizeof checks
        #skipping max value checks related to size_t

        #ensure LID is a subset of GID (for positive numbers)
        if !(LID <: GID) && (GID != BigInt) && (GID != Integer)
            # all ints are assumed to be able to handle 1, up to their max
            if LID == BigInt || LID == Integer || typemax(LID) > typemax(GID)
                throw(InvalidArgumentError("The positive values of GID must "
                        * "be a superset of the positive values of LID"))
            end
        end

        graph
    end
end


#### Constructors #####

#TODO add plist constructors that passes kwargs
function CRSGraph(rowMap::BlockMap{GID, PID, LID}, maxNumEntriesPerRow::LID,
        pftype::ProfileType, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    CRSGraph(rowMap, Nullable{BlockMap{GID, PID, LID}}(), maxNumEntriesPerRow, pftype, plist)
end
function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        maxNumEntriesPerRow::LID, pftype::ProfileType, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    CRSGraph(rowMap, Nullable(colMap), maxNumEntriesPerRow, pftype, plist)
end
    
function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::Nullable{BlockMap{GID, PID, LID}},
        maxNumEntriesPerRow::LID, pftype::ProfileType, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    graph = CRSGraph(
        rowMap,
        colMap,
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),

        LocalCRSGraph{LID, LID}(), #lclGraph
        
        LID(0), #nodeNumEntries

        pftype,

        (pftype == STATIC_PROFILE ?
              STORAGE_1D_UNPACKED 
            : STORAGE_2D),
        
        LOCAL_INDICES,
        plist
    )
        
    allocateIndices(graph, LOCAL_INDICES, maxNumEntriesPerRow)
    
    resumeFill(graph, plist)
    checkInternalState(graph)
    
    graph
end


function CRSGraph(rowMap::BlockMap{GID, PID, LID}, numEntPerRow::Array{LID, 1},
        pftype::ProfileType, plist::Dict{Symbol})  where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    CRSGraph(rowMap, Nullable{BlockMap{GID, PID, LID}}(), numEntPerRow, pftype, plist)
end

function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        numEntPerRow::Array{LID, 1}, pftype::ProfileType,
        plist::Dict{Symbol})  where {GID <: Integer, PID <: Integer, LID <: Integer}
    CRSGraph(rowMap, Nullable(colMap), numEntPerRow, pftype, plist)
end

function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::Nullable{BlockMap{GID, PID, LID}},
        numEntPerRow::Array{LID, 1}, pftype::ProfileType,
        plist::Dict{Symbol})  where {GID <: Integer, PID <: Integer, LID <: Integer}
    debug = get(plist, :debug, false)
    graph = CRSGraph(
        rowMap,
        colMap,
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),

        LocalCRSGraph{LID, LID}(), #lclGraph
        
        LID(0), #nodeNumEntries

        #Whether the graph was allocated with static or dynamic profile.
        pftype,
        
        (pftype == STATIC_PROFILE ?
              STORAGE_1D_UNPACKED 
            : STORAGE_2D),

        LOCAL_INDICES,
        plist
    )
    
    lclNumRows = numLocalElements(rowMap)
    if length(numEntPerRow) != lclNumRows
        throw(InvalidArgumentError("numEntPerRows has length $(length(numEntPerRow)) " *
                "!= the local number of rows $lclNumRows as spcified by the input row Map"))
    end
    
    if debug
        for r = 1:lclNumRows
            curRowCount = numEntPerRow[r]
            if curRowCount <= 0
                throw(InvalidArgumentError("numEntPerRow[$r] = $curRowCount is not valid"))
            end
        end
    end
    
    allocateIndices(graph, LOCAL_INDICES, numEntPerRow)
    
    resumeFill(graph, plist)
    checkInternalState(graph)
    
    graph
end


function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        rowPointers::Array{LID, 1}, columnIndices::Array{LID, 1},
        plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    graph = CRSGraph(
        rowMap,
        Nullable(colMap),
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),

        LocalCRSGraph{LID, LID}(), #lclGraph
        
        LID(0), #nodeNumEntries
        
        STATIC_PROFILE,
        
        STORAGE_1D_PACKED,

        LOCAL_INDICES,
        plist
    )
    #seems to be already taken care of
    #allocateIndices(graph, LOCAL_INDICES, numEntPerRow)
    
    setAllIndicies(graph, rowPointers, columnIndicies)
    checkInternalState(graph)
    
    graph
end


function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localGraph::LocalCRSGraph{LID, LID}, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    mapRowCount = numMyElements(rowMap)
    graph = CRSGraph(
        rowMap,
        Nullable(colMap),
        Nullable{}(rowMap),
        Nullable{}(colMap),
        
        localGraph,
        
        localGraph.rowMap[mapRowCount+1], #nodeNumEntries
        
        STATIC_PROFILE,
        
        STORAGE_1D_PACKED,
        
        LOCAL_INDICES,
        plist
    )
    
    if numRows(localGraph) != numMyElements(rowMap)
        throw(InvalidArgumentError("input row map and input local "
                * "graph need to have the same number of rows.  The "
                * "row map claims $(numMyElements(rowMap)) row(s), "
                * "but the local graph claims $(numRows(localGraph)) "
                * "row(s)."))
    end
    
    #seems to be already taken care of
    #allocateIndices(graph, LOCAL_INDICES, numEntPerRow)
    
    makeImportExport(graph)
    
    d_inds = localGraph.entries
    graph.localIndices1D = d_inds
    
    d_ptrs = localGraph.rowMap
    graph.rowOffsets = d_ptrs
    
    #reset local properties
    graph.upperTriangle = true
    graph.lowerTriangle = true
    graph.nodeMaxNumRowEntries = 0
    graph.nodeNumDiags
    
    for localRow = 1:mapRowCount
        globalRow = gid(rowMap, localRow)
        rowLID = lid(colMap, globalRow)
        
        #possible that the local matrix has no entries in the column
        #corrisponding to the current row, in that case, the column map
        #might not contain that GID.  Hence, the index validity check
        if rowLID != 0
            if rowLID +1 > length(d_ptrs)
                throw(InvalidArgumentError("The given row Map and/or column Map "
                        * "is/are not compatible with the provided local graphs."))
            end
            if d_ptrs[rowLID] != d_ptr[rowLID_1]
                const smallestCol = d_inds[d_ptrs[rowLID]]
                const largestCol  = d_inds[d_ptrs[rowLID+1]-1]
                
                if smallestCol < localRow
                    graph.upperTriangle = false
                end
                if localRow < largestCol
                    graph.lowerTriangle = false
                end
                for i = d_ptrs[rowLID]:d_ptrs[rowLID]-1
                    if d_inds[i] == rowLID
                        graph.nodeNumDiags += 1
                        break #can only be 1 diagonal per row
                    end
                end
            end
            
            graph.nodeMaxNumRowEntries = max((d_ptrs[rowLID + 1] - d_ptrs[rowLID]),
                                            graph.nodeMaxNumRowEntries)
        end
    end
    
    graph.hasLocalConstants = true
    computeGlobalConstants(graph)
    
    graph.fillComplete = true
    checkInternalState(graph)
    
    graph
end
    
    
#### internal methods ####
function computeGlobalConstants(graph::CRSGraph{GID, PID, LID}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    
    #short circuit if already computed
    graph.haveGlobalConstants && return
    
    if get(graph.plist, :debug, false)
        @assert !null(graph.colMap) "The graph must have a column map at this point"
    end
    
    computeLocalConstants()
    
    #if graph.haveGlobalConstants == false  #short circuited above
    graph.globalNumEntries, graph.globalNumDiags = sumAll(comm(graph.map),
        [GID(graph.nodeNumEntries), GID(graph.nodeNumDiags)])
        
    graph.globalMaxNumRowEntries = maxAll(comm(graph.map), GID(graph.nodeMaxNumRowEntries))
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
    
    if get(graph.plist, :debug, false)
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
        for localRow = 1:numLocalRows
            const globalRow = gid(rowMap, localRow)
            const rowLID = lid(colMap, globalRow)
            
            const rowInfo = getRowInfo(graph, localRow)
            rowView = getLocalView(rowInfo)

            for i = rowView
                if rowLID == i
                    graph.nodeNumDiags += 1
                    break #should only have 1 diag per row
                end
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
              

function hasRowInfo(graph::CRSGraph)
    #indicesAreAllocated => true
    (getProfileType(graph) != STATIC_PROFILE
        || length(graph.rowOffsets) != 0)
end

#implement getRowInfoFromGlobalRow
function getRowInfo(graph::CRSGraph{GID, PID, LID}, row::LID)::RowInfo{LID} where {GID, PID, LID <: Integer}
    if get(graph.plist, :debug, false)
        @assert hasRowInfo(graph) "Graph does not have row info anymore.  Should have been caught earlier"
    end
    
    if !hasRowInfo(graph) || !myLID(graph.rowMap, row)
        return RowInfo{LID}(graph, 0, 0, 0, 0)
    end
    
    if getProfileType(graph) == STATIC_PROFILE
        if length(graph.rowOffsets) == 0
            ofset1D = 0
            allocSize = 0
        else
            offset1D  = graph.rowOffsets[row]
            allocSize = graph.rowOffsets[row+1] - graph.rowOffsets[row]
        end
        numEntries = (length(graph.numRowEntries) == 0 ?
            allocSize : numRowEntries[row])
    else #dynamic profile
        offset1D = 0
        if isLocallyIndexed(graph)
            allocSize = (length(graph.localIndices2D) == 0 ?
                0 : length(graph.localIndices2D[row]))
        elseif isGloballyIndexed(graph)
            allocSize = (length(graph.globalIndices2D) == 0 ?
                0 : length(graph.globalIndices2D[row]))
        else
            allocSize = 0
        end
        numEntries = (length(graph.numRowEntries) == 0 ?
            0 : graph.numRowEntries[row])
    end
    RowInfo{LID}(graph, row, allocSize, numEntries, offset1D)
end

#DECISION put this somewhere else?  Its only an internal grouping
struct RowInfo{LID <: Integer}
    graph::CRSGraph{<:Integer, <:Integer, LID}
    localRow::LID
    allocSize::LID
    numEntries::LID
    offset1D::LID
end

function getLocalView(rowInfo::RowInfo{LID})::AbstractArray{LID, 1} where LID <: Integer
    graph = rowInfo.graph
    if rowInfo.allocSize == 0
        LID[]
    elseif length(graph.localIndices1D) != 0
        start = rowInfo.offset1D
        len = rowInfo.allocSize
        
        view(graph.localIndices1D, start:len)
    elseif length(graph.localIndices2D[rowInfo.localRow]) != 0
        graph.localIndices2D[rowInfo.localRow]
    else
        LID[]
    end
end

function allocateIndices(graph::CRSGraph{GID, <:Integer, LID},
        lg::IndexType, numAllocPerRow::Array{<:Integer, 1}) where {
        GID <: Integer, LID <: Integer}
    #using function call style, since otherwise it needs to be 1 line or julia funcalls a Bool
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
    @assert isLocallyIndexed(graph) == (lg == LOCAL_INDICES) "Graph is locally indexed, but lg=$lg"
    @assert isGloballyIndexed(graph) == (lg == GLOBAL_INDICES) "Graph is globally indexed but lg=$lg"
    
    numRows = getNodeNumRows(graph)

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
            graph.localInds2D = Array{Array{LID, 1}, 1}(numRows)
            for i = 1:numRows
                const howMany = numAllocPerRow(i)
                if howMany > 0
                    resize(graph.localInds2D[i], howMany)
                end
            end
        else #lg == GLOBAL_INDICES
            graph.globalInds2D = Array{Array{GID, 1}, 1}(numRows)
            for i = 1:numRows
                const howMany = numAllocPerRow(i)
                if howMany > 0
                    resize(graph.globalInds2D[i], howMany)
                end
            end
        end
        graph.storageStatus = STORAGE_2D
    end

    graph.indicesType = lg
    
    if numRows > 0
        numRowEntries = zeros(LID, 1)
        graph.numRowEntries = numRowEntries
    end
    
    #let the calling constructor take care of this
    #checkInternalState(graph)
end
    
    
function makeImportExport(graph::CRSGraph{GID, PID, LID}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    @assert !null(graph.colMap) "Cannot make imports and exports without a column map"
    
    if null(graph.importer)
        if graph.domainMap != graph.colMap && !isSameAs(graph.domainMap, graph.colMap)
            graph.importer = Import(graph.domainMap, graph.colMap, graph.plist)
        end
    end
    
    if null(graph.exporter)
        if graph.rangeMap != graph.rowMap && !isSameAs(graph.rangeMap, graph.rowMap)
            graph.exporter = Export(graph.rowMap, graph.rangeMap, graph.plist)
        end
    end
end
    
function checkInternalState(graph::CRSGraph)
    if get(graph.plist, :debug, false)
        const localNumRows = getNodeNumRows(graph)
        
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
                * "!= getNodeNumRows(graph) = $localNumRows")
            @assert(!isGloballyIndexed(graph)
                || length(graph.globalIndices2D) == localNumRows,
                "Graph calims to be globally index and have 2D storage, "
                * "but length(graph.globalIndices2D) = $(length(graph.globalIndices2D)) "
                * "!= getNodeNumRows(graph) = $localNumRows")
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
                && length(graph.globalIndices2D) == 0,
                "Indices are global, but local allocations are present")
            @assert(graph.localNumEntries == 0
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


#TODO implement globalAssemble(::CRSGraph)
#TODO implement setDomainRangeMaps(graph, domainMap, rangeMap)
#TODO implement makeIndicesLocal(graph)
#TODO implement sortAndMergeAllIndices(graph, isSorted, isMerged)

#### external API ####

function getNumEntriesInGlobalRow(graph::CRSGraph{GID}, globalRow::GID)::Integer where {GID <: Integer}
    localRow = lid(graph.rowMap, globalRow)
    if hasRowInfo(graph) && localRow != 0
        getRowInfo(localRow).numEntries
    else
        -1
    end
end

function getNumEntriesInLocalRow(graph::CRSGraph{GID, PID, LID}, localRow::LID)::Integer where {GID, PID, LID <: Integer}
    if hasRowInfo(graph) && myLID(graph.rowMap, localRow)
        getRowInfo(localRow).numEntries
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

function getGlobalRowCopy(graph::CRSGraph{GID}, globalRow::GID)::Array{GID, 1} where {GID <: Integer}
    Array{GID, 1}(getGlobalRowView(graph, globalRow))
end

function getGlobalRowView(graph::CRSGraph{GID}, globalRow::GID)::AbstractArray{GID, 1} where {GID <: Integer}
    debug = get(graph.plist, :debug, false)
    if isLocallyIndexed(graph)
        throw(InvalidArgumentError("The graph's indices are currently stored as local indices, so a view with global column indices cannot be returned.  Use getGlobalRowCopy(::CRSGraph) instead"))
    end

    if debug
        @assert hasRowInfo() "Graph row information was deleted"
    end
    rowInfo = getRowInfoFromGlobalRowIndex(globalRow)
    
    if rowInfo.localRow != 0 && rowInfo.numEntries > 0
        indices = view(getGlobalView, 1:rowInfo.numEntries)
        if debug
            @assert(length(indices) == getNumEntriesInGlobalRow(globalRow),
                "length(indices) = $(length(indices)) "
                * "!= getNumEntriesInGlobalRow(graph, $globalRow) "
                * "= $(getNumEntriesInGlobalRow(globalRow))")
        end
        indices
    else
        GID[]
    end
end
    


resumeFill(graph::CRSGraph; plist...) = resumeFill(graph, Dict(Array{Tuple{Symbol, Any}, 1}(plist)))

function resumeFill(graph::CRSGraph, plist::Dict{Symbol})
    if !hasRowInfo(graph)
        throw(InvalidStateException("Cannot resume fill of the CRSGraph, "
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
        throw(InvalidStateException("Graph fill state must be active to call fillComplete(...)"))
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
            throw(InvalidStateException("Only one process, but nonlocal entries are present"))
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
    debug = get(graph.plist, :debug, false)
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

#internal implementation of makeColMap, needed to handle some return and debuging stuff
#returns Tuple(errCode, colMap)
function __makeColMap(graph::CRSGraph{GID, PID, LID},domMap::BlockMap{GID, PID, LID},
        sortEachProcsGIDs::Bool) where {GID, PID, LID}
    errCode = 0#TODO improve from int error code

    if isnull(domMap)
        colMap = Nullable{BlockMap{GID, PID, LID}}()
    else
        myColumns = GID[]

        if isLocallyIndexed(graph)
            wrappedColMap = graph.colMap

            if isnull(wrappedColMap)
                warn("$(myPid(comm(graph))): The graph is locally indexed, but does not have a column map")

                errCode = -1
            else
                colMap = get(wrappedColMap)
                if linearMap(colMap) #i think isContiguous(map) <=> linearMap(map)?
                    numCurGIDs = numMyElements(colMap)
                    resize!(myColumns, numCurGIDs)
                    myFirstGlobalIndex = minMyGIDs(colMap)

                    for k = 1:numCurGIDs
                        myColumns[k] = myFirstGlobalIndex + k
                    end
                else
                    curGIDs = myGlobalElements(colMap)
                    numCurGIDs = length(curGIDs)
                    resize!(myColumns, numCurGIDs)
                    for k = 1:numCurGIDs
                        myColumns[k] = curGIDs[k]
                    end
                end
            end
        else #if graph.isGloballyIndexed
            numLocalColGIDs = 0
            numRemoteColGIDs = 0

            gidIsLocal = zeros(Bool, localNumRows)
            remoteGIDSet = Set()
            remoteGIDUnorderedVector = GID[]

            #if rowMap != null
            const rowMap = graph.rowMap

            for localRow = 1:localNumRows
                globalRow = gid(rowMap, localRow)
                rowGIDs = getGlobalRowView(graph)

                numEnt = length(rowGIDs)
                if numEnt != 0
                    for k = 1:numEnt
                        gid = rowGIDs[k]
                        lid = julia_petra.lid(domMap, gid)
                        if lid != 0
                            if !gidIsLocal[lid]
                                gidIsLocal = true
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
        else
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
                if debug
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
                    if debug
                        warn("$(myPid(comm(graph))): numLocalCount = $numLocalCount "
                            * "!= numLocalColGIDs = $numLocalColGIDs.  "
                            * "This should not happen.")
                    end
                    errCode = -4
                end
            end
        end

        #TODO look into FIXME on line 393
    end
    return(errCode, BlockMap(-1, -1, myColumns, comm(domMap)))
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
        throw(InvalidStateException("The graph must have a "
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
    if isOpt && get(graph.plist, :debug, false)
        @assert(getProfileType(graph) == STATIC_PROFILE,
            "Matrix claims optimized storage by profile type "
            * "is dynamic.  This shouldn't happend.")
    end
    isOpt
end

"""
    isFillActive(graph)::Bool

Whether the graph is in fill mode
"""
isFillActive(g::CRSGraph) = !g.fillComplete
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

#TODO implement DistObject methods
#TODO implement methods similar to RowMatrix
