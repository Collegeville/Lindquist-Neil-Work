
export CRSGraph, isLocallyIndexed, isGloballyIndexed, getProfileType


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
    indiciesAreSorted::Bool
    noRedundancies::Bool
    haveLocalConstants::Bool
    haveGlobalConstants::Bool

    debug::Bool

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
        debug::Bool
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

            debug,
        
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
    debug = get(plist, "debug", false)
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
        
        LOCAL,
        debug
    )
        
    allocateIndices(graph, LOCAL, maxNumEntriesPerRow)
    
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
    debug = get(plist, "debug", false)
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

        LOCAL,
        debug
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
    
    allocateIndices(graph, LOCAL, numEntPerRow)
    
    resumeFill(graph, plist)
    checkInternalState(graph)
    
    graph
end


function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        rowPointers::Array{LID, 1}, columnIndices::Array{LID, 1},
        plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    debug = get(plist, "debug", false)
    graph = CRSGraph(
        rowMap,
        Nullable(colMap),
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),

        LocalCRSGraph{LID, LID}(), #lclGraph
        
        LID(0), #nodeNumEntries
        
        STATIC_PROFILE,
        
        STORAGE_1D_PACKED,

        LOCAL,
        debug
    )
    #seems to be already taken care of
    #allocateIndices(graph, LOCAL, numEntPerRow)
    
    setAllIndicies(graph, rowPointers, columnIndicies)
    checkInternalState(graph)
    
    graph
end


function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localGraph::LocalCRSGraph{LID, LID}, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    debug = get(plist, "debug", false)
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
        
        LOCAL,
        debug
    )
    
    if numRows(localGraph) != numMyElements(rowMap)
        throw(InvalidArgumentError("input row map and input local "
                * "graph need to have the same number of rows.  The "
                * "row map claims $(numMyElements(rowMap)) row(s), "
                * "but the local graph claims $(numRows(localGraph)) "
                * "row(s)."))
    end
    
    #seems to be already taken care of
    #allocateIndices(graph, LOCAL, numEntPerRow)
    
    makeImportExport(graph)
    
    d_inds = localGraph.entries
    graph.localIndices1D = d_inds
    
    d_ptrs = localGraph.rowMap
    graph.rowOffsets = d_ptrs
    
    #reset local properties
    graph.upperTriangular = true
    graph.lowerTriangular = true
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
                    graph.upperTriangular = false
                end
                if localRow < largestCol
                    graph.lowerTriangular = false
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
    
    if graph.debug
        @assert !null(graph.colMap) "The graph must have a column map at this point"
    end
    
    computeLocalConstants()
    
    #if graph.haveGlobalConstants == false  #short circuited above
    graph.globalNumEntries, graph.globalNumDiags = sumAll(comm(graph.map),
        [GID(graph.nodeNumEntries), GID(graph.nodeNumDiags)])
        
    graph.globalMaxNumRowEntries = maxAll(comm(graph.map), GID(graph.nodeMaxNumRowEntries))
    graph.haveGlobalConstants = true
end


function computeLocalConstants(graph::CRSGraph{GID, PID, LID}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    
    #short circuit if already computed
    graph.haveLocalConstants && return
    
    if graph.debug
        @assert !null(graph.colMap) "The graph must have a column map at this point"
    end
    
    #if graph.haveLocalConstants == false  #short circuited above
    
    graph.upperTriangular = true
    graph.lowerTriangular = true
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
                graph.upperTriangular = false
            end
            if localRow < largestCol
                graph.lowerTriangular = false
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

function getRowInfo(graph::CRSGraph{GID, PID, LID}, row::LID)::RowInfo{LID} where {GID, PID, LID <: Integer}
    if graph.debug
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
        lg::IndexType, numAllocPerRow::Union{Integer, Array{<:Integer, 1}}) where {
        GID <: Integer, LID <: Integer}
    @assert isLocallyIndexed(graph) == (lg == LOCAL) "Graph is locally indexed, by lg=$lg"
    @assert isGloballyIndexed(graph) == (lg == GLOBAL) "Graph is globally indexed by lg=$lg"
    
    numRows = getNodeNumRows(graph)

    if getProfileType(graph) == STATIC_PROFILE
        rowPtrs = Array{LID, 1}(numRows + 1)

        #using function call style, since otherwise it needs to be 1 line or julia funcalls a Bool
        @assert(!isa(numAllocPerRow, Array{<:Integer}) 
            || length(numAllocPerRow) == numRows,
            "numAllocRows has length = $(length(numAllocPerRow)) "
            * "!= numRows = $numRows")
        computeOffsets(rowPtrs, numAllocPerRow)
        
        graph.rowOffsets = rowPtrs
        numInds = rowPtrs[numRows+1]
        
        if lg == LOCAL
            graph.localIndices1D = Array{LID, 1}(numInds)
        else
            graph.globalIndices1D = Array{GID, 1}(numInds)
        end
        graph.storageStatus = STORAGE_1D_UNPACKED
    else
        numAllocIsArray = isa(numAllocPerRow, Array{<: Integer})
        if lg == LOCAL
            graph.localInds2D = Array{Array{LID, 1}, 1}(numRows)
            for i = 1:numRows
                const howMany = (numAllocIsArray ?
                    numAllocPerRow[i] : numAllocPerRow)
                if howMany > 0
                    resize(graph.localInds2D[i], howMany)
                end
            end
        else #lg == GLOBAL
            graph.globalInds2D = Array{Array{GID, 1}, 1}(numRows)
            for i = 1:numRows
                const howMany = (numAllocIsArray ?
                    numAllocPerRow[i] : numAllocPerRow)
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
    
    checkInternalState(graph)
end
    
    
#TODO implement makeImportExport(::CRSGraph)
#TODO implement resumeFill(::CRSGraph, ::Dict{Symbol})
function resumeFill(g::CRSGraph, d::Dict{Symbol}) end
#TODO implement checkInternalState(::CRSGraph)
function checkInternalState(g::CRSGraph) end
#TODO implement setAllIndices(::CRSGraph, ::Array{LID, 1}, ::Array{LID, 1}) 


#### external API ####

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
    graph.indicesType == GLOBAL
end
   
"""
    isLocallyIndexed(::CRSGraph)

Whether the graph uses local indexes
"""
function isLocallyIndexed(graph::CRSGraph)
    graph.indicesType == LOCAL
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
