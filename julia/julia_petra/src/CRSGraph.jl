
export CRSGraph

#=
k_numAllocPerRow_ and numAllocForAllRows_ are not copied to julia

/// FIXME (mfh 07 Aug 2014) We want graph's constructors to
/// allocate, rather than doing lazy allocation at first insert.
/// This will make both k_numAllocPerRow_ and numAllocForAllRows_
/// obsolete.
=#
#TODO remember to do allocations in constructor, not lazily

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
    localIndices2D::Array{LID, 2}
    globalIndices2D::Array{LID, 2}
    #may exist in 1-D storage if not packed
    numRowEntries::Array{LID, 1}

    storageStatus::StorageStatus

    indiciesAllowed::Bool
    indiciesType::IndexType
    fillComplete::Bool

    lowerTriangle::Bool
    upperTriangle::Bool
    indiciesAreSorted::Bool
    noRedundancies::Bool
    haveLocalConstants::Bool
    haveGlobalConstants::Bool

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

        indiciesType::IndexType
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
            Array{LID, 2}(0, 0),
            Array{LID, 2}(0, 0),
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
        
        UNKNOWN
    )
        
    #TODO do allocations
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
    graph = CRSGraph(
        rowMap,
        colMap,
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),

        LocalCRSGraph{LID, LID}(), #lclGraph
        
        LID(0), #nodeNumEntries

        #Whether the graph was allocated with static or dynamic profile.
        pftype,
        
        #TODO find numAllocForAllRows (see line 248)

        (pftype == STATIC_PROFILE ?
              STORAGE_1D_UNPACKED 
            : STORAGE_2D),

        UNKNOWN
    )
    
    #DECISION allow rowMap to be null?
    lclNumRows = numLocalElements(rowMap)
    if length(numEntPerRow) != lclNumRows
        throw(InvalidArgumentError("numEntPerRows has length $(length(numEntPerRow)) " *
                "!= the local number of rows $lclNumRows as spcified by the input row Map"))
    end
    
    if get(plist, "debug", false)
        for r = 1:lclNumRows
            curRowCount = numEntPerRow[r]
            if curRowCount <= 0
                throw(InvalidArgumentError("numEntPerRow[$r] = $curRowCount is not valid"))
            end
        end
    end
    #TODO do allocations
    
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
        
        #TODO figure out numAllocForAllRows

        STORAGE_1D_PACKED,

        LOCAL
    )
    #TODO do allocations
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
        
        #TODO figure out numAllocForAllRows

        STORAGE_1D_PACKED,
        
        LOCAL
    )
    
    if numRows(localGraph) != numMyElements(rowMap)
        throw(InvalidArgumentError("input row map and input local "
                * "graph need to have the same number of rows.  The "
                * "row map claims $(numMyElements(rowMap)) row(s), "
                * "but the local graph claims $(numRows(localGraph)) "
                * "row(s)."))
    end
    
    makeImportExport(graph)
    
    d_inds = localGraph.entries
    graph.localIndices1D = d_inds
    
    d_ptrs = localGraph.rowMap
    graph.rowOffsets = d_ptrs
    
    
    #TODO figure out if these can be computed pre-inner constructor call
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
    
    
    

#TODO implement computeGlobalConstants(::CRSGraph)
#TODO implement makeImportExport(::CRSGraph)
#TODO implement resumeFill(::CRSGraph, ::Dict{Symbol})
function resumeFill(g::CRSGraph, d::Dict{Symbol}) end
#TODO implement checkInternalState(::CRSGraph)
function checkInternalState(g::CRSGraph) end
#TODO implement setAllIndices(::CRSGraph, ::Array{LID, 1}, ::Array{LID, 1}) 

function map(graph::CRSGraph)
    graph.rowMap
end

#TODO implement DistObject methods
#TODO implement methods similar to RowMatrix