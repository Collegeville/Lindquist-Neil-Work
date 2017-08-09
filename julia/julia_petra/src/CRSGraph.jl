
#=
k_numAllocPerRow_ and numAllocForAllRows_ are not copied to julia

/// FIXME (mfh 07 Aug 2014) We want graph's constructors to
/// allocate, rather than doing lazy allocation at first insert.
/// This will make both k_numAllocPerRow_ and numAllocForAllRows_
/// obsolete.
=#
#TODO remember to do allocations in constructor, not lazily

#TODO figure out type of Data
mutable struct CRSGraph{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::Nullable{BlockMap{GID, PID, LID}}
    rangeMap::Nullable{BlockMap{GID, PID, LID}}
    domainMap::Nullable{BlockMap{GID, PID, LID}}
    
    #may be null if domainMap and colMap are the same
    importer::Nullable{Import{GID, PID, LID}}
    #may be null if rangeMap and rowMap are the same
    exporter::Nullable{Export{GID, PID, LID}}

    lclGraph::Array{Data, 2}

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
end


#### Constructors #####

#TODO add plist call that passes kwargs
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
        
        Nullable{Import{GID, PID, LID}}(),
        Nullable{Export{GID, PID, LID}}(),

        [], #lclGraph
        
        0, #nodeNumEntries
        #using -1 to indicate uninitiallized, likely to cause an error if used
        -1, #nodeNumDiags
        -1, #nodeMaxNumRowEntries
        -1, #globalNumEntries
        -1, #globalNumDiags
        -1, #globalMaxNumRowEntries

        pftype,

        ## 1-D storage (Static profile) data structures ##
        [],
        [],
        [],

        ## 2-D storage (Dynamic profile) data structures ##
        [],
        [],
        [],

        (pftype == STATIC_PROFILE ?
              STORAGE_1D_UNPACKED 
            : STORAGE_2D),
        
        false,
        UNKNOWN,
        false,

        false,
        false,
        true,
        true,
        false,
        false,

        Dict{GID, Array{GID, 1}}()
    )
        
        
    staticAssertions(graph)
    #TODO do allocations
    resumueFill(graph, params)
    checkInternalState(graph)
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

        Nullable{Import{GID, PID, LID}}(),
        Nullable{Export{GID, PID, LID}}(),

        [], #lclGraph
        
        0, #nodeNumEntries
        #using -1 to indicate uninitiallized, likely to cause an error if used
        -1, #nodeNumDiags
        -1, #nodeMaxNumRowEntries
        -1, #globalNumEntries
        -1, #globalNumDiags
        -1, #globalMaxNumRowEntries

        #Whether the graph was allocated with static or dynamic profile.
        pftype,
        
        #TODO find numAllocForAllRows (see line 248)

        ## 1-D storage (Static profile) data structures ##
        [],
        [],
        [],

        ## 2-D storage (Dynamic profile) data structures ##
        [],
        [],
        [],

        (pftype == STATIC_PROFILE ?
              STORAGE_1D_UNPACKED 
            : STORAGE_2D),

        false,
        UNKNOWN,
        false,

        false,
        false,
        true,
        true,
        false,
        false,

        Dict{GID, Array{GID, 1}}()
    )
    
    staticAssertions(graph);
    
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
end


function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        rowPointers::Array{LID, 1}, columnIndices::Array{LID, 1},
        plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    graph = CRSGraph(
        rowMap,
        Nullable(colMap),
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),
        
        Nullable{Import{GID, PID, LID}}(),
        Nullable{Export{GID, PID, LID}}(),

        [], #lclGraph
        
        0, #nodeNumEntries
        #using -1 to indicate uninitiallized, likely to cause an error if used
        -1, #nodeNumDiags
        -1, #nodeMaxNumRowEntries
        -1, #globalNumEntries
        -1, #globalNumDiags
        -1, #globalMaxNumRowEntries
        
        STATIC_PROFILE,
        
        #TODO figure out numAllocForAllRows
        
        ## 1-D storage (Static profile) data structures ##
        [],
        [],
        [],

        ## 2-D storage (Dynamic profile) data structures ##
        [],
        [],
        [],

        STORAGE_1D_PACKED,

        false,
        LOCAL,
        false,

        false,
        false,
        true,
        true,
        false,
        false,

        Dict{GID, Array{GID, 1}}()
    )
    
    staticAssertions(graph)
    #TODO do allocations
    setAllIndicies(graph, rowPointers, columnIndicies)
    checkInternalState(graph)
end

#This method appears to require Kokkos.StaticCRSGraph
#=
function CRSGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localGraph::Kokkos.StaticCRSGraph, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    graph = CRSGraph(
        rowMap,
        Nullable(colMap),
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),
        
        Nullable{Import{GID, PID, LID}}(),
        Nullable{Export{GID, PID, LID}}(),

        localGraph,
        
        length(localGraph), #nodeNumEntries
        #using -1 to indicate uninitiallized, likely to cause an error if used
        -1, #nodeNumDiags
        -1, #nodeMaxNumRowEntries
        -1, #globalNumEntries
        -1, #globalNumDiags
        -1, #globalMaxNumRowEntries
        
        STATIC_PROFILE,
        
        #TODO figure out numAllocForAllRows
        
        ## 1-D storage (Static profile) data structures ##
        [],
        [],
        [],

        ## 2-D storage (Dynamic profile) data structures ##
        [],
        [],
        [],

        STORAGE_1D_PACKED,
        
        false,
        LOCAL,
        false,

        false,
        false,
        true,
        true,
        false,
        false,

        Dict{GID, Array{GID, 1}}()
    )
    =#
    
#TODO group duplicate constructor code and staticAssertions into inner constructor    
    

#TODO implement staticAssertions(::CRSGraph)
#TODO implement resumeFill(::CRSGraph, ::Dict{Symbol})
#TODO implement checkInternalState(::CRSGraph)
#TODO implement setAllIndices(::CRSGraph, ::Array{LID, 1}, ::Array{LID, 1}) 

function map(graph::CRSGraph)
    graph.rowMap
end

#TODO implement Constructors
#TODO implement DistObject methods
#TODO implement methods similar to RowMatrix