
#=
k_numAllocPerRow_ and numAllocForAllRows_ are not copied to julia

/// FIXME (mfh 07 Aug 2014) We want graph's constructors to
/// allocate, rather than doing lazy allocation at first insert.
/// This will make both k_numAllocPerRow_ and numAllocForAllRows_
/// obsolete.
=#
#TODO remember to do allocations in constructor, not lazily

mutable struct CSRGraph{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::Nullable{BlockMap{GID, PID, LID}}
    rangeMap::Nullable{BlockMap{GID, PID, LID}}
    domainMap::Nullable{BlockMap{GID, PID, LID}}
    
    #may be null if domainMap and colMap are the same
    importer::Nullable{Import{GID, PID, LID}}
    #may be null if rangeMap and rowMap are the same
    exporter::Nullable{Export{GID, PID, LID}}

    #TODO is this needed? what type should it be?
    #local_graph_type lclGraph_;

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
function CSRGraph(rowMap::BlockMap{GID, PID, LID}, maxNumEntriesPerRow::LID,
        pftype::ProfileType, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    CSRGraph(rowMap, Nullable{BlockMap{GID, PID, LID}}(), maxNumEntriesPerRow, pftype, plist)
end
function CSRGraph(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        maxNumEntriesPerRow::LID, pftype::ProfileType, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    CSRGraph(rowMap, Nullable(colMap), maxNumEntriesPerRow, pftype, plist)
end
    
function CSRGraph(rowMap::BlockMap{GID, PID, LID}, colMap::Nullable{BlockMap{GID, PID, LID}}(),
        maxNumEntriesPerRow::LID, pftype::ProfileType, plist::Dict{Symbol}) where {
        GID <: Integer, PID <: Integer, LID <: Integer}
    graph = CSRGraph(
        rowMap,
        colMap,
        Nullable{BlockMap{GID, PID, LID}}(),
        Nullable{BlockMap{GID, PID, LID}}(),
        
        Nullable{Import{GID, PID, LID}}(),
        Nullable{Export{GID, PID, LID}}(),

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
    resumueFill(graph, params)
    checkInternalState(graph)
end



function map(graph::CSRGraph)
    graph.rowMap
end

#TODO implement Constructors
#TODO implement DistObject methods
#TODO implement methods similar to RowMatrix