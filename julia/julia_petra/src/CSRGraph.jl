
#=
k_numAllocPerRow_ and numAllocForAllRows_ are not copied to julia

/// FIXME (mfh 07 Aug 2014) We want graph's constructors to
/// allocate, rather than doing lazy allocation at first insert.
/// This will make both k_numAllocPerRow_ and numAllocForAllRows_
/// obsolete.
=#
#TODO remember to do allocations in constructor, not lazily

mutable struct CSRGraph{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{Data, GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::BlockMap{GID, PID, LID}
    rangeMap::BlockMap{GID, PID, LID}
    domainMap::BlockMap{GID, PID, LID}
    
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


function map(graph::CGSGraph)
    graph.rowMap
end

#TODO implement DistObject methods
#TODO implement methods similar to RowMatrix