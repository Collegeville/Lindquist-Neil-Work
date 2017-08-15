export SrcDistRowGraph, DistRowGraph, RowGraph
#required methods
export getRowMap, getColMap, getDomainMap, getRangeMap, getImporter, getExporter
export getGlobalNumRows, getGlobalNumCols, getGlobalNumEntries, getGlobalNumDiags
export getNodeNumRows, getNodeNumCols, getNodeNumEntries, getNodeNumDiags
export getNumEntriesInGlobalRow, getNumEntriesInLocalRow
export getGlobalMaxNumRowEntries, getNodeMaxNumRowEntries
export hasColMap, isLowerTriangular, sUpperTriangular
export isLocallyIndexed, isGloballyIndexed, isFillComplete
export getGlobalRowCopy, getLocalRowCopy, pack
#implemented methods
export isFillActive

"""
The version of RowMatrix that isn't a subtype of DistObject
"""
abstract type SrcDistRowGraph{GID <: Integer, PID <: Integer, LID <: Integer} <: SrcDistObject{GID, PID, LID} 
end

"""
The version of RowMatrix that is a subtype of DistObject
"""
abstract type DistRowGraph{GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{GID, PID, LID} 
end

#TODO change function names that reference "Node" to reference "Local"

"""
RowGraph is the base "type" for all row oriented storage graphs

RowGraph is actually a type union of SrcDistRowGraph and DistRowGraph,
which are (direct) subtypes of SrcDistObject and DistObject, respectively.

Instances of these types are required to implement the following submethods

    getRowMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}
Gets the row map for the graph

    getColMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}
Gets the column map for the graph

    getDomainMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}
Gets the domain map for the graph

    getRangeMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}
Gets the range map for the graph

    getImporter(::RowGraph{GID, PID, LID})::Import{GID, PID, LID}
Gets the graph's Import object

    getExporter(::RowGraph{GID, PID, LID})::Export{GID, PID, LID}
Gets the graph's Export object

    getGlobalNumRows(::RowGraph{GID})::GID
Returns the number of global rows in the graph

    getGlobalNumCols(::RowGraph{GID})::GID
Returns the number of global columns in the graph

    getNodeNumRows(::RowGraph{GID, PID, LID})::LID
Returns the number of rows owned by the calling process

    getNodeNumCols(::RowGraph{GID, PID, LID})::LID
Returns the number of columns owned by teh calling process

    getGlobalNumEntries(::RowGraph{GID, PID, LID})::GID
Returns the global number of entries in the graph

    getNodeNumEntries(::RowGraph{GID, PID, LID})::LID
Returns the local number of entries in the graph

    getNumEntriesInGlobalRow(::RowGraph{GID, PID, LID}, row::GID)::LID
Returns the current number of local entries in the given row

    getNumEntriesInLocalRow(::RowGraph{GID, PID, LID}, row::LID)::LID
Returns the current number of local entries in the given row

    getGlobalNumDiags(::RowGraph{GID, PID, LID})::GID
Returns the global number of diagonal entries

    getNodeNumDiags(::RowGraph{GID, PID, LID})::LID
Returns the local number of diagonal entries

    getGlobalMaxNumRowEntries(::RowGraph{GID, PID, LID})::LID
Returns the maximum number of entries across all rows/columns on all processors

    getNodeMaxNumRowEntries(::RowGraph{GID, PID, LID})::LID
Returns the maximum number of entries across all rows/columns on this processor

    hasColMap(::RowGraph{GID, PID, LID})::Bool
Whether the graph has a well-defined column map

    isLowerTriangular(::RowGraph{GID, PID, LID})::Bool
Whether the graph is lower trianguluar

    isUpperTriangular(::RowGraph{GID, PID, LID})::Bool
Whether the graph is upper trianguluar

    isLocallyIndexed(::RowGraph)::Bool
Whether the graph is using local indices

    isGloballyIndexed(::RowGraph)::Bool
Whether the graph is using global indices

    isFillComplete(::RowGraph)
Whether `fillComplete()` has been called

    getGlobalRowCopy(::RowGraph{GID, PID, LID}, row::GID)::Array{GID, 1}
Extracts a copy of the given row of the graph

    getLocalRowCopy(::RowGraph{GID, PID, LID}, row::LID)::Array{LID, 1}
Extracts a copy of the given row of the graph

    pack(::RowGraph{GID, PID, LID}, exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID})::Array{Array{LID, 1}}
Packs this object's data for import or export
"""
const RowGraph{GID <: Integer, PID <: Integer, LID <: Integer} = Union{SrcDistRowGraph{GID, PID, LID}, DistRowGraph{GID, PID, LID}}

"""
    isFillActive(::RowGraph)

Whether the graph is being built
"""
isFillActive(graph::RowGraph) = !isFillComplete(graph)


"""
    getNumEntriesInGlobalRow(graph::RowGraph{GID, PID, LID}, row::Integer)::LID

Returns the current number of local entries in the given row
"""
function getNumEntriesInGlobalRow(graph::RowGraph{GID, PID, LID},
        row::Integer)::LID where{GID, PID, LID}
    getNumEntriesInGlobalRow(graph, GID(row))
end

"""
    getNumEntriesInLocalRow(::RowGraph{GID, PID, LID}, row::Integer)::LID

Returns the current number of local entries in the given row
"""
function getNumEntriesInLocalRow(graph::RowGraph{GID, PID, LID},
        row::Integer)::LID where{GID, PID, LID}
    getNumEntriesInLocalRow(graph, LID(row))
end

"""
    getGlobalRowCopy(::RowGraph{GID, PID, LID}, row::Integer)::Array{GID, 1}

Extracts a copy of the given row of the graph
"""
function getGlobalRowCopy(graph::RowGraph{GID, PID, LID},
        row::Integer)::Array{GID, 1} where{GID, PID, LID}
    getGlobalRowCopy(graph, GID(row))
end

"""
    getLocalRowCopy(::RowGraph{GID, PID, LID}, row::Integer)::Array{LID, 1}

Extracts a copy of the given row of the graph
"""
function getLocalRowCopy(graph::RowGraph{GID, PID, LID},
        row::Integer)::Array{LID, 1} where{GID, PID, LID}
    getLocalRowCopy(graph, LID(row))
end


#### SrcDistObject methods ####
map(graph::RowGraph) = getRowMap(graph)


#### documentation for required methods ####
"""
    getRowMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}

Gets the row map for the graph
"""
function getRowMap end

"""
    getColMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}

Gets the column map for the graph
"""
function getColMap end

"""
    getDomainMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}

Gets the domain map for the graph
"""
function getDomainMap end

"""
    getRangeMap(::RowGraph{GID, PID, LID})::BlockMap{GID, PID, LID}

Gets the range map for the graph
"""
function getRangeMap end

"""
    getImporter(::RowGraph{GID, PID, LID})::Import{GID, PID, LID}

Gets the graph's Import object
"""
function getImporter end

"""
    getExporter(::RowGraph{GID, PID, LID})::Export{GID, PID, LID}

Gets the graph's Export object
"""
function getExporter end

"""
    getGlobalNumRows(::RowGraph{GID})::GID

Returns the number of global rows in the graph
"""
function getGlobalNumRows end

"""
    getGlobalNumCols(::RowGraph{GID})::GID

Returns the number of global columns in the graph
"""
function getGlobalNumCols end

"""
    getNodeNumRows(::RowGraph{GID, PID, LID})::LID

Returns the number of rows owned by the calling process
"""
function getNodeNumRows end

"""
    getNodeNumCols(::RowGraph{GID, PID, LID})::LID

Returns the number of columns owned by teh calling process
"""
function getNodeNumCols end

"""
    getGlobalNumEntries(::RowGraph{GID, PID, LID})::GID

Returns the global number of entries in the graph
"""
function getGlobalNumEntries end

"""
    getNodeNumEntries(::RowGraph{GID, PID, LID})::LID

Returns the local number of entries in the graph
"""
function getNodeNumEntries end

"""
    getGlobalNumDiags(::RowGraph{GID, PID, LID})::GID

Returns the global number of diagonal entries
"""
function getGlobalNumDiags end

"""
    getNodeNumDiags(::RowGraph{GID, PID, LID})::LID

Returns the local number of diagonal entries
"""
function getNodeNumDiags end

"""
    getGlobalMaxNumRowEntries(::RowGraph{GID, PID, LID})::LID

Returns the maximum number of entries across all rows/columns on all processors
"""
function getGlobalMaxNumRowEntries end

"""
    getNodeMaxNumRowEntries(::RowGraph{GID, PID, LID})::LID

Returns the maximum number of entries across all rows/columns on this processor
"""
function getNodeMaxNumRowEntries end

"""
    hasColMap(::RowGraph{GID, PID, LID})::Bool

Whether the graph has a well-defined column map
"""
function hasColMap end

"""
    isLowerTriangular(::RowGraph{GID, PID, LID})::Bool

Whether the graph is lower trianguluar
"""
function isLowerTriangular end

"""
    isUpperTriangular(::RowGraph{GID, PID, LID})::Bool

Whether the graph is upper trianguluar
"""
function isUpperTriangular end

"""
    isLocallyIndexed(::RowGraph)::Bool

Whether the graph is using local indices
"""
function isLocallyIndexed end

"""
    isGloballyIndexed(::RowGraph)::Bool

Whether the graph is using global indices
"""
function isGloballyIndexed end

"""
    isFillComplete(::RowGraph)

Whether `fillComplete()` has been called
"""
function isFillComplete end

"""
    pack(::RowGraph{GID, PID, LID}, exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID})::Array{Array{LID, 1}}

Packs this object's data for import or export
"""
function pack end
