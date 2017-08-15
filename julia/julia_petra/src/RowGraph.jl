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

    getRowMap(::RowGraph)
Gets the row map for the graph

    getColMap(::RowGraph)
Gets the column map for the graph

    getDomainMap(::RowGraph)
Gets the domain map for the graph

    getRangeMap(::RowGraph)
Gets the range map for the graph

    getImporter(::RowGraph)
Gets the graph's Import object

    getExporter(::RowGraph)
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

isFillActive(graph::RowGraph) = !isFillComplete(graph)

#### SrcDistObject methods ####
map(graph::RowGraph) = getRowMap(graph)
