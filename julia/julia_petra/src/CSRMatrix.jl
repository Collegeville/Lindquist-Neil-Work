export CSRMatrix

mutable struct CSRMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistRowMatrix{Data, GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::Nullable{BlockMap{GID, PID, LID}}
    
    #TODO look into how nessacery these are
    importMV::Nullable{MultiVector{Data, GID, PID, LID}}
    exportMV::Nullable{MultiVector{Data, GID, PID, LID}}

    myGraph::CRSGraph{GID, PID, LID}
    
    #TODO figure out
    localMatrix::LocalCSRMatrix{Data, LID}

    values2D::Array{Array{Data, 1}, 1}

    #pull storageStatus and fillComplete  from graph

    #Dict keys are row indices
    #first element of each tuple is a column index
    #second element of each tuple is the matching entry
    nonlocals::Dict{GID, Array{Tuple{GID, Data}, 1}}

    function CSRMatrix{Data, GID, PID, LID}(rowMap::BlockMap{GID, PID, LID}, colMap::Nullable{BlockMap{GID, PID, LID}}, myGraph::CRSGraph{GID, PID, LID}, localMatrix::LocalCSRMatrix{Data, LID}, values2D::Array{Array{Data, 1}, 1}) where {Data, GID, PID, LID}
        new(rowMap,
            colMap,
            Nullable{MultiVector{Data, GID, PID, LID}}(),
            Nullable{MultiVector{Data, GID, PID, LID}}(),
            myGraph,
            values2D,
            Dict{GID, Array{Tuple{GID, Data}, 1}}())
    end
        
end

#### Constructors ####
#TODO implement Constructors
#TODO document Constructors
#TODO make sure to add constructor to convert from Julia's data types

function CSRMatrix{Data}(rowMap::BlockMap{GID, PID, LID},
        maxNumEntriesPerRow::Union{Integer, Array{Integer, 1}}, 
        pftype::ProfileType, plist::Dict{Symbol}) where {Data, GID, PID, LID}
    CSRMatrix{Data}(rowMap, Nullable{BlockMap{GID, PID, LID}}(),
        maxNumEntriesPerRow, pftype, plist)
end
function CSRMatrix{Data}(rowMap::BlockMap{GID, PID, LID},
        colMap::BlockMap{GID, PID, LID},
        maxNumEntriesPerRow::Union{Integer, Array{Integer, 1}}, 
        pftype::ProfileType, plist::Dict{Symbol}) where {Data, GID, PID, LID}
    CSRMatrix{Data}(rowMap, Nullable(colMap), maxNumEntriesPerRow,
        pftype, plist)
end
function CSRMatrix{Data}(rowMap::BlockMap{GID, PID, LID},
        colMap::Nullable{BlockMap{GID, PID, LID}},
        maxNumEntriesPerRow::Union{Integer, Array{Integer, 1}}, 
        pftype::ProfileType, plist::Dict{Symbol}) where {Data, GID, PID, LID}
    graph = CRSGraph(rowMap, maxNumEntriesPerRow, pftype, plist)
    
    matrix = CSRMatrix{Data, GID, PID, LID}(rowMap, colMap, 
        graph, LocalCSRMatrix{Data, LID}(), [])
    
    resumeFill(matrix, plist)
    checkInternalState(matrix)
    
    matrix
end

function CSRMatrix{Data}(graph::CRSGraph{GID, PID, LID},plist::Dict{Symbol}
        ) where {Data, GID, PID, LID}
    numCols = numMyElements(getColMap(graph))
    localGraph = getLocalGraph(graph)
    val = Array{Data, 1}(length(localGraph.entries))
    localMatrix = LocalCSRMatrix(numCols, val, localGraph)
    
    matrix = CSRMatrix(graph.rowMap, graph.colMap, graph, localMatrix, [])
    
    checkInternalState(matrix)
    
    matrix
end

function CSRMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        rowOffsets::Array{LID, 1}, colIndices::Array{LID, 1}, values::Array{Data, 1},
        plist::Dict{Symbol}) where {Data, GID, PID, LID}
    
    #check user's input.  Might throw on only some processes, causing deadlock
    if length(values) != length(colIndices)
        throw(InvalidArgumentError("values and columnIndices must "
                * "have the same length"))
    end
    #TODO add debug from line 328?
    #if so, use a minAll on the previous error when debuging to propagate the error
    
    graph = CRSGraph(rowMap, colMap, rowOffsets, columnIndices, plist)
    localGraph = getLocalGraph(graph)
    
    numCols = numMyElements(getColMap(graph))
    localMatrix = LocalCSRMatrix(numCols, values, localGraph)
    
    matrix = CSRMatrix(rowMap, Nullable(colMap), graph, localMatrix, [])
    
    checkInternalState(matrix)
    matrix
end


function CSRMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localMatrix::LocalCSRMatrix{Data, LID}, plist::Dict{Symbol}
        ) where {Data, GID, PID, LID}
    
    graph = CRSGraph(rowMap, colMap, localMatrix.graph, plist)
    
    matrix = CSRMatrix(rowMap, colMap, graph, localMatrix, [])
    
    computeGlobalConstants(matrix)
    checkInternalState(matrix)
    matrix
end

function CSRMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localMatrix::AbstractArray{Data, 2}, plist::Dict{Symbol}
        ) where {Data, GID, PID, LID}
    linearIndices = find(x -> x!=0, localMatrix)
    rowIndicesIter, colIndicesIter, valuesIter = zip(
        sort!(collect(zip(ind2sub(size(localMatrix), linearIndices)...,
                          localMatrix[linearIndices])))...)
    rowIndices = collect(rowIndicesIter)
    rowOffsets = Array{LID, 1}(size(localMatrix, 1)+1)
    row = 1
    j = 1
    for i in 1:length(rowIndices)
        if rowIndices[i] > row
            row += 1
            rowOffsets[row] = i
        end
    end
    rowOffsets[length(rowOffsets)] = length(rowIndices)+1
    
    CSRMatrix(rowMap, colMap, rowOffsets, 
        collect(colIndicesIter), collect(valuesIter), plist)
end


#### Internal methods ####
#TODO implement checkInternalState(matrix)

#### External methods ####
#TODO document external methods

#TODO implement resumeFill(matrix, plist)

#TODO implement RowMatrix methods
#TODO implement DistObject methods
#TODO implement SrcDistObject methods