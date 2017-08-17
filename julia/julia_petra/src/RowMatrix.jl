
export SrcDistRowMatrix, DistRowMatrix, RowMatrix
export isFillActive, isLocallyIndexed
export isFillComplete, getRowMap, hasColMap, getColMap, isGloballyIndexed, getGraph, getGlobalNumRows, getGlobalNumCols, getLocalNumRows, getLocalNumCols, getGlobalNumEntries, getLocalNumEntries, getNumEntriesInGlobalRow, getNumEntriesInLocalRow, getGlobalNumDiags, getLocalNumDiags, getGlobalMaxNumRowEntries, getLocalMaxNumRowEntries, isLowerTriangular, isUpperTriangular, getGlobalRowCopy, getLocalRowCopy, getGlobalRowView, getLocalRowView, getLocalDiagCopy, leftScale!, rightScale!

"""
The version of RowMatrix that isn't a subtype of DestObject
"""
abstract type SrcDistRowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: SrcDistObject{GID, PID, LID} 
end

"""
The version of RowMatrix that is a subtype of DestObject
"""
abstract type DistRowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{GID, PID, LID} 
end

#DECISION are any other mathmatical operations needed?

"""
RowMatrix is the base "type" for all row oriented matrices

RowMatrix is actually a type union of SrcDestRowMatrix and DestRowMatrix,
which are (direct) subtypes of SrcDestObject and DestObject, respectively.

All subtypes must have the following methods, with Impl standing in for the subtype:

    isFillComplete(mat::RowMatrix)
Whether `fillComplete(...)` has been called
    getRowMap(mat::RowMatrix)
Returns the BlockMap associated with the rows of this matrix
    hasColMap(mat::RowMatrix)
Whether the matrix has a column map
    getColMap(mat::RowMatrix)
Returns the BlockMap associated with the columns of this matrix
    isGloballyIndexed(mat::RowMatrix)
Whether the matrix stores indices with global indexes
    getGraph(mat::RowMatrix)
Returns the graph that represents the structure of the row matrix
    getGlobalNumRows(mat::RowMatrix)
Returns the number of rows across all processors
    getGlobalNumCols(mat::RowMatrix)
Returns the number of columns across all processors
    getLocalNumRows(mat::RowMatrix)
Returns the number of rows on the calling processor
    getLocalNumCols(mat::RowMatrix)
Returns the number of columns on the calling processor
    getGlobalNumEntries(mat::RowMatrix)
Returns the number of entries across all processors
    getLocalNumEntries(mat::RowMatrix)
Returns the number of entries on the calling processor
    getNumEntriesInGlobalRow(mat::RowMatrix, globalRow)
Returns the number of entries on the local processor in the given row
    getNumEntriesInLocalRow(mat::RowMatrix, localRow)
Returns the number of entries on the local processor in the given row
    getGlobalNumDiags(mat::RowMatrix)
Returns the number of diagonal elements across all processors
    getLocalNumDiags(mat::RowMatrix)
Returns the number of diagonal element on the calling processor
    getGlobalMaxNumRowEntries(mat::RowMatrix)
Returns the maximum number of row entries across all processors
    getLocalMaxNumRowEntries(mat::RowMatrix)
Returns the maximum number of row entries on the calling processor
    isLowerTriangular(mat::RowMatrix)
Whether the matrix is lower triangular
    isUpperTriangular(mat::RowMatrix)
Whether the matrix is upper triangular

    getGlobalRowCopy(matrix::RowMatrix{Data, GID, PID, LID}, globalRow::Integer)::Tuple{Array{GID, 1}, Array{Data, 1}}
Returns a copy of the given row using global indices
    getLocalRowCopy(matrix::RowMatrix{Data, GID, PID, LID},localRow::Integer)::Tuple{AbstractArray{LID, 1}, AbstractArray{Data, 1}}
Returns a copy of the given row using local indices
    getGlobalRowView(matrix::RowMatrix{Data, GID, PID, LID},globalRow::Integer)::Tuple{AbstractArray{GID, 1}, AbstractArray{Data, 1}} 
Returns a view to the given row using global indices
    getLocalRowView(matrix::RowMatrix{Data, GID, PID, LID},localRow::Integer)::Tuple{AbstractArray{GID, 1}, AbstractArray{Data, 1}}
Returns a view to the given row using local indices
    getLocalDiagCopy(matrix::RowMatrix{Data, GID, PID, LID})::MultiVector{Data, GID, PID, LID}
Returns a copy of the diagonal elements on the calling processor

    leftScale!(matrix::Impl{Data, GID, PID, LID}, X::Array{Data})
Scales matrix on the left with X

    rightScale!(matrix::Impl{Data, GID, PID, LID}, X::Array{Data})
Scales matrix on the right with X


Additionally, the following method must be implemented to fufil the operator interface:

    apply!(matrix::RowMatrix{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, mode::TransposeMode, alpha::Data, beta::Data)

    domainMap(operator::RowMatrix{Data, GID, PID, LID})::BlockMap{GID, PID, LID}

    rangeMap(operator::RowMatrix{Data, GID, PID, LID})::BlockMap{GID, PID, LID}

Finally, the required methods from DistObject must also be implemented.  `map(...)`, as required by SrcDistObject, is implemented to forward the call to `rowMap(...)`
"""
const RowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} = Union{SrcDistRowMatrix{Data, GID, PID, LID}, DistRowMatrix{Data, GID, PID, LID}}



function leftScale!(matrix::RowMatrix{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}) where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    if numVectors(X) != 1
        throw(InvalidArgumentError("Can only scale CRS matrix with column vector, not multi vector"))
    end
    leftScale!(matrix, X.data)
end

function rightScale!(matrix::RowMatrix{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}) where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    if numVectors(X) != 1
        throw(InvalidArgumentError("Can only scale CRS matrix with column vector, not multi vector"))
    end
    rightScale!(matrix, X.data)
end

isFillActive(matrix::RowMatrix) = !isFillComplete(matrix)
isLocallyIndexed(matrix::RowMatrix) = !isGloballyIndexed(matrix)

#for SrcDistObject
function map(matrix::RowMatrix)
    rowMap(matrix)
end


#TODO document
function getLocalDiagCopyWithoutOffsetsNotFillComplete(A::RowMatrix{Data, GID, PID, LID})::MultiVector{Data, GID, PID, LID} where {Data, GID, PID, LID}

    localRowMap = getLocalMap(getRowMap(A))
    localColMap = getLocalMap(getColMap(A))
    sorted = isSorted(A.myGraph)
    
    localNumRows = getLocalNumRows(A)
    diag = MultiVector{Data, GID, PID, LID}(getRowMap(A), 1)
    diagLocal1D = getVectorView(diag, 1)
    
    range = 1:localNumRows
    for localRowIndex in range
        diagLocal1D[localRowIndex] = 0
        globalIndex = gid(localRowMap, localRowIndex)
        localColIndex = lid(localColMap, globalIndex)
        if localColIndex != 0
            indices, values = getLocalRowView(A, localRowIndex)
            
            if !sorted
                offset = findfirst(indices, localColumnIndex)
            else
                offset = searchsorted(indices, localColumnIndex)
            end

            if offset <= length(indices)
                diagLocal1D[localRowIndex] = values[offset]
            end
        end
    end
    diag
end



#### required method documentation stubs ####
#TODO look into the fact that RowGraph and RowMatrix have methods that share function objects, thus causing documentation overwritting
"""
    isFillComplete(mat::RowMatrix)

Whether `fillComplete(...)` has been called
"""
function isFillComplete end

"""
    getRowMap(mat::RowMatrix)

Returns the BlockMap associated with the rows of this matrix
"""
function getRowMap end

"""
    hasColMap(mat::RowMatrix)

Whether the matrix has a column map
"""
function hasColMap end

"""
    getColMap(mat::RowMatrix)

Returns the BlockMap associated with the columns of this matrix
"""
function getColMap end

"""
    isGloballyIndexed(mat::RowMatrix)

Whether the matrix stores indices with global indexes
"""
function isGloballyIndexed end

"""
    getGraph(mat::RowMatrix)

Returns the graph that represents the structure of the row matrix
"""
function getGraph end

"""
    getGlobalNumRows(mat::RowMatrix)

Returns the number of rows across all processors
"""
function getGlobalNumRows end

"""
    getGlobalNumCols(mat::RowMatrix)

Returns the number of columns across all processors
"""
function getGlobalNumCols end

"""
    getLocalNumRows(mat::RowMatrix)

Returns the number of rows on the calling processor
"""
function getLocalNumRows end

"""
    getLocalNumCols(mat::RowMatrix)

Returns the number of columns on the calling processor
"""
function getLocalNumCols end

"""
    getGlobalNumEntries(mat::RowMatrix)

Returns the number of entries across all processors
"""
function getGlobalNumEntries end

"""
    getLocalNumEntries(mat::RowMatrix)

Returns the number of entries on the calling processor
"""
function getLocalNumEntries end

"""
    getNumEntriesInGlobalRow(mat::RowMatrix, globalRow)

Returns the number of entries on the local processor in the given row
"""
function getNumEntriesInGlobalRow end

"""
    getNumEntriesInLocalRow(mat::RowMatrix, localRow)

Returns the number of entries on the local processor in the given row
"""
function getNumEntriesInLocalRow end

"""
    getGlobalNumDiags(mat::RowMatrix)

Returns the number of diagonal elements across all processors
"""
function getGlobalNumDiags end

"""
    getLocalNumDiags(mat::RowMatrix)

Returns the number of diagonal element on the calling processor
"""
function getLocalNumDiags end

"""
    getGlobalMaxNumRowEntries(mat::RowMatrix)

Returns the maximum number of row entries across all processors
"""
function getGlobalMaxNumRowEntries end

"""
    getLocalMaxNumRowEntries(mat::RowMatrix)

Returns the maximum number of row entries on the calling processor
"""
function getLocalMaxNumRowEntries end

"""
    isLowerTriangular(mat::RowMatrix)

Whether the matrix is lower triangular
"""
function isLowerTriangular end

"""
    isUpperTriangular(mat::RowMatrix)

Whether the matrix is upper triangular
"""
function isUpperTriangular end

"""
    getGlobalRowCopy(matrix::RowMatrix{Data, GID, PID, LID}, globalRow::Integer)::Tuple{Array{GID, 1}, Array{Data, 1}}

Returns a copy of the given row using global indices
"""
function getGlobalRowCopy end

"""
    getLocalRowCopy(matrix::RowMatrix{Data, GID, PID, LID},localRow::Integer)::Tuple{AbstractArray{LID, 1}, AbstractArray{Data, 1}}

Returns a copy of the given row using local indices
"""
function getLocalRowCopy end

"""
    getGlobalRowView(matrix::RowMatrix{Data, GID, PID, LID},globalRow::Integer)::Tuple{AbstractArray{GID, 1}, AbstractArray{Data, 1}} 

Returns a view to the given row using global indices
"""
function getGlobalRowView end

"""
    getLocalRowView(matrix::RowMatrix{Data, GID, PID, LID},localRow::Integer)::Tuple{AbstractArray{GID, 1}, AbstractArray{Data, 1}}

Returns a view to the given row using local indices
"""
function getLocalRowView end

"""
    getLocalDiagCopy(matrix::RowMatrix{Data, GID, PID, LID})::MultiVector{Data, GID, PID, LID}

Returns a copy of the diagonal elements on the calling processor
"""
function getLocalDiagCopy end

"""
    leftScale!(matrix::Impl{Data, GID, PID, LID}, X::Array{Data})

Scales matrix on the left with X
"""
function leftScale! end

"""
    rightScale!(matrix::Impl{Data, GID, PID, LID}, X::Array{Data})

Scales matrix on the right with X
"""
function rightScale! end
