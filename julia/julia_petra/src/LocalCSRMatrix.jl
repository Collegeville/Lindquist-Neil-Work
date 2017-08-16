export LocalCSRMatrix, numRows, numCols, getRowView

struct LocalCSRMatrix{Data, IndexType <: Integer}
    graph::LocalCRSGraph{IndexType, IndexType}
    values::Array{Data, 1}

    numCols::IndexType
end

function LocalCSRMatrix{Data, IndexType}() where {Data, IndexType}
    LocalCSRMatrix(LocalCRSGraph{IndexType, IndexType}(), Data[], IndexType(0))
end

function LocalCSRMatrix(nRows::Integer, nCols::Integer,
        vals::Array{Data, 1}, rows::Array{IndexType, 1},
        cols::Array{IndexType, 1}) where {Data, IndexType}
    if length(rows) != nRows + 1
        throw(InvalidArgumentError("length(rows) = $(length(rows)) != nRows+1 "
                * "= $(nRows + 1)"))
    end
    LocalCSRMatrix(LocalCRSGraph(cols, rows), vals, IndexType(nCols))
end

"""
    nnz(::LocalCRSMatrix{Data, IndexType})::IndexType

Gets the number of nonzero elements in the matrix
"""
function Base.nnz(matrix::LocalCSRMatrix{Data, IndexType}) where {Data, IndexType}
    IndexType(length(matrix.values))
end

"""
    numRows(::LocalCSRMatrix{Data, IndexType})::IndexType

Gets the number of rows in the matrix
"""
numRows(matrix::LocalCSRMatrix) = numRows(matrix.graph)

"""
    numCols(::LocalCSRMatrix{Data, IndexType})::IndexType

Gets the number of columns in the matrix
"""
numCols(matrix::LocalCSRMatrix) = matrix.numCols

"""
    getRowView((matrix::LocalCSRMatrix{Data, IndexType}, row::Integer)::SparseRowView{Data, IndexType}

Gets a view of the requested row
"""
function getRowView(matrix::LocalCSRMatrix{Data, IndexType},
        row::Integer)::SparseRowView{Data, IndexType} where {Data, IndexType}
    start = matrix.graph.rowMap[row]
    count = matrix.graph.rowMap[row+1] - start
    
    if count == 0
        SparseRowView(Data[], IndexType[])
    else
        SparseRowView(matrix.values, matrix.graph.entries, count, start)
    end
end



