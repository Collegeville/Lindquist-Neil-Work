

"""
The version of RowMatrix that isn't a subtype of DestObject
"""
abstract type SrcDestRowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: SrcDestObject{GID, PID, LID} 
end

"""
The version of RowMatrix that is a subtype of DestObject
"""
abstract type DestRowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DestObject{GID, PID, LID} 
end

#TODO figure out the needed operation functions (multiple, solve, norms, ect), (lines 123-221)

"""
RowMatrix is actually a type union of SrcDestRowMatrix and DestRowMatrix,
which are (direct) subtypes of SrcDestObject and DestObject, respectively.

All subtypes must have the following methods, with Impl standing in for the subtype:
numMyRowEntries(matrix::Impl{Data, GID, PID, LID}, myRow::LID}::LID
    - Returns the number of nonzero entries in row `myRow`

maxNumEntries(matrix::Impl{Data, GID, PID, LID})::LID
    - Returns the maximum of `numMyRowEntries()` over all rows

extractMyRowCopy(matrix::Impl{Data, GID, PID, LID}, myRow::LID)::Tuple{Array{Data}, Array{LID}}
    - Returns the contents of the row `myRow` as a Tuple containing
        1 - the values of the row
        2 - the indices where each element belongs

extractDiagonalCopy(matrix::Impl{Data, GID, PID, LID})::Array{Data}
    - Returns a copy of the main diagonal

numGlobalNonzeros(matrix::Impl{Data, GID})::GID
    - Returns the number of nonzero elements on all processors.
        Note that depending on the matrix implementation, it is sometimes
        possible to have some nonzeros that appear on multiple processors.
        In that case, those nonzeros may be counted multiple times (also
        depending on the matrix implementation)

numGobalRows(matrix::Impl{Data, GID})::GID
    - Returns the number of rows on all processors

numGlobalCols(matrix::Impl{Data, GID})::GID
    - Returns the number of columns on all processors

numGlobalDiagonals(matrix::Impl{Data, GID})::GID
    - Returns the number of nonzero diagonal entries on all processors

numMyNonzeros(matrix::Impl{Data, GID})::GID
    - Returns the number of nonzero elements on the calling processor.

numMyRows(matrix::Impl{Data, GID})::GID
    - Returns the number of rows on calling processor

numMyCols(matrix::Impl{Data, GID})::GID
    - Returns the number of columns on calling processor

numMyDiagonals(matrix::Impl{Data, GID})::GID
    - Returns the number of nonzero diagonal entries on calling processor

isLowerTriangle(matrix::Impl)::Bool
    - Returns true if the matrix is lower triagular in local index space, otherwise false

isUpperTriangle(matrix::Impl)::Bool
    - Returns true if the matrix is upper triagular in local index space, otherwise false

rowMap(matrix::Impl{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    - Returns the BlockMap associated with the rows of this matrix

colMap(matrix::Impl{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    - Returns the BlockMap associated with the columns of this matrix
"""
const RowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} = Union{SrcDestRowMatrix{Data <: Number, GID, PID, LID}, DestRowMatrix{GID, PID, LID}}