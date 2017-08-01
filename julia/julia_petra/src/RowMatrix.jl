

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
RowMatrix is the base of 

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

leftScale!(matrix::Impl{Data, GID, PID, LID}, X::PetraVector{Data, GID, PID, LID})
    Scales matrix on the left with X

rightScale!(matrix::Impl{Data, GID, PID, LID}, X::PetraVector{Data, GID, PID, LID})
    Scales matrix on the right with X


Additionally, the following method must be implemented to fufil the operator interface:

apply!(matrix::Impl{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, mode::TransposeMode, alpha::Data, beta::Data)
    Computes ``Y = α\cdot A^{mode}\cdot X + β\cdot Y``, with the following exceptions
        If beta == 0, apply MUST overwrite Y, so that any values in Y (including NaNs) are ignored.
        If alpha == 0, apply MAY short-circuit the operator, so that any values in X (including NaNs) are ignored

domainMap(operator::Op{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    Returns the BlockMap associated with the domain of this operation

rangeMap(operator::Op{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    Returns the BlockMap associated with the range of this operation
"""
const RowMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} = Union{SrcDistRowMatrix{Data, GID, PID, LID}, DistRowMatrix{Data, GID, PID, LID}}

@operatorFunctions RowMatrix
