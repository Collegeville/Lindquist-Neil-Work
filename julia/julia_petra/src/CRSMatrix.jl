
#TODO get a good understanding of whats going on for CRSMatrix and CRSGraph
#TODO make sure to add constructor to convert from Julia's sparce matrix

type CRSMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistRowMatrix{Data, GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::BlockMap{GID, PID, LID}
    
    #row start indexes (in local index space)
    IA::Array{LID, 1}
    #non zero entries
    A::Array{Data, 1}
    #columns (in local index space)
    JA::Array{LID, 1}
    
    numMyCols::LID
    numMyRows::LID
    numMyNonzeros::LID
    numMyDiagonals::LID
    
    numGlobalCols::GID
    numGlobalRows::GID
    numGlobalNonzeros::GID
    numGlobalDiagonals::GID
end

function CRSMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        A::Array{Data, 1}, IA::Array{LID, 1}, JA::Array{LID, 1}, numMyCols::Integer,
        numMyRows::Integer, numMyNonzeros::Integer, numMyDiagonals::Integer) where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    values = sumAll(comm(rowMap), [numMyRows, numMyNonzeros, numMyDiagonals])
    numGlobalRows = values[1]
    numGlobalNonZeros = values[2]
    numGlobalDiagonals = values[3]
end


#TODO constructors

function map(matrix::CRSMatrix)
    matrix.map
end

function numGlobalNonzeros(matrix::CRSMatrix{Data, GID})::GID where {
        Data <: Number, GID <: Integer}
    matrix.numGlobalNonzeros
end

function numGlobalCols(matrix::CRSMatrix{Data, GID})::GID where {
        Data <: Number, GID <: Integer}
    matrix.numGlobalCols
end

function numGlobalRows(matrix::CRSMatrix{Data, GID})::GID where {
        Data <: Number, GID <: Integer}
    matrix.numGlobalRows
end

function numGlobalDiagonals(matrix::CRSMatrix{Data, GID})::GID where {
        Data <: Number, GID <: Integer}
    return matrix.numGlobalDiagonals
end

function numLocalNonzeros(matrix::CRSMatrix{Data, GID, PID, LID})::LID where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    matrix.numLocalNonzeros
end

function numLocalCols(matrix::CRSMatrix{Data, GID, PID, LID})::LID where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    matrix.numLocalCols
end

function numLocalRows(matrix::CRSMatrix{Data, GID, PID, LID})::LID where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    matrix.numLocalRows
end

function numMyRowEntries(matrix::CRSMatrix{Data, GID, PID, LID},
        myRow::LID)::LID where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    matrix.IA[myRow+1] - matrix.IA[myRow]
end

function numMyDiagonals(matrix::CRSMatrix{Data, GID})::GID where {
        Data <: Number, GID <: Integer}
    return matrix.numMyDiagonals
end

function maxNumEnrtries(matrix::CRSMatrix{Data, GID, PID, LID})::LID where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    best = -1
    start = matrix.IA[1]
    for i = 2:length(matrix.IA)
        last = matrix.IA[i]
        dif = last-start
        if best < dif
            best = dif
        end
        start = last
    end
    best
end

function extractMyRowCopy(matrix::CRSMatrix{Data, GID, PID, LID}, myRow::LID
        )::Tuple{Array{Data}, Array{LID}} where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    start = matrix.IA[myRow]
    last = matrix.IA[myRow+1]-1
    
    (matrix.A[start, last], matrix.JA[start, last])
end

function rowMap(matrix::CRSMatrix{Data, GID, PID, LID})::BlockMap{GID, PID, LID} where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    matrix.rowMap
end

function colMap(matrix::CRSMatrix{Data, GID, PID, LID})::BlockMap{GID, PID, LID} where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    matrix.colMap
end
    
function leftScale!(matrix::CRSMatrix{Data, GID, PID, LID}, X::Array{Data}) where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    start = matrix.IA[1]
    for i in 1:numMyRows(matrix)
        last = matrix.IA[myRow+1]
        matrix.A[start:last-1] *= X[i]
    end
    matrix
end

function rightScale!(matrix::CRSMatrix{Data, GID, PID, LID}, X::Array{Data}) where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    for i = 1:numMyNonzeros(matrix)
        matrix.A[i] *= X[matrix.JA[i]]
    end
    matrix
end

function checkSize(source::CRSMatrix{Data, GID, PID, LID},
        target::CRSMatrix{Data, GID, PID, LID})::Bool where {
        Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    #DECISION are these checks really nessacery? can easily resize arrays
    (source.numMyRows == target.numMyRows
        && source.numMyRows == target.numMyRows
        && source.numGlobalRows == target.numGlobalRows
        && source.numGlobalCols == target.numGlobalCols)
end

#TODO implement all theses methods
#=
extractDiagonalCopy(matrix::Impl{Data, GID, PID, LID})::Array{Data}
    - Returns a copy of the main diagonal

isLowerTriangle(matrix::Impl)::Bool
    - Returns true if the matrix is lower triagular in local index space, otherwise false

isUpperTriangle(matrix::Impl)::Bool
    - Returns true if the matrix is upper triagular in local index space, otherwise false


apply!(matrix::Impl{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, mode::TransposeMode, alpha::Data, beta::Data)
    Computes ``Y = α\cdot A^{mode}\cdot X + β\cdot Y``, with the following exceptions
        If beta == 0, apply MUST overwrite Y, so that any values in Y (including NaNs) are ignored.
        If alpha == 0, apply MAY short-circuit the operator, so that any values in X (including NaNs) are ignored

domainMap(operator::Op{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    Returns the BlockMap associated with the domain of this operation

rangeMap(operator::Op{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    Returns the BlockMap associated with the range of this operation

copyAndPermute(source::???{GID, PID, LID}, target::Impl{GID, PID, LID},
        numSameIDs::LID, permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1})
    Perform copies and permutations that are local to this process.

packAndPrepare(source::???{GID, PID, LID}, target::Impl{GID, PID, LID},
        exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID}
        )::Array
    Perform any packing or preparation required for communications.  The
    method returns the array of objects to export

unpackAndCombine(target::Impl{GID, PID, LID}, importLIDs::Array{LID, 1},
        imports::Array, distor::Distributor{GID, PID, LID},
        cm::CombineMode)
    Perform any unpacking and combining after communication
=#
