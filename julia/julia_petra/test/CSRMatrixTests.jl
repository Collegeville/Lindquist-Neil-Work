
#TODO write serial tests
#TODO write MPI tests

#TODO ensure result of CSRMatrix(rowMap, colMap, localMatrix, plist) is fill complete

#TODO make testing version of checkInternalState to run during testing

n = 8
m = 6

Data = Float32
GID = UInt16
PID = Bool
LID = UInt8

commObj = SerialComm{GID, PID, LID}()
rowMap = BlockMap(n, n, commObj)

mat = CSRMatrix{Data}(rowMap, m, STATIC_PROFILE, Dict{Symbol, Any}())
@test isa(mat, CSRMatrix{Data, GID, PID, LID})
@test STATIC_PROFILE == getProfileType(mat)
@test isFillActive(mat)
@test !isFillComplete(mat)
@test isGloballyIndexed(mat)
@test !isLocallyIndexed(mat)
@test rowMap == getRowMap(mat)
@test !hasColMap(mat)
@test n == getGlobalNumRows(mat)
@test n == getLocalNumRows(mat)


@test 0 == getNumEntriesInLocalRow(mat, 2)
@test 0 == getNumEntriesInGlobalRow(mat, 2)
@test 0 == getLocalNumEntries(mat)
@test 0 == getGlobalNumEntries(mat)
@test 0 == getGlobalNumDiags(mat)
@test 0 == getLocalNumDiags(mat)
@test 0 == getGlobalMaxNumRowEntries(mat)
@test 0 == getLocalMaxNumRowEntries(mat)
insertGlobalValues(mat, 2, LID[1, 3, 4], Data[2.5, 6.21, 77])
@test 3 == getNumEntriesInLocalRow(mat, 2)
@test 3 == getNumEntriesInGlobalRow(mat, 2)
@test 3 == getLocalNumEntries(mat)
@test 0 == getLocalNumDiags(mat)
@test 3 == julia_petra.getRowInfo(mat.myGraph, LID(2)).numEntries
#skipped many of the global methods because that require re-generating and may not be up to date

row = getGlobalRowCopy(mat, 2)
@test isa(row, Tuple{<: AbstractArray{GID, 1}, <: AbstractArray{Data, 1}})
@test GID[1, 3, 4] == row[1]
@test Data[2.5, 6.21, 77] == row[2]

row = getGlobalRowView(mat, 2)
@test isa(row, Tuple{<: AbstractArray{GID, 1}, <: AbstractArray{Data, 1}})
@test GID[1, 3, 4] == row[1]
@test Data[2.5, 6.21, 77] == row[2]


#=

getGlobalNumCols(mat::CSRMatrix) = -1#TODO figure out
getLocalNumCols(mat::CSRMatrix) = numCols(mat.localMatrix)
=#
