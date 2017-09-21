
#TODO write MPI tests

#TODO ensure result of CSRMatrix(rowMap, colMap, localMatrix, plist) is fill complete

#TODO make testing version of checkInternalState to run during testing

#### Serial Tests####s

n = 8
m = 6

Data = Float32
GID = UInt16
PID = Bool
LID = UInt8

commObj = SerialComm{GID, PID, LID}()
rowMap = BlockMap(n, n, commObj)


mat = CSRMatrix{Data}(rowMap, m, STATIC_PROFILE)
@test isa(mat, CSRMatrix{Data, GID, PID, LID})

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
#skipped many of the global methods because those require re-generating and may not be up to date

row = getGlobalRowCopy(mat, 2)
@test isa(row, Tuple{<: AbstractArray{GID, 1}, <: AbstractArray{Data, 1}})
@test GID[1, 3, 4] == row[1]
@test Data[2.5, 6.21, 77] == row[2]

row = getGlobalRowView(mat, 2)
@test isa(row, Tuple{<: AbstractArray{GID, 1}, <: AbstractArray{Data, 1}})
@test GID[1, 3, 4] == row[1]
@test Data[2.5, 6.21, 77] == row[2]

println("\n\nBefore filling:")
println("mat.myGraph.localIndices1D = $(mat.myGraph.localIndices1D)")
println("mat.myGraph.globalIndices1D = $(mat.myGraph.globalIndices1D)")
println("mat.myGraph.rowOffsets = $(mat.myGraph.rowOffsets)")

println("mat.myGraph.localIndices2D = $(mat.myGraph.localIndices2D)")
println("mat.myGraph.globalIndices2D = $(mat.myGraph.globalIndices2D)")
println("mat.myGraph.numRowEntries = $(mat.myGraph.numRowEntries)")

fillComplete(mat)
println("\n\nafter filling:")
println("mat.myGraph.localIndices1D = $(mat.myGraph.localIndices1D)")
println("mat.myGraph.globalIndices1D = $(mat.myGraph.globalIndices1D)")
println("mat.myGraph.rowOffsets = $(mat.myGraph.rowOffsets)")

println("mat.myGraph.localIndices2D = $(mat.myGraph.localIndices2D)")
println("mat.myGraph.globalIndices2D = $(mat.myGraph.globalIndices2D)")
println("mat.myGraph.numRowEntries = $(mat.myGraph.numRowEntries)")

println("getting local row copy")
row = getLocalRowCopy(mat, 2)
@test isa(row, Tuple{<: AbstractArray{LID, 1}, <: AbstractArray{Data, 1}})
@test LID[1, 2, 3] == row[1]
@test Data[2.5, 6.21, 77] == row[2]

println("getting local row view")
row = getLocalRowView(mat, 2)
@test isa(row, Tuple{<: AbstractArray{LID, 1}, <: AbstractArray{Data, 1}})
@test LID[1, 2, 3] == row[1]
@test Data[2.5, 6.21, 77] == row[2]

#=

getGlobalNumCols(mat::CSRMatrix) = -1#TODO figure out
getLocalNumCols(mat::CSRMatrix) = numCols(mat.localMatrix)
=#

map = BlockMap(2, 2, commObj)

mat = CSRMatrix{Data}(map, 2, STATIC_PROFILE)
insertGlobalValues(mat, 1, LID[1, 2], Data[2, 3])
insertGlobalValues(mat, 2, LID[1, 2], Data[5, 7])
fillComplete(mat)

Y = MultiVector(map, diagm(Data(1):2))
X = MultiVector(map, fill(Data(2), 2, 2))

apply!(Y, mat, X, NO_TRANS, Float32(3), Float32(.5))

@assert fill(2, 2, 2) == X #ensure X isn't mutated
exp = Array{Data, 2}(2, 2)
exp[:, 1] = [30.5, 30]
exp[:, 2] = [72,   73]
@assert exp == Y.data



Y = MultiVector(map, diagm(Data(1):2))
#X = MultiVector(map, fill(Data(2), 2, 2))

apply!(Y, mat, X, NO_TRANS, Float32(3), Float32(.5))

@assert fill(2, 2, 2) == X #ensure X isn't mutated
exp = Array{Data, 2}(2, 2)
exp[:, 1] = [42.5, 42]
exp[:, 2] = [20,   21]
@assert exp == Y.data
