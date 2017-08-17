
#TODO write serial tests
#TODO write MPI tests

#TODO ensure result of CSRMatrix(rowMap, colMap, localMatrix, plist) is fill complete

#TODO make testing version of checkInternalState to run during testing

n = 8
m = 6

commObj = SerialComm{UInt16, Bool, UInt8}()
rowMap = BlockMap(n, n, commObj)

mat = CSRMatrix{Float32}(rowMap, m, STATIC_PROFILE, Dict{Symbol, Any}())
@test isa(mat, CSRMatrix{Float32, UInt16, Bool, UInt8})
@test STATIC_PROFILE == getProfileType(mat)
@test isFillActive(mat)
@test !isFillComplete(mat)
@test isGloballyIndexed(mat)
@test !isLocallyIndexed(mat)
@test rowMap == getRowMap(mat)
@test !hasColMap(mat)

@test 0 == mat.myGraph.nodeNumEntries
insertGlobalValues(mat, 2, UInt8[1, 3, 4], Float32[2.5, 6.21, 77])
@test 3 == mat.myGraph.nodeNumEntries


