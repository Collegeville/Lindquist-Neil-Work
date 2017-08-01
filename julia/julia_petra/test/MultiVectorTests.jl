#number of elements in vectors
n = 8

# use serial comm for local tests
serialComm = SerialComm{Int32, Bool, Int16}()
curMap = BlockMap(n, n, serialComm)

# test basic construction with setting data to zeros
vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, 3, true)
@test n == localLength(vect)
@test n == globalLength(vect)
@test 3 == numVectors(vect)
@test curMap == julia_petra.map(vect)
@test zeros(Float64, (n, 3)) == vect.data

# test basic construction without setting data to zeros
vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, 3, false)
@test n == localLength(vect)
@test n == globalLength(vect)
@test 3 == numVectors(vect)
@test curMap == julia_petra.map(vect)

# test wrapper constructor
arr = Array{Float64, 2}(n, 3)
vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, arr)
@test n == localLength(vect)
@test n == globalLength(vect)
@test 3 == numVectors(vect)
@test curMap == julia_petra.map(vect)
@test arr === vect.data

# test copy
vect2 = copy(vect)
@test n == localLength(vect)
@test n == globalLength(vect)
@test 3 == numVectors(vect)
@test curMap == julia_petra.map(vect)
@test vect.data == vect2.data
@test vect.data !== vect2.data #ensure same contents, but different address


# test scale and scale!
vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, ones(Float64, n, 3))
@test vect == scale!(vect, 5.0)
@test 5*ones(Float64, (n, 3)) == vect.data

vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, ones(Float64, n, 3))
vect2 = scale(vect, 5.0)
@test vect !== vect2
@test 5*ones(Float64, (n, 3)) == vect2.data

vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, ones(Float64, n, 3))
@test vect == scale!(vect, [2.0, 3.0, 4.0])
@test hcat(2*ones(Float64, n), 3*ones(Float64, n), 4*ones(Float64, n))  == vect.data

vect = MultiVector{Float64, Int32, Bool, Int16}(curMap, ones(Float64, n, 3))
vect2 = scale(vect, [2.0, 3.0, 4.0])
@test vect !== vect2
@test hcat(2*ones(Float64, n), 3*ones(Float64, n), 4*ones(Float64, n))  == vect2.data

#TODO create MPIMultiVectorTests that get included with the other MPI tests