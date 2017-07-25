
#TODO stop using assert for routine validation
#TODO test errors

#### Test BlockMap with SerialComm ####

macro SerialMapTests()
    quote
        mapCopy = BlockMap(map)
        
        @test uniqueGIDs(map)

        for i = 1:5
            @test myLID(map, i)
            @test myGID(map, i)
            @test i == lid(map, i)
            @test i == gid(map, i)
        end

        @test !myLID(map, -1)
        @test !myLID(map, 0)
        @test !myLID(map, 6)
        @test !myLID(map, 30)
        @test !myGID(map, -1)
        @test !myGID(map, 0)
        @test !myGID(map, 6)
        @test !myGID(map, 30)

        @test 0 == lid(map, -1)
        @test 0 == lid(map, 0)
        @test 0 == lid(map, 6)
        @test 0 == lid(map, 30)
        @test 0 == gid(map, -1)
        @test 0 == gid(map, 0)
        @test 0 == gid(map, 6)
        @test 0 == gid(map, 30)

        @test !distributedGlobal(map)

        @test 5 == numGlobalElements(map)
        @test 5 == numMyElements(map)

        @test 1 == minMyGID(map)
        @test 5 == maxMyGID(map)
        @test 1 == minAllGID(map)
        @test 5 == maxAllGID(map)
        @test 1 == minLID(map)
        @test 5 == maxLID(map)

        @test ([1, 1, 1, 1, 1], [1, 2, 3, 4, 5]) == remoteIDList(map, [1, 2, 3, 4, 5])

        @test [1, 2, 3, 4, 5] == myGlobalElements(map)

        @test sameBlockMapDataAs(map, mapCopy)
        @test sameBlockMapDataAs(mapCopy, map)
        @test !sameBlockMapDataAs(map, map2)
        @test !sameBlockMapDataAs(map2, map)
        @test !sameBlockMapDataAs(map, diffMap)
        @test !sameBlockMapDataAs(diffMap, map)

        @test sameAs(map, mapCopy)
        @test sameAs(mapCopy, map)
        @test sameAs(map, map2)
        @test sameAs(map2, map)
        @test !sameAs(map, diffMap)
        @test !sameAs(diffMap, map)
        @test !sameAs(map2, diffMap)
        @test !sameAs(diffMap, map2)

        @test linearMap(map)

        @test [1, 2, 3, 4, 5] == myGlobalElementIDs(map)

        @test commVal == comm(map)
    end
end

commVal = SerialComm()


## constructor 1 ##
map = BlockMap(5, commVal)
map2 = BlockMap(5, commVal)
diffMap = BlockMap(6, commVal)

@SerialMapTests

## constructor 2 ##
map = BlockMap(5, 5, commVal)
map2 = BlockMap(5, 5, commVal)
diffMap = BlockMap(6, 6, commVal)

@SerialMapTests

## constructor 3 ##
map = BlockMap(5, 5, [1, 2, 3, 4, 5], commVal)
map2 = BlockMap(5, 5, [1, 2, 3, 4, 5], commVal)
diffMap = BlockMap(6, 6, [1, 2, 3, 4, 5, 6], commVal)

@SerialMapTests

## constructor 4 ##

map = BlockMap(5, 5, [1, 2, 3, 4, 5], false, 1, 5, commVal)
map2 = BlockMap(5, 5, [1, 2, 3, 4, 5], false, 1, 5, commVal)
diffMap = BlockMap(6, 6, [1, 2, 3, 4, 5, 6], false, 1, 6, commVal)

@SerialMapTests

#TODO test BlockMap with parallel Comm