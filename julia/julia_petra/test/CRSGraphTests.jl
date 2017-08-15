
#a few light tests to catch basic issues

commObj = SerialComm{UInt32, UInt8, UInt16}()
map = BlockMap(20, commObj)

function basicTests(graph)
    @test isLocallyIndexed(graph)
    @test !isGloballyIndexed(graph)
    @test isFillActive(graph)
    @test !isFillComplete(graph)
    @test !hasColMap(graph)
end

graph = CRSGraph(map, UInt16(15), STATIC_PROFILE, Dict{Symbol, Any}())
@test map == julia_petra.map(graph)
@test STATIC_PROFILE == getProfileType(graph)
basicTests(graph)

graph = CRSGraph(map, UInt16(15), STATIC_PROFILE, Dict{Symbol, Any}(:debug=>true))
@test map == julia_petra.map(graph)
@test STATIC_PROFILE == getProfileType(graph)
basicTests(graph)

graph = CRSGraph(map, UInt16(15), DYNAMIC_PROFILE, Dict{Symbol, Any}(:debug=>true))
@test map == julia_petra.map(graph)
@test DYNAMIC_PROFILE == getProfileType(graph)
basicTests(graph)


commObj = SerialComm{UInt8, Int8, UInt16}()
map = BlockMap(20, commObj)

@test_throws InvalidArgumentError CRSGraph(map, UInt16(15), STATIC_PROFILE, Dict{Symbol, Any}())
