
#a few light tests to catch basic issues

commObj = SerialComm{UInt32, UInt8, UInt16}()
map = BlockMap(20, commObj)

graph = CRSGraph(map, UInt16(15), STATIC_PROFILE, Dict{Symbol, Any}())
@test julia_petra.map(graph) == map
@test graph.pftype == STATIC_PROFILE

commObj = SerialComm{UInt8, Int8, UInt16}()
map = BlockMap(20, commObj)

@test_throws InvalidArgumentError CRSGraph(map, UInt16(15), STATIC_PROFILE, Dict{Symbol, Any}())
