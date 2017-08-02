n = 4

srcMap = BlockMap(4*n, comm)
desMap = BlockMap(4*n, n, collect((1:n) + n*(pid%4)), comm)

function basicMPITest(impor, debugSetting=false)
    if isa(impor, Import)
        data = impor.importData
    else
        data = impor.exportData
    end
    @test debugSetting == impor.debug
    @test srcMap == data.source
    @test desMap == data.target
    @test 0 == data.numSameIDs
    @test isa(data.distributor, Distributor{UInt64, UInt16, UInt32})
    @test [] == data.permuteToLIDs
    @test [] == data.permuteFromLIDs
    
    @test true == data.isLocallyComplete
end

# basic import
basicMPITest(Import(srcMap, desMap))
basicMPITest(Import(srcMap, desMap, debug=true), true)
basicMPITest(Import(srcMap, desMap, Nullable{Array{UInt16}}()))
basicMPITest(Import(srcMap, desMap, Nullable{Array{UInt16}}(), debug=true), true)

# import using Dicts
basicMPITest(Import(srcMap, desMap, Dict{Symbol, Any}()))
basicMPITest(Import(srcMap, desMap, Dict([(:debug, true)])), true)
basicMPITest(Import(srcMap, desMap, Nullable{Array{UInt16}}(), Dict{Symbol, Any}()))
basicMPITest(Import(srcMap, desMap, Nullable{Array{UInt16}}(), Dict([(:debug, true)])), true)


# basic export
basicMPITest(Export(srcMap, desMap))
basicMPITest(Export(srcMap, desMap, debug=true), true)
basicMPITest(Export(srcMap, desMap, Nullable{Array{UInt16}}()))
basicMPITest(Export(srcMap, desMap, Nullable{Array{UInt16}}(), debug=true), true)

#export using Dicts
basicMPITest(Export(srcMap, desMap, Dict{Symbol, Any}()))
basicMPITest(Export(srcMap, desMap, Dict([(:debug, true)])), true)
basicMPITest(Export(srcMap, desMap, Nullable{Array{UInt16}}(), Dict{Symbol, Any}()))
basicMPITest(Export(srcMap, desMap, Nullable{Array{UInt16}}(), Dict([(:debug, true)])), true)