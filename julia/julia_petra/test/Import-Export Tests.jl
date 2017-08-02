#mainly just exercies the basic constructors

n = 8

serialComm = SerialComm{Int32, Bool, Int16}()
srcMap = BlockMap(n, n, serialComm)
desMap = BlockMap(n, n, serialComm)

function basicTest(impor, debugSetting=false)
    if isa(impor, Import)
        data = impor.importData
    else
        # basic exports are about the same anyways
        data = impor.exportData
    end
    @test debugSetting == impor.debug
    @test srcMap == data.source
    @test desMap == data.target
    @test n == data.numSameIDs
    @test isa(data.distributor, Distributor{Int32, Bool, Int16})
    @test [] == data.permuteToLIDs
    @test [] == data.permuteFromLIDs
    @test [] == data.remoteLIDs
    @test [] == data.exportLIDs
    @test [] == data.exportPIDs
    @test true == data.isLocallyComplete
end

# basic import
basicTest(Import(srcMap, desMap))
basicTest(Import(srcMap, desMap, debug=true), true)
basicTest(Import(srcMap, desMap, Nullable{Array{Bool}}()))
basicTest(Import(srcMap, desMap, Nullable{Array{Bool}}(), debug=true), true)

# import using Dicts
basicTest(Import(srcMap, desMap, Dict{Symbol, Any}()))
basicTest(Import(srcMap, desMap, Dict([(:debug, true)])), true)
basicTest(Import(srcMap, desMap, Nullable{Array{Bool}}(), Dict{Symbol, Any}()))
basicTest(Import(srcMap, desMap, Nullable{Array{Bool}}(), Dict([(:debug, true)])), true)


# basic export
basicTest(Export(srcMap, desMap))
basicTest(Export(srcMap, desMap, debug=true), true)
basicTest(Export(srcMap, desMap, Nullable{Array{Bool}}()))
basicTest(Export(srcMap, desMap, Nullable{Array{Bool}}(), debug=true), true)

#export using Dicts
basicTest(Export(srcMap, desMap, Dict{Symbol, Any}()))
basicTest(Export(srcMap, desMap, Dict([(:debug, true)])), true)
basicTest(Export(srcMap, desMap, Nullable{Array{Bool}}(), Dict{Symbol, Any}()))
basicTest(Export(srcMap, desMap, Nullable{Array{Bool}}(), Dict([(:debug, true)])), true)