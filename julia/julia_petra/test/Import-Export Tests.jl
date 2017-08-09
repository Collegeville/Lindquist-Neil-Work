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

#for scoping purposes
impor = Array{Import, 1}(1)
expor = Array{Export, 1}(1)

#ensure at least a few lines, each starting with the PID
debugregex = "^(?:$(myPid(serialComm)): .+\n){3,}\$"

# basic import
@test_nowarn impor[1] = Import(srcMap, desMap)
basicTest(impor[1])
@test_warn debugregex impor[1] = Import(srcMap, desMap, debug=true)
basicTest(impor[1], true)
@test_nowarn impor[1] = Import(srcMap, desMap, Nullable{Array{Bool}}())
basicTest(impor[1])
@test_warn debugregex impor[1] = Import(srcMap, desMap, Nullable{Array{Bool}}(), debug=true)
basicTest(impor[1], true)


# import using Dicts
@test_nowarn impor[1] = Import(srcMap, desMap, Dict{Symbol, Any}())
basicTest(impor[1])
@test_warn debugregex impor[1] = Import(srcMap, desMap, Dict([(:debug, true)]))
basicTest(impor[1], true)
@test_nowarn impor[1] = Import(srcMap, desMap, Nullable{Array{Bool}}(), Dict{Symbol, Any}())
basicTest(impor[1])
@test_warn debugregex impor[1] = Import(srcMap, desMap, Nullable{Array{Bool}}(), Dict([(:debug, true)])) 
basicTest(impor[1], true)


# basic export
@test_nowarn expor[1] = Export(srcMap, desMap)
basicTest(expor[1])
@test_warn debugregex expor[1] = Export(srcMap, desMap, debug=true)
basicTest(expor[1], true)
@test_nowarn expor[1] = Export(srcMap, desMap, Nullable{Array{Bool}}())
basicTest(expor[1])
@test_warn debugregex expor[1] = Export(srcMap, desMap, Nullable{Array{Bool}}(), debug=true)
basicTest(expor[1], true)

#export using Dicts
@test_nowarn expor[1] = Export(srcMap, desMap, Dict{Symbol, Any}())
basicTest(expor[1])
@test_warn debugregex expor[1] = Export(srcMap, desMap, Dict([(:debug, true)]))
basicTest(expor[1], true)
@test_nowarn expor[1] = Export(srcMap, desMap, Nullable{Array{Bool}}(), Dict{Symbol, Any}())
basicTest(expor[1])
@test_warn debugregex expor[1] = Export(srcMap, desMap, Nullable{Array{Bool}}(), Dict([(:debug, true)]))
basicTest(expor[1], true)