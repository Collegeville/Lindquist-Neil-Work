
export Export 

"""
Communication plan for data redistribution from a (possibly) multipl-owned to a uniquely owned distribution
"""
type Export{GID <: Integer, PID <:Integer, LID <: Integer}
    debug::Bool
    exportData::ImportExportData{GID, PID, LID}
end

## Constructors ##

function Export{GID, PID, LID}(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID}, remotePIDs::Nullable{Array{PID}}=Nullable{Array{PID}}; plist...) where {GID <: Integer, PID <: Integer, LID <: Integer}
    Export{GID, PID, LID}(source, target, remotePIDs, Dict(plist))
end

function Export{GID, PID, LID}(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID}, plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    const debug = get(plist, :debug, false)

    if debug
        print("$(myPid(comm(source))): Export ctor\n")
    end

    expor = Export{GID, PID, LID}(debug, ImportExportData(source, target))

    exportGIDs = setupSamePermuteExport(expor)

    if debug
        print("$(myPid(comm(source))): Export ctor: setupSamePermuteExport done\n")
    end
    if isDistributed(source)
        setupRemote(expor, remoteGIDs, remotePIDs)
    end
    if debug
        print("$(myPid(comm(source))): Export ctor: done\n")
    end

    impor
end


## internal construction methods ##

function setupSamePermuteExport(expor::Export{GID, PID, LID})::Array{GID} where {GID <: Integer, PID <: Integer, LID <: Integer}
    
    data = expor.exportData
    
    source = sourceMap(expor)
    target = targetMap(expor)
    
    sourceGIDs = myGlobalElements(source)
    targetGIDs = myGlobalElements(target)
    
    numSrcGIDs = length(sourceGIDs)
    numTgtGIDs = length(targetGIDs)
    numGIDs = min(numSrcGIDs, numTgtGIDs)
    
    numSameGIDs = 1
    while numSameGIDs <= numGIDs && rawSrcGIDs[numSameGIDs] == rawTgtGIDs[numSameGIDs]
        numSameGIDs += 1
    end
    numSameGIDs -= 1
    numSameIDs(data, numSameGIDs)
    
    exportGIDs = Array{GID, 1}(0)
    permuteToLIDs = permuteToLIDs(data)
    permuteFromLIDs = permuteFromLIDs(data)
    exportLIDs = exportLIDs(data)
    
    for srcLID = (numSameGIDs+1):numSrcLIDs
        const curSrcGID = sourceGIDs[srcLID]
        const tgtLID = lid(target, curSrcGID)
        if tgtLID != 0
            push!(permuteToLIDs, tgtLID)
            push!(permuteFromLIDs, srcLID)
        else
            push!(exportGIDs, curSrcGID)
            push!(exportLIDs, srcLID)
        end
    end
    
    if length(exportLIDs) != 0 && !isDistributed(source)
        isLocallyComplete(data, false)
        warn("Source has export LIDs but source not distributed globally.  " *
            "Exporting to a submap of the target map.")
    end
    
    if isDistributed(source)
        resize!(exportPIDs(data), length(exportGIDs))
        
        (exportPIDs, exportLIDs) = remoteIDList(target, exportGIDs)
        julia_petra.exportPIDs(data, exportPIDs)
        missingGIDs = 0
        for i = 1:length(exportPIDs)
            if exportPIDs[i] == 0
                missingGIDs += 1
            end
        end
        
        if missingGIDs != 0
            warn("The source Map has GIDs not found in the target Map")
            
            isLocallyComplete(data, false)
            numInvalidExports = missingGIDs
            totalNumExports = length(exportPIDs)
            
            if numInvalidExports == totalNumExports
                # all exports invalid, can delete all exports
                resize!(exportGIDs, 0)
                resize!(exportLIDs, 0)
                resize!(exportPIDs, 0)
            else
                #some exports are valid, need to keep the valid exports
                numValidExports = 1
                for e = 1:totalNumExports
                    if exportPIDs[e] != 0
                        exportGIDs[numValidExports] = exportGIDs[e]
                        exportLIDs[numValidExports] = exportLIDs[e]
                        exportPIDs[numValidExports] = exportPIDs[e]
                        numValidExports += 1
                    end
                end
                numValidExports -= 1

                resize!(exportGIDs, numValidExports)
                resize!(exportLIDs, numValidExports)
                resize!(exportPIDs, numValidExports)
            end
        end
    end
    exportGIDs
end
    
function setupRemote(expor::Export{GID, PID, LID}, exportGIDs::Array{GID, 1}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    
    data = expor.exportData
    
    target = targetMap(data)

    if expor.debug
        print("$(myPid(comm(target))): setupRemote\n")
    end
    
    order = sortperm(exportPIDs(data))
    permute!(exportPIDs(data), order)
    permute!(exportLIDs(data), order)
    permute!(exportGIDs, order)
    
    if expor.debug
        print("$(myPid(comm(target))): setupRemote: Calling createFromSends\n")
    end
    
    numRemoteIDs = createFromSends(distributor(data), exportPIDs(data))
    
    if expor.debug
        print("$(myPid(comm(target))): setupRemote: Calling doPostsAndWaits\n")
    end
    
    remoteGIDs = resolve(distributor(data), exportGIDs)
    
    remoteLIDs = remoteLIDs(data)
    
    resize!(remoteLIDs, numRemoteIDs)
    
    for (i, j) in zip(length(remoteGIDs), length(remoteLIDs(data)))
        remoteLIDs[j] = lid(target, remoteGIDs[i])
    end
    
    if expor.debug
        print("$(myPid(comm(target))): setupRemote: done\n")
    end
end
    
## Getters ##

"""
Get the source map for the given ImportExportData
"""
function sourceMap(impor::Export{GID, PID, LID})::BlockMap{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.importData.source
end

"""
Get the target map for the given ImportExportData
"""
function targetMap(impor::Export{GID, PID, LID})::BlockMap{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.importData.target
end

"""
List of elements in the target map that are permuted.
"""
function permuteToLIDs(impor::Export{GID, PID, LID})::Array{LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.importData.permuteToLIDs
end

"""
List of elements in the source map that are permuted.
"""
function permuteFromLIDs(impor::Export{GID, PID, LID})::Array{LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.importData.permuteFromLIDs
end

"""
List of elements in the target map that are coming from other processors
"""
function remoteLIDs(impor::Export{GID, PID, LID})::Array{LID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.importData.remoteLIDs
end

"""
List of elements that will be sent to other processors
"""
function exportLIDs(impor::Export{GID, PID, LID})::Array{LID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.importData.exportLIDs
end

"""
List of processors to which elements will be sent `exportLID[i]` will be sent to processor `exportPIDs[i]`
"""
function exportPIDs(impor::Export{GID, PID, LID})::Array{PID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.importData.exportPIDs
end

"""
Returns the number of elements that are identical between the source and target maps, up to the first different ID
"""
function numSameIDs(impor::Export{GID, PID, LID})::LID where GID <: Integer where PID <: Integer where LID <: Integer
    impor.importData.numSameIDs
end


"""
Returns the distributor being used
"""
function distributor(impor::Export{GID, PID, LID})::Distributor{GID, PID, LID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.importData.distributor
end

"""
Returns whether the import or export is locally complete
"""
function isLocallyComplete(impor::Export{GID, PID, LID})::Bool where GID <: Integer where PID <: Integer where LID <: Integer
    impor.importData.isLocallyComplete
end