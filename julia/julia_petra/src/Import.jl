
"""
Communication plan for data redistribution from a uniquely-owned to a (possibly) multiply-owned distribution.
"""
type Import{GID <: Integer, PID <:Integer, LID <: Integer}
    debug::Bool
    importData::ImportExportData{GID, PID, LID}
end

## Constructors ##

function Import(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID},
        userRemotePIDs::Array{PID}, remoteGIDs::Array{GID}, userExportLIDs::Array{LID},
        userExportPIDs::Array{PID}, useRemotePIDGID::Bool,
        ; plist...) where {GID <: Integer, PID <: Integer, LID <: Integer}
    Import{GID, PID, LID}(source, target, userRemotePIDs, remoteGIDs,
        userExportLIDs, userExportPIDs, useREmotePIDGID, Dict(plist))
end

function Import(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID},
        userRemotePIDs::Array{PID}, remoteGIDs::Array{GID}, userExportLIDs::Array{LID},
        userExportPIDs::Array{PID}, useRemotePIDGID::Bool,
        plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    sourceGIDs = myGlobalElements(source)
    targetGIDs = myGlobalElements(target)
    numSrcGIDs = length(sourceGIDs)
    numTgtGIDs = length(targetGITs)
    numGIDS = min(numSrcGIDs, numTgtGIDs)
    
    numSameGIDs = 1
    while numSameGIDS < numGIDs && sourceGIDs[numSameGIDs] == targetGIDs[numSameGIDs]
        numSameGIDs += 1
    end
    numSameGIDs -= 1
    
    const debug = get(plist, :debug, false)
    
    if debug
        println("$(myPid(comm(source))): Import ctor expert\n")
    end
    
    importData = ImportExportData(source, target)
    numSameIDs(importData, numSameGIDS)
    
    permuteToLIDs = permuteToLIDs(importData)
    permuteFromLIDs = permuteFromLIDs(importData)
    remoteLIDs = remoteLIDs(importData)
    
    if !userRemotePIDGID
        empty!(remoteGIDs)
        empty!(remoteLIDs)
    end
    
    for tgtLID =  (numSameGIDs+1):numTgtLIDs
        const curTgtGID = targetGIDs[tgtLID]
        const srcLID = lid(curTgtGID)
        if srcLID != 0
            push!(permuteToLID, tgtLID)
            push!(permuteFromLID, srcLID)
        else
            if !userRemotePIDGID
                push!(remoteGIDs, curTgtGID)
                push!(remoteLIDs, tgtLID)
            end
        end
    end
    
    if length(remoteGIDs) > 0 && !isDistributed(source)
        throw(InvalidArgumentError("Target has remote LIDs but source is not distributed globally"))
    end
    
    (remotePIDs, _) = remoteIDList(source, remoteGIDs)
    
    remoteProcIDs = (useRemotePIDGID) ? userRemotePIDs : remotePIDs
    
    if !(length(remoteProcIDs) == length(remoteGIDs) && length(remoteGIDs) == length(remoteLIDs))
        throw(InvalidArgumentError("Size miss match on remoteProcIDs, remoteGIDs and remoteLIDs"))
    end
    
    # ensure remoteProcIDs[i], remoteGIDs[i] and remoteLIDs[i] refer to the same thing
    order = sortperm(remoteProcIDs)
    permute!(remoteProcIDs, order)
    permute!(remoteGIDs, order)
    permute!(remoteLIDs, order)
    
    exportPIDs = Array{PID, 1}(length(userExportPIDs))
    exportLIDs = Array{PID, 1}(length(userExportPIDs))
    
    #need the funcitons with these names, not the variables
    julia_petra.remoteLIDs(importData, remoteLIDs)
    julia_petra.exportPIDs(importData, exportPIDs)
    julia_petra.exportLIDs(importData, exportLIDs)
    
    locallyComplete = true
    for i = 1:length(userExportPIDs)
        if userExportPIDs[i] == 0
            locallyComplete = false
        end
        
        exportPIDs[i] = userExportPIDs[i]
        exportLIDs[i] = userExportLIDs[i]
    end
    
    isLocallyComplete(importData, locallyComplete)
    createFromSendsAndRecvs(distributor(importData), exportPIDs, remoteProcIDs)
    
    Import(debug, importData)
end

function Import{GID, PID, LID}(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID}, remotePIDs::Nullable{Array{PID}}=Nullable{Array{PID}}; plist...) where {GID <: Integer, PID <: Integer, LID <: Integer}
    Import{GID, PID, LID}(source, target, remotePIDs, Dict(plist))
end

function Import{GID, PID, LID}(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID}, plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    Import{GID, PID, LID}(source, target, Nullable{Array{PID}}(), plist)
end

function Import{GID, PID, LID}(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID}, remotePIDs::Nullable{Array{PID}}, plist::Dict{Symbol}) where {GID <: Integer, PID <: Integer, LID <: Integer}
        const debug = get(plist, :debug, false)
        
        if debug
            print("$(myPid(comm(source))): Import ctor\n")
        end
        
        impor = Import{GID, PID, LID}(debug, importData(source, target))
        
        remoteGIDs = setupSamePermuteRemote(impor)
        
        if debug
            print("$(myPid(comm(source))): Import ctor: setupSamePermuteRemote done\n")
        end
        if isDistributed(source)
            setupExport(impor, remoteGIDs, remotePIDs)
        end
        if debug
            print("$(myPid(comm(source))): Import ctor: done\n")
        end
        
        impor
    end

## internal construction methods ##

function setupSamePermuteRemote(impor::Import{GID, PID, LID}) where {GID <: Integer, PID <: Integer, LID <: Integer} 
    
    data = impor.data
    
    remoteGIDs = Array{GID, 1}(0)
    
    source = sourceMap(impor)
    target = targetMap(impor)
    
    sourceGIDs = myGlobalElements(source)
    targetGIDs = myGlobalElements(target)
    
    numSrcGIDs = length(sourceGIDs)
    numTgtGIDs = length(targetGIDs)
    numGIDs = min(numSrcGIDs, numTgtGIDs)
    
    numSameGIDs = 1
    while numSameGIDs < numGIDs && sourceGIDs[numSameGIDs] == targetGIDs[numSameGIDs]
        numSameGIDs += 1
    end
    numSameGIDs -= 1
    numSameIDs(data, numSameGIDs)
    
    #TODO refacter into a seperatate method, duplicated in the expert constructor
    permuteToLIDs = permuteToLIDs(data)
    permuteFromLIDs = permuteFromLIDs(data)
    remoteLIDs = remoteLIDs(data)
    for tgtLID = (numSameGIDs+1):numTgtLIDs
        const curTargetGID = targetGIDs[tgtLID]
        const srcLID = lid(source, curTargetLID)
        if srcLID != 0
            push!(permuteToLIDs, tgtLID)
            push!(permuteFromLIDs, srcLID)
        else
            push!(remoteGIDs, curTargetGID)
            push!(remoteLIDs, tgtLID)
        end
    end
    
    if length(remoteLID) != 0 && !isDistributed(source)
        isLocallyComplete(data, false)
        
        warn("Target has remote LIDs but source is not distributed globally.  " *
            "Importing a submap of the target map")
    end
    
    remoteGIDs
end

#TODO implement setupExport
function setupExpert(impor::Import{GID, PID, LID}, remoteGIDs::Array{GID}, useRemotePIDs::Bool, userRemotePIDs::Array{PID}) where {GID <: Integer, PID <: Integer, LID <: Integer}
    
    const source = sourceMap(impor)
    
    # Sanity Checks
    if !useRemotePIDs && length(userRemotePIDs) > 0
        throw(InvalidArgumentError("remotePIDs are non-empty, but " *
                "their use has not been requested"))
    end
    if length(userRemotePIDs) > 0 && length(userRemotePIDs) != length(remoteGIDs)
        throw(InvalidArgumentError("remotePIDs must either be of size zero " *
                "or match the size of remoteGIDs."))
    end
    
    
    missingGID = 0
    
    if !useRemotePIDs
        newRemotePIDs = Array{PID, 1}(length(remoteGIDs))
        if impor.debug
            print("$(myPid(comm(source))): setupExport(Import): about to call " * 
                "getRemoteIndexList on sourceMap")
        end
        (remoteProcIds, remoteLIDs) = remoteIDList(source, remoteGIDs)
        for e in remoteLIDs
            if e == 0
                missingGID += 1
            end
        end
    else
        remoteProcIds = userRemotePIDs
    end
    
    #line 688
    
    if missingGID != 0
        isLocallyComplete(data, false)
        
        warn("Source map was un-able to figure out which process owns one " *
            "or more of the GIDs in the list of remote GIDs.  This probably " *
            "means that there is at least one GID owned by some process in " *
            "the target map which is not owned by any process in the source " *
            "Map.  (That is, the source and target maps do not contain the " *
            "same set of GIDs globally")
        
        #ignore remote GIDs that aren't owned by any process in the source Map
        numInvalidRemote = missingGID
        totalNumRemote = length(remoteGIDs)
        if numInvalidRemote == totalNumRemote
            #if all remotes are invalid, can delete them all
            empty!(remoteProcIDs)
            empty!(remoteGIDs)
            empty!(remoteLIDs(data))
        else
            numValidRemote = 1
            
            remoteLIDs = remoteLIDs(data)
            
            for r = 1:totalNumRemote
                if remoteProcIds[r] != 0
                    remoteProcIds[numValidRemote] = remoteProcIDs[r]
                    remoteGIDs[numValidRemote] = remoteGIDs[r]
                    remoteLIDs[numValidRemote] = remoteLIDs[r]
                    numValidRemote += 1
                end
            end
            numValidRemote -= 1
            
            if numValidRemote != totalNumRemote - numInvalidRemote
                throw(InvalidStateException("numValidRemote = $numValidRemote " *
                        "!= totalNumRemote - numInvalidRemote " *
                        "= $(totalNumRemote - numInvalidRemote)"))
            end
            
            resize!(remoteProcIds, numValidRemote)
            resize!(remoteGIDs, numValidRemote)
            resize!(remoteLIDs, numValidRemote)
        end
    end
    
    order = sortperm(remoteProcIDs)
    permute!(remoteProcIDs, order)
    permute!(remoteGIDs, order)
    permute!(remoteLIDs, order)
    
    (exportGIDs, exportPIDs) = createFromRecvs(distributor(data), remoteGIDs, remoteProcs)
    
    julia_petra.exportPIDs(data, exportPIDs)
    
    numExportIDs = length(exportGIDs)
    
    if numExportIDs > 0
        resize!(exportLIDs(data), numExportIDs)
        exportLIDs = exportLIDs(data)
        for k in 1:numExportIDs
            exportLIDs[k] = lid(source, exportGIDs[k])
        end
    end
    
    if impor.debug
        print("$(myPid(comm(source))): setupExport: done\n")
    end
end

## Getters ##

"""
Get the source map for the given ImportExportData
"""
function sourceMap(impor::Import{GID, PID, LID})::BlockMap{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.data.source
end

"""
Get the target map for the given ImportExportData
"""
function targetMap(impor::Import{GID, PID, LID})::BlockMap{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.data.target
end

"""
List of elements in the target map that are permuted.
"""
function permuteToLIDs(impor::Import{GID, PID, LID})::Array{LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.data.permuteToLIDs
end

"""
List of elements in the source map that are permuted.
"""
function permuteFromLIDs(impor::Import{GID, PID, LID})::Array{LID} where GID <: Integer where PID <:Integer where LID <: Integer
    impor.data.permuteFromLIDs
end

"""
List of elements in the target map that are coming from other processors
"""
function remoteLIDs(impor::Import{GID, PID, LID})::Array{LID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.data.remoteLIDs
end

"""
List of elements that will be sent to other processors
"""
function exportLIDs(impor::Import{GID, PID, LID})::Array{LID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.data.exportLIDs
end

"""
List of processors to which elements will be sent `exportLID[i]` will be sent to processor `exportPIDs[i]`
"""
function exportPIDs(impor::Import{GID, PID, LID})::Array{PID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.data.exportPIDs
end

"""
Returns the number of elements that are identical between the source and target maps, up to the first different ID
"""
function numSameIDs(impor::Import{GID, PID, LID})::LID where GID <: Integer where PID <: Integer where LID <: Integer
    impor.data.numSameIDs
end


"""
Returns the distributor being used
"""
function distributor(impor::Import{GID, PID, LID})::Distributor{GID, PID, LID} where GID <: Integer where PID <: Integer where LID <: Integer
    impor.data.distributor
end

"""
Returns whether the import or export is locally complete
"""
function isLocallyComplete(impor::Import{GID, PID, LID})::Bool where GID <: Integer where PID <: Integer where LID <: Integer
    impor.data.isLocallyComplete
end