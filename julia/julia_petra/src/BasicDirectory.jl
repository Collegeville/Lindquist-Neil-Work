

"""
    BasicDirectory(map::BlockMap)
Creates a BasicDirectory, which implements the methods of Directory with
basic implmentations
"""
type BasicDirectory <: Directory
    map::BlockMap
    
    directoryMap::Nullable{BlockMap}
    
    procList::Array{Integer}
    procListLists::Array{Integer}
    
    entryOnMultipleProcs::Bool
    
    localIndexList::Array{Integer}
    allMinGIDs::Array{Integer}
    
    function BasicDirectory(map::BlockMap)
        if !(distributedGlobal(map))
            new(map, Nullable{BlockMap}(), [], [], false, [], [])
        elseif linearMap(map)
            commObj = comm(map)
            
            minMyGID = minMyGID(map)
            allMinGIDs = gatherAll(commObj, minMyGID)
            allMinGIDs = vcat(allMinGIDs, [1+maxAllGID(map)])
            
            new(map, Nullabl{BlockMap}, [], [], false, [], allMinGIDs)
        else
            generateContent(
                new(map, Nullable{BlockMap}, [], [], false, [], []),
                map)
        end
    end
end

function generateContent(dir::BasicDirectory, map::BlockMap)
    minAllGID = minAllGID(map)
    maxAllGID = maxAllGID(map)

    dirNumGlobalElements = maxAllGID - minAllGID + 1

    directoryMap = BlockMap(dirNumGlobalElements, minAllGID, commObj)
    dir.directoryMap = Nullable(directoryMap)

    dirNumMyElements = numMyElements(dir.directoryMap)

    if dirNumMyElements > 0
        dir.procList = Array{Integer}(dirNumMyElements)
        dir.localIndexList = Array{Integer}(dirNumMyElements)

        fill!(dir.procList, -1)
        fill!(dir.localIndexList, -1)
    else
        dir.procList = []
        dir.localIndexList = []
    end

    map_numMyElements = numMyElements(map)
    map_myGlobalElements = myGlobalElements(map)

    sendProcs = remoteIDList(directoryMap, map_numMyElements, myGlobalElements(map))

    distributor = createDistributor(commObj)
    numRecvs = createFromSends(distributor, map_numMyElements, sendProcs)

    exportElements = Array{Tuple{Integer, Integer, Integer}}(numMyElements)

    myPIDVal = myPID(commObj)
    for i = 1:numMyElements
        exportElements[i] = (map_myGlobalElements[i], myPIDVal, i)
    end

    importElements = resolve(distributor, exportElements)


    for i = 1:numRecvs
        currLID = lid(directoryMap, importElements[i][1])
        @assert currLID > 0 //internal error

        proc = importElements[i][2]
        if dir.procList[currLID] >= 0
            if dir.procList[currLID] != proc
                if dir.procLists == []
                    numProcLists = numMyElements(directoryMap)
                    procListLists = Array{Array{Integer}}(numProcLists)
                    fill!(procListLists, [])
                end

                l = procListLists[currLID]

                index = searchsortedfirst(l, procList[currLID])
                insert(l, index, procList[currLID])

                index = searchsortedfirst(l, proc)
                insert(l, index, proc)

                dir.procList[currLID] = dir.procListLists[curr_LID][1]
            end
        else
            dir.procList[currLID] = proc
        end

        dir.localIndexList[currLID] = importElements[i][3]
    end

    localVal = numProcLists
    globalVal = maxAll(commObj, localVal)
    dir.entryOnMultipleProcs = globalval > 0 ? true : false;

    dir
end


#    - Returns processor and local id infor for non-local map entries.  Returns a tuple
#        containing
#            1 - an Array of processors owning the global ID's in question
#            2 - an Array of local IDs of the global on the owning processor
function getDirectoryEntries(directory::BasicDirectory, map::BlockMap, globalEntries::Array{GID},
        high_rank_sharing_procs::Bool)::Tuple{Array{Integer}, Array{Integer}} where GID <: Integer
        #where PID <: Integer where LID <:Integer
    numEntries = length(globalEntries)
    procs = Array{Integer}(numEntries)
    localEntries = Array{Integer}(numEntries)
    
    if !distributedGlobal(map)
        myPIDVal = myPid(comm(map))
        
        for i = 1:numEntries
            lidVal = lid(map, globalEntries[i])
            
            if lidVal == 0
                procs[i] = 0
                warn("GID $(globalEntries[i]) is not part of this map")
            else
                procs[i] = myPIDVal
            end
            localEntries[i] = lidVal
        end
    elseif linearMap(map)
        minAllGIDVal = minAllGID(map)
        maxAllGIDVal = maxAllGID(map)
        
        numProcVal = numProc(comm(map))
        
        n_over_p = numGlobalElements(map)/numProcVal
        
        for i = 1:numEntries
            LID  = 0
            Proc = 0
            
            gid = globalEntries[i]
            if gid < minAllGID || gid > maxAllGID
                throw(InvalidArgumentError("GID out of valid range"))
            end
            #guess uniform distribution and start a little above it
            proc1 = min(gid/max(n_over_p, 1) + 2, numProcVal)
            found = false
            allMinGIDs_list = directory.allMinGIDs
            
            while proc1 >= 0 && proc1 < numProcVal
                if allMinGIDs_list[proc1] <= gid
                    if (gid < allMinGIDS_list[proc1 + 1])
                        found = true
                        break
                    else
                        proc1 += 1
                    end
                else
                    proc1 -= 1
                end
            end
            if found
                proc = proc1
                lid = gid - AllMinGIDs_list[proc]
            end
            
            procs[i] = proc
            localEntries[i] = lid
        end
    else #general case
        distributor = createDistributor(comm(map))
        
        dirProcs = remoteIDList(map, numEntries, globalEntries)
        
        numMissing = 0
        for i = 1:numEntries
            if dirProcs[i] == 0
                procs[i] = 0
                localEntries[i] = 0
                numMissing += 1
            end
        end
        
        (sendGIDs, sendPIDs) = createFromRecvs(distrbutor, globalEntries, dirProcs)
        numSends = length(sendGIDs)
        
        if numSends > 0
            exports = Array{Tuple{GID, Integer, Integer}}(numSends)
            for i = 1:numSends
                currGID = sendGIDs[i]
                exports[i][1] = currGID
                
                currLID = lid(map, currGID)
                @assert currLID > 0 #internal error
                if !high_rank_sharing_procs
                    exports[i][2] = procList[currLID]
                else
                    numProcLists = numMyElements(directory.directoryMap)
                    if numProcLists > 0
                        num = length(directory.procListLists[currLID])
                        if num > 1
                            exports[i][2] = directory.procListLists[currLID][num]
                        else
                            exports[i][2] = directory.procList[currLID]
                        end
                    else
                        exports[i][2] = directory.procList[currLID]
                    end
                end
                exports[i][3] = directory.localIndexList[currLID]
            end
        end
        
        numRecv = numEntries - numMissing
        imports = resolve(distributor, exports)
        
        offsets = sortperm(globalEntries)
        sortedGE = globalEntries[offsets]
        
        for i = 1:numRecv
            currLID = imports[i][1]
            j = searchsortedfirst(sortedGE, currLID)
            if j > 0
                procs[offsets[j]] = imports[i][2]
                localEntries[offsets[j]] = imports[i][3]
            end
        end
    end
    (procs, localEntries)
end

function gidsAllUniquelyOwned(directory::BasicDirectory)
    !directory.entryOnMultiplProcs
end
