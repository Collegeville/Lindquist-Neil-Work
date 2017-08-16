export CSRMatrix

mutable struct CSRMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistRowMatrix{Data, GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::Nullable{BlockMap{GID, PID, LID}}
    
    #TODO look into how nessacery these are
    importMV::Nullable{MultiVector{Data, GID, PID, LID}}
    exportMV::Nullable{MultiVector{Data, GID, PID, LID}}

    myGraph::CRSGraph{GID, PID, LID}
    
    localMatrix::LocalCSRMatrix{Data, LID}

    values2D::Array{Array{Data, 1}, 1}

    #pull storageStatus and fillComplete  from graph

    #Dict keys are row indices
    #first element of each tuple is a column index
    #second element of each tuple is the matching entry
    nonlocals::Dict{GID, Tuple{Array{Data, 1}, Array{GID, 1}}}

    plist::Dict{Symbol}

    function CSRMatrix{Data, GID, PID, LID}(rowMap::BlockMap{GID, PID, LID}, colMap::Nullable{BlockMap{GID, PID, LID}}, myGraph::CRSGraph{GID, PID, LID}, localMatrix::LocalCSRMatrix{Data, LID}, values2D::Array{Array{Data, 1}, 1}, plist::Dict{Symbol}) where {Data, GID, PID, LID}
        new(rowMap,
            colMap,
            Nullable{MultiVector{Data, GID, PID, LID}}(),
            Nullable{MultiVector{Data, GID, PID, LID}}(),
            myGraph,
            values2D,
            Dict{GID, Tuple{Array{Data, 1}, Array{GID, 1}}}(),
            plist)
    end
        
end

#### Constructors ####
#TODO document Constructors

function CSRMatrix{Data}(rowMap::BlockMap{GID, PID, LID},
        maxNumEntriesPerRow::Union{Integer, Array{Integer, 1}}, 
        pftype::ProfileType, plist::Dict{Symbol}) where {Data, GID, PID, LID}
    CSRMatrix{Data}(rowMap, Nullable{BlockMap{GID, PID, LID}}(),
        maxNumEntriesPerRow, pftype, plist)
end
function CSRMatrix{Data}(rowMap::BlockMap{GID, PID, LID},
        colMap::BlockMap{GID, PID, LID},
        maxNumEntriesPerRow::Union{Integer, Array{Integer, 1}}, 
        pftype::ProfileType, plist::Dict{Symbol}) where {Data, GID, PID, LID}
    CSRMatrix{Data}(rowMap, Nullable(colMap), maxNumEntriesPerRow,
        pftype, plist)
end
function CSRMatrix{Data}(rowMap::BlockMap{GID, PID, LID},
        colMap::Nullable{BlockMap{GID, PID, LID}},
        maxNumEntriesPerRow::Union{Integer, Array{Integer, 1}}, 
        pftype::ProfileType, plist::Dict{Symbol}) where {Data, GID, PID, LID}
    graph = CRSGraph(rowMap, maxNumEntriesPerRow, pftype, plist)
    
    matrix = CSRMatrix{Data, GID, PID, LID}(rowMap, colMap, 
        graph, LocalCSRMatrix{Data, LID}(), [], plist)
    
    resumeFill(matrix, plist)
    
    matrix
end

function CSRMatrix{Data}(graph::CRSGraph{GID, PID, LID},plist::Dict{Symbol}
        ) where {Data, GID, PID, LID}
    numCols = numMyElements(getColMap(graph))
    localGraph = getLocalGraph(graph)
    val = Array{Data, 1}(length(localGraph.entries))
    localMatrix = LocalCSRMatrix(numCols, val, localGraph)
    
    CSRMatrix(graph.rowMap, graph.colMap, graph, localMatrix, [], plist)
end

function CSRMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        rowOffsets::Array{LID, 1}, colIndices::Array{LID, 1}, values::Array{Data, 1},
        plist::Dict{Symbol}) where {Data, GID, PID, LID}
    
    #check user's input.  Might throw on only some processes, causing deadlock
    if length(values) != length(colIndices)
        throw(InvalidArgumentError("values and columnIndices must "
                * "have the same length"))
    end
    #TODO add debug from line 328?
    #if so, use a minAll on the previous error when debuging to propagate the error
    
    graph = CRSGraph(rowMap, colMap, rowOffsets, columnIndices, plist)
    localGraph = getLocalGraph(graph)
    
    numCols = numMyElements(getColMap(graph))
    localMatrix = LocalCSRMatrix(numCols, values, localGraph)
    
    CSRMatrix(rowMap, Nullable(colMap), graph, localMatrix, [], plist)
end


function CSRMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localMatrix::LocalCSRMatrix{Data, LID}, plist::Dict{Symbol}
        ) where {Data, GID, PID, LID}
    
    graph = CRSGraph(rowMap, colMap, localMatrix.graph, plist)
    
    matrix = CSRMatrix(rowMap, colMap, graph, localMatrix, [], plist)
    
    computeGlobalConstants(matrix)
    matrix
end

function CSRMatrix(rowMap::BlockMap{GID, PID, LID}, colMap::BlockMap{GID, PID, LID},
        localMatrix::AbstractArray{Data, 2}, plist::Dict{Symbol}
        ) where {Data, GID, PID, LID}
    linearIndices = find(x -> x!=0, localMatrix)
    rowIndicesIter, colIndicesIter, valuesIter = zip(
        sort!(collect(zip(ind2sub(size(localMatrix), linearIndices)...,
                          localMatrix[linearIndices])))...)
    rowIndices = collect(rowIndicesIter)
    rowOffsets = Array{LID, 1}(size(localMatrix, 1)+1)
    row = 1
    j = 1
    for i in 1:length(rowIndices)
        if rowIndices[i] > row
            row += 1
            rowOffsets[row] = i
        end
    end
    rowOffsets[length(rowOffsets)] = length(rowIndices)+1
    
    CSRMatrix(rowMap, colMap, rowOffsets, 
        collect(colIndicesIter), collect(valuesIter), plist)
end


#### Internal methods ####
#does nothing, exists only to be a parallel to CRSGraph
computeGlobalConstants(matrix::CSRMatrix) = nothing

#Tpetra's only clears forbNorm, exists only to be a parallel to CRSGraph
clearGlobalConstants(matrix::CSRMatrix) = nothing

function globalAssemble(matrix::CSRMatrix)
    comm = julia_petra.comm(matrix)
    if !isFillActive(matrix)
        throw(InvalidStateError("Fill must be active to call this method"))
    end
    
    numNonLocalRows = length(matrix.nonlocals)
    nooneHasNonLocalRows = maxAll(comm, myNumNonLocalRows == 0)
    if nooneHasNonLocalRows
        #no process has nonlocal rows, so nothing to do
        return
    end
    
    #nonlocalRowMap = BlockMap{GID, PID, LID}()
    
    numEntPerNonlocalRow = Array{Integer, 1}(numNonLocalRows)
    myNonlocalGlobalRows = Array{GID, 1}(numNonLocalRows)
    
    curPos = 1
    for (key, val) in matrix.nonlocal
        myNonlocalGlobalRows[curPos] = key
        
        values = val[1]
        globalColumns = val[2]
        
        order = sortperm(globalColumns)
        permute!(globalColumns, order)
        permute!(values, order)
        
        
        curPos += 1
    end
    
    #merge2
    if length(globalColumns) > 0
        setIndex = 1
        valuesSum = 0
        currentIndex = globalColumns[1]
        for i = 1:length(globalColumns)
            if currentIndex != globalColumns[i]
                values[setIndex] = valuesSum
                globalColumns[setIndex] = currentIndex
                setIndex += 1
                valuesSum = 0
            end
            valuesSum += values[i]
            i += 1
        end
        values[setIndex] = valuesSum
        globalColumns[setIndex] = currentIndex
        
        resize!(values, setIndex)
        resize!(globalColumns, setIndex)
    end
    
    numEntPerNonloalRow[curPose] = length(globalColumns)
    
    
    #don't need to worry about finding the indexBase
    nonlocalRowMap = BlockMap(0, myNonlocalGlobalRows, comm)
    
    nonlocalMatrix = CSRMatrix(nonlocalRowMap, numEntPernonlocalRow, STATIC_PROFILE)
    
    curPos = 1
    for (key, val) in matrix.nonlocals
        globalRow = key
        
        vals = val[1]
        globalCols = val[2]
        
        insertGlobalValues(nonlocalMatrix, globalRow, globalCols, vals)
    end
    
    origRowMap = rowMap(matrix)
    origRowMapIsOneToOne = isOneToOne(origRowMap)
    
    
    if origRowMapIsOneToOne
        exportToOrig = Export(nonlocalRowMap, origRowMap)
        isLocallyComplete = isLocallyComplete(exportToOrig)
        doExport(nonlocalMatrix, matrix, exportToOrig, ADD)
    else
        oneToOneRowMap = createOneToOne(origRowMap)
        exportToOneToOne = Export(nonlocalRowMap, oneToOneRowMap)
        
        isLocallyComplete = isLocallyComplete(exportToOneToOne)
    
        oneToOneMatrix = CSRMatrix{Data}(oneToOneRowMap, 0)
        
        doExport(nonlocalMatrix, onToOneMatrix, exportToOneToOne, ADD)
        
        #nonlocalMatrix = null
        
        importToOrig = Import(oneToOneRowMap, origRowMap)
        doImport(oneToOneMatrix, matrix, importToOrig, ADD)
    end
    
    empty!(matrix.nonlocals)
    
    globallyComplete = minAll(comm, isLocallyComplete)
    if !globallyComplete
        throw(InvalidArgumentError("On at least one process, insertGlobalValues "
                * "was called with a global row index which is not in the matrix's "
                * "row map on any process in its communicator."))
    end
end

    
function fillLocalGraphAndMatrix(matrix::CSRMatrix{Data, GID, PID, LID},
        plist::Dict{Symbol}) where {Data, GID, PID, LID}
    localNumRows = getNodeNumRows(matrix)
    
    myGraph = matrix.myGraph
    localMatrix = matrix.localMatrix
    
    #TODO figure out if correct
    matrix.localMatrix.indices = myGraph.localIndices1D
    
    #most of the debug sections were taken out, since they're for debuging Petra itself, and julia doesn't have a compiler option to enable macros
    if getProfileType(matrix) == DYNAMIC_PROFILE
        numRowEntries = myGraph.numRowEntries
        
        ptrs = Array{LID, 1}(localNumRows+1)
        localTotalNumEntries = computeOffsetsFromCounts(ptrs, numRowEntries)
        
        inds = Array{LID, 1}(localTotalNumEntries)
        vals = Array{Data, 1}(localTotalNumEntries)
        
        localIndices2D = myGraph.localIndices2D
        for row = 1:localNumRows
            numEnt = numRowEnt[row]
            dest = range(ptrs[row], 1, numEnt)
            
            inds[dest] = localIndices2D[row][:]
            vals[dest] = matrix.values2D[row][:]
        end
    elseif getProfileType(matrix) == STATIC_PROFILE
        curRowOffsets = myGraph.rowOffsets
        
        if myGraph.STORAGE_STATUS == STORAGE_1D_UNPACKED
            #pack row offsets into ptrs
            
            localTotalNumEntries = 0
            
            packedRowOffsets = Array{LID, 1}(localNumRows + 1)
            numRowEnt = myGraph.numRowEntries
            localTotalNumEntries = computeOffsetsFromCounts(packedRowOffsets, numRowEnt)
            ptrs = packedRowOffsets
            
            inds = Array{LID, 1}(localTotalNumEntries)
            vals = Array{Data, 1}(localTotalNumEntries)
            
            #line 1234
            for row in 1:localNumRows
                srcPos = curRowOffsets[row]
                dstPos = ptrs[row]
                dstEnd = ptrs[row+1]
                dst = dstPos:dstEnd
                src = range(srcPos, 1, dstEnd-dstPos)
                
                vals[dst] = myGraph.localIndices1D[src]
                inds[dst] = localMatrix.values[src]
            end
        else
            #dont have to pack, just set pointers
            ptrs = myGraph.rowOffsets
            inds = myGraph.localIndices1D
            vals = localMatrix.values
        end
    end
    
    if get(plist, :optimizeStorage, true)
        empty!(myGraph.localIndices2D)
        empty!(myGraph.numRowEntries)
        
        empty!(matrix.values2D)
        
        myGraph.rowOffsets = ptrs
        myGraph.localIndices1D = inds
        
        myGraph.pftype = STATIC_PROFILE
        myGraph.storageStatus = STORAGE_1D_PACKED
    end
    
    myGraph.localGraph = LocalCRSGraph(inds, ptrs)
    matrix.localMatrix = LocalCSRMatrix(getNodeNumCols(matrix), vals, graph.localGraph)
end
        
function insertNonownedGlobalValues(matrix::CSRMatrix{Data, GID, PID, LID},
        globalRow::GID, indices::Array{GID, 1}, values::Array{Data, 1}
        ) where {Data, GID, PID, LID}
    
    curRow = matrix.nonlocals[globalRow]
    curRowVals = curRow[1]
    curRowInds = curRow[2]
    
    newCapacity = length(curRowInds) + length(indices)
    
    append!(curRowVals, values)
    append!(curRowInds, indices)
end
    

#### External methods ####
#TODO document external methods
    
function insertGlobalValues(matrix::CSRMatrix{Data, GID, PID, LID}, globalRow::Integer,
        indices::AbstractArray{LID, 1}, values::AbstractArray{Data, 1}
        ) where {Data, GID, PID, LID}
    myGraph = matrix.myGraph
    
    localRow = lid(getRowMap(matrix), globalRow)
    
    if localRow == 0
        insertNonownedGlobalValues(matrix, globalRow, indices, values)
    else
        numEntriesToInsert = length(indices)
        if hasColMap(matrix)
            colMap = getColMap(matrix)
            
            for k = 1:numEntriesToInsert
                if !myGID(colMap, indices[k])
                    throw(InvalidArgumentError("Attempted to insert entries into "
                            * "owned row $globalRow, at the following column indices "
                            * "$indices.  At least one of those indices ($(indices[k])"
                            * ") is not in the column map on this process"))
                end
            end
        end
        
        #TODO is this ok? from lines 1965-1969
        inds_view = indices
        vals_view = values
        
        rowInfo = getRowInfo(myGraph, localRow)
        curNumEntries = rowInfo.numEntries
        newNumEntries = curNumEntries = length(numEntriesToInsert)
        if newNumEntries > rowInfo.allocSize
            if(getProfileType(matrix) == STATIC_PROFILE 
                    && newNumEntries > rowInfo.allocSize)
                throw(InvalidArgumentException("number of new indices exceed "
                        * "statically allocated graph structure"))
            end
            
            rowInfo = updateGlobalAllocAndValues(myGraph, rowInfo, newNumEntries,
                            matrix.values2D[localRow])
        end

        if isGloballyIndexed(matrix)
            insertIndicesAndValues(myGraph, rowInfo, indices, getView(matrix, rowInfo),
                values, GLOBAL_INDICES, GLOBAL_INDICES)
        else
            insertIndicesAndValues(myGraph, rowInfo, indices, getView(matrix, rowInfo),
                values, GLOBAL_INDICES, LOCAL_INDICES)
        end
    end
end
            
    

function resumeFill(matrix::CSRMatrix, plist::Dict{Symbol})
    resumeFill(matrix.myGraph, plist)
    
    clearClobalConstants(matrix)
    #graph handles fillComplete variable
end

function fillComplete(matrix::CSRMatrix, plist::Dict{Symbol})
    #TODO figure out if the second arg should be getColMap(matrix)
    fillComplete(getRowMap(matrix), getRowMap(matrix), plist)
end

function fillComplete(matrix::CSRMatrix{Data, GID, PID, LID},
        domainMap::BlockMap{GID, PID, LID}, rangeMap::BlockMap{GID, PID, LID},
        plist::Dict{Symbol}) where {Data, GID, PID, LID}
    if isFillComplete(matrix)
        throw(InvalidStateError(
                "Matrix cannot be fill complete when fillComplete(...) is called"))
    end
    
    const myGraph = matrix.myGraph
    
    assertNoNonlocalInserts = get(plist, :noNonlocalChanges, false)
    #skipping sort ghosts stuff
    
    needGlobalAssemble = !assertNoNonlocalInserts && numProcs > 1
    if needGlobalAssemble
        globalAssemble(matrix)
    else
        if numProcs == 1 && length(nonlocals) != 0
            throw(InvalidStateError("Cannot have nonlocal entries on a serial run.  An invalid entry is present."))
        end
    end
    
    setDomainRanageMaps(myGraph, domainMap, rangeMap)
    if !hasColMap(myGraph)
        makeColMap(myGraph)
    end
    
    makeIndicesLocal(myGraph)
    
    sortAndMergeIndicesAndValues(myGraph, isSorted(myGraph), isMerged(myGraph))
    
    makeImportExport(myGraph)
    computeGlobalConstants(myGraph)
    graph.fillComplete = true
    checkInternalState(myGraph)
    fillLocalGraphAndMatrix(matrix, plist)
end
    


#TODO implement CSRMatrix methods
#TODO implement RowMatrix methods
#TODO implement DistObject methods
#TODO implement SrcDistObject methods