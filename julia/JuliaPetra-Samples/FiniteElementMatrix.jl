
importall JuliaPetra

struct FE2dMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: RowMatrix{Data, GID, PID, LID}
    graph::FE2dGraph{GID, PID, LID}

    problemWidth::GID
    problemHeight::GID
    rangeMap::BlockMap{GID, PID, LID}
    domainMap::BlockMap{GID, PID, LID}
    colMap::BlockMap{GID, PID, LID}
    importer::Nullable{Import{GID, PID, LID}}

    linearIndices::Bool

    function FE2dMatrix{Data, GID, PID, LID}(problemWidth::GID, problemHeight::GID,
            rangeMap::BlockMap{GID, PID, LID}, domainMap::BlockMap{GID, PID, LID},
            colMap::BlockMap{GID, PID, LID}, importer::Nullable{Import{GID, PID, LID}},
            linearIndices::Bool
            ) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}

        new{Data, GID, PID, LID}(FE2dGraph(problemWidth, problemHeight, rangeMap, domainMap, colMap, importer),
            problemWidth, problemHeight, rangeMap, domainMap, colMap, importer, linearIndices)
    end
end

function FE2dMatrix{Data}(problemWidth::Integer, problemHeight::Integer, rangeMap::BlockMap{GID, PID, LID},
            domainMap::BlockMap{GID, PID, LID} = rangeMap) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    FE2dMatrix{Data}(GID(problemWidth), GID(problemHeight), rangeMap::BlockMap{GID, PID, LID}, domainMap::BLockMap{GID, PID, LID})
end

function FE2dMatrix{Data}(problemWidth::GID, problemHeight::GID, rangeMap::BlockMap{GID, PID, LID},
            domainMap::BlockMap{GID, PID, LID} = rangeMap) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}

    if !uniqueGIDs(rangeMap)
        error("Can only build a finite element matrix for a one to one map")
    end

    colMap, linearIndices = buildFEColMap(problemWidth, problemHeight, rangeMap, domainMap)

    #make the resulting matrix
    if colMap != domainMap && !sameAs(colMap, domainMap)
        importer = Nullable(Import(domainMap, colMap))
    else
        importer = Nullable{Import{GID, PID, LID}}()
    end
    FE2dMatrix{Data, GID, PID, LID}(problemWidth, problemHeight, rangeMap, domainMap, colMap, importer, linearIndices)
end

function buildFEColMap(problemWidth::GID, problemHeight::GID, rangeMap::BlockMap{GID, PID, LID},
            domainMap::BlockMap{GID, PID, LID}) where {GID <: Integer, PID <: Integer, LID <: Integer}

    @assert minLID(rangeMap) == 1
    @assert minLID(domainMap) == 1

    comm = getComm(rangeMap)

    if numProc(comm) == 1
        return domainMap, linearMap(domainMap)
    end
    if linearMap(rangeMap)
        # start with main diagonal
        baseRange = minMyGID(rangeMap) : maxMyGID(rangeMap)

        range = baseRange
        # add tridiagonal, unless at edge of grid
        if minMyGID(rangeMap)%problemWidth != 1 && minMyGID(rangeMap) != minAllGID(rangeMap)
            range = minMyGID(rangeMap)-1 : maxMyGID(rangeMap)
        end
        if maxMyGID(rangeMap)%problemWidth != 0 && maxMyGID(rangeMap) != maxAllGID(rangeMap)
            range = first(range) : maxMyGID(rangeMap)+1
        end

        # add other diagonals
        # work around needed negatives if GID is unsigned
        if first(baseRange) <= problemWidth
            if last(baseRange) <= problemWidth
                rightDiag = 1 : 0
            else
                rightDiag = 1 : last(baseRange) - problemWidth
            end
        else
            rightDiag = baseRange - problemWidth
        end
        leftDiag = baseRange+problemWidth

        globalRange = minAllGID(rangeMap) : maxAllGID(rangeMap)


        if last(rightDiag)+1 >= first(range)
            # then last(range) >= first(leftDiag)
            myColumns = intersect(first(rightDiag) : last(leftDiag),
                                  globalRange)
            return BlockMap(numGlobalElements(domainMap), myColumns, comm), true
        else
            #then last(range) < first(leftDiag)
            rightDiag = intersect(rightDiag, globalRange)
            leftDiag = intersect(leftDiag, globalRange)
            myColumns = vcat(rightDiag, range, leftDiag)
        end
    else
        const localNumRows = getLocalNumRows(graph)
        numLocalColGIDs = 0
        gidIsLocal = falses(localNumRows)
        remoteGIDSet = Set()

        for localRow = 1:localNumRows
            globalRow = gid(rowMap, localRow)
            for neighbor in getNeighbors(globalRow)
                localNeighbor = lid(domainMap, neighbor)
                if localNeighbor != 0
                    @inbounds if !gidIsLocal[lid]
                        gidIsLocal[lid] = true
                        numLocalColGIDs += 1
                    end
                else
                    push!(remoteGIDSet, neighbor)
                end
            end
        end
        myColumns = Vector{GID}(numLocalColGIDs+numRemoteColGIDs)
        localColGIDs  = view(myColumns, 1:numLocalColGIDs)
        remoteColGIDs = view(myColumns, numLocalColGIDs+(1:numRemoteColGIDs))

        remoteColGIDs[:] = collect(remoteGIDSet)

        remotePIDs = Array{PID, 1}(numRemoteColGIDs)

        remotePIDs = remoteIDList(domainMap, remoteColGIDs)[1]
        order = sortperm(remotePIDs)
        permute!(remotePIDs, order)
        permute!(remoteColGIDs, order)

        # PID's of 0 will be sorted to the front
        #if any(remotePIDs .== 0)
        if remotePIDs[1] == 0
            error("Some column indices are not in the domain Map")
        end

        #line 333
        numDomainElts = numMyElements(domainMap)
        if numLocalColGIDs == numDomainElts
            if linearMap(domainMap)
                localColGIDs[1:numLocalColGIDs] = minMyGID(domainMap) : max
            else
                domElts = myGlobalElements(domainMap)
                localColGIDs[1:length(domElts)] = domElts
            end
        else
            numLocalCount = 0
            if linearMap(domainMap) #I think isContiguous() <=> linearMap()
                for i = 1:numDomainElts
                    if gidIsLocal[i]
                        numLocalCount += 1
                        localColGIDs[numLocalCount] = curColMapGID
                    end
                end
            else
                domainElts = myGlobalElement(domainMap)
                for i = 1:numDomainElts
                    if gidIsLocal[domainElts[i]]
                        numLocalCount += 1
                        localColGIDs[numLocalCount] = domainElts[i]
                    end
                end
            end
            if numLocalCount != numLocalColGIDs
                error("$(myPid(getComm(graph))): numLocalCount = $numLocalCount "
                    * "!= numLocalColGIDs = $numLocalColGIDs.  "
                    * "This should not happen.")
            end
        end
    end
    return BlockMap(numGlobalElements(domainMap), myColumns, comm), false
end

function getNeighbors(problemWidth::GID, problemHeight::GID, i::GID) where GID <: Integer
    x = ((i-1) % problemWidth) + 1
    y = fld((i-1), problemWidth)+1
    arr = GID[]
    if x > 1
        push!(arr, x-1 + y*problemWidth)
    end
    if x < problemWidth
        push!(arr, x+1 + y*problemWidth)
    end
    if y > 1
        push!(arr, x + (y-1)*problemWidth)
    end
    if y < problemHeight
        push!(arr, x + (y-1)*problemWidth)
    end
    arr
end


##### general public API

getProblemWidth(mat::FE2dMatrix) = mat.problemWidth
getProblemHeight(mat::FE2dMatrix) = mat.problemHeight

##### Row graph/row matrix methods

getRowMap(mat::FE2dMatrix) = mat.rangeMap
getColMap(mat::FE2dMatrix) = mat.colMap
hasColMap(mat::FE2dMatrix) = true
getImporter(mat::FE2dMatrix) = mat.importer

getGraph(mat::FE2dMatrix) = mat.graph

leftScale!(::FE2dMatrix, ::AbstractArray) = error("Finite Element matrices can't be modified")
rightScale!(::FE2dMatrix, ::AbstractArray) = error("Finite Element matrices can't be modified")

#the rows are never stored, so need to copy the row
function getGlobalRowView(matrix::FE2dMatrix{Data, GID}, row::GID) where {Data, GID}
    getGlobalRowCopy(matrix, row)
end

function getGlobalRowCopy(matrix::FE2dMatrix{Data, GID}, row::GID) where {Data, GID}
    globalWidth = matrix.problemWidth
    globalHeight = matrix.problemHeight

    globalX = mod(row-1, globalWidth) + 1
    globalY = fld(row-1, globalWidth) + 1

    inds = GID[]
    vals = Data[]

    if globalY > 1
        push!(inds, row-globalWidth)
        push!(vals, Data(-1))
    end
    if globalX > 1
        push!(inds, row-1)
        push!(vals, Data(-1))
    end
    push!(inds, row)
    push!(vals, Data(4))
    if globalX < globalWidth
        push!(inds, row+1)
        push!(vals, Data(-1))
    end
    if globalY < globalHeight
        push!(inds, row+globalWidth)
        push!(vals, Data(-1))
    end
    (inds, vals)
end

#the rows are never stored, so need to copy the row
function getLocalRowView(matrix::FE2dMatrix{Data, GID}, row::GID) where {Data, GID}
    getGlobalRowCopy(matrix, row)
end

function getLocalRowCopy(mat::FE2dMatrix{Data, GID, PID, LID}, row::LID) where {Data, GID, PID, LID}
    rowMap = getRowMap(mat)
    inds, vals = getGlobalRowCopy(mat, gid(rowMap, row))
    colMap = getColMap(mat)
    map(x->lid(colMap, x), inds), vals
end

function getLocalDiagCopy(mat::FE2dMatrix{Data}) where Data
    @assert isa(mat, RowMatrix)
    fill(Data(4), getLocalNumDiags(mat))
end

##### Operator methods

getRangeMap(mat::FE2dMatrix) = mat.rangeMap
getDomainMap(mat::FE2dMatrix) = mat.domainMap


#=
TypeStability.@stable_function [(MultiVector{D, G, P, L}, FE2dMatrix{D, G, P, L},
                        MultiVector{D, G, P, L}, TransposeMode, D, D)
                    for (D, G, P, L) in Base.Iterators.product(
                        [Float64, Complex64], #Data
                        [UInt64, Int64, UInt32], #GID
                        [UInt8, Int8, UInt32], #PID
                        [UInt32, Int32]) #LID
] nothing begin
function localApply(Y::MultiVector{Data, GID, PID, LID},
        A::FE2dMatrix{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID},
        mode::TransposeMode, alpha::Data, beta::Data) where {Data, GID, PID, LID}

    const rawY = Y.data
    const rawX = X.data

    const problemWidth = A.problemWidth
    const problemHeight = A.problemHeight
    @assert problemWidth > 1
    @assert problemHeight > 1

    const numRows = getLocalNumRows(A)

    if A.linearIndices
        const firstGlobalRow = gid(getRowMap(A), 1)
        const firstX = mod(firstGlobalRow-1, problemWidth)+1
        const firstY = fld(firstGlobalRow-1, problemWidth)+1
        const firstCol = lid(getColMap(A), firstX+(firstY-1)*problemWidth)

        for vect = LID(1):numVectors(Y)

            x::GID = firstX
            y::GID = firstY
            row::LID = 1
            col::LID = firstCol


            while row <= numRows
                sum::Data = Data(4)*rawX[col, vect];
                if x > 1
                    sum -= rawX[col-1, vect]
                end
                if x < problemWidth
                    sum -= rawX[col+1, vect]
                end
                if y > 1
                    sum -= rawX[col-problemWidth, vect]
                end
                if y < problemHeight
                    sum -= rawX[col+problemWidth, vect]
                end

                sum = applyConjugation(mode, sum*alpha)
                @inbounds rawY[row, vect] = sum + beta*rawY[row, vect]

                row += 1
                col += 1
                x += 1
                if x > problemWidth
                    x = GID(1)
                    y += 1
                end
            end
        end
    else
        for vect = LID(1):numVectors(Y)
            for row = LID(1):numRows
                globalRow = gid(getRowMap(A), row)
                col = lid(getColMap(A), globalRow)
                x = mod(globalRow-1, problemWidth)+1
                y = fld(globalRow-1, problemWidth)+1

                sum::Data = Data(4)*rawX[col]
                if x > 1
                    sum -= rawX[lid(getColMap(A), globalRow-1)]
                end
                if x < problemWidth
                    sum -= rawX[lid(getColMap(A), globalRow+1)]
                end
                if y > 1
                    sum -= rawX[lid(getColMap(A), globalRow-problemWidth)]
                end
                if y < problemHeight
                    sum -= rawX[lid(getColMap(A), globalRow+problemWidth)]
                end

                sum = applyConjugation(mode, sum*alpha)
                @inbounds rawY[row, vect] = sum + beta*rawY[row, vect]
            end
        end
    end
end
end

# =#
