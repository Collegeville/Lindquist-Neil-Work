
importall JuliaPetra

struct FE2dGraph{GID <: Integer, PID <: Integer, LID <: Integer} <: RowGraph{GID, PID, LID}
    problemWidth::GID
    problemHeight::GID
    rangeMap::BlockMap{GID, PID, LID}
    domainMap::BlockMap{GID, PID, LID}
    colMap::BlockMap{GID, PID, LID}
    importer::Nullable{Import{GID, PID, LID}}
    exporter::Nullable{Export{GID, PID, LID}}
end

getRowMap(graph::FE2dGraph) = graph.rangeMap
getColMap(graph::FE2dGraph) = graph.colMap
hasColMap(graph::FE2dGraph) = true
getImporter(graph::FE2dGraph) = graph.importer
function getExporter(graph::FE2dGraph)
    if isnull(graph.exporter)
        graph.exporter = Nullable(Export(rangeMap, rangeMap))
    end
    get(graph.exporter)
end

getGlobalNumDiags(graph::FE2dGraph) = numAllElements(getRowMap(graph))
getLocalNumDiags(graph::FE2dGraph) = numMyElements(getRowMap(graph))

isFillComplete(::FE2dGraph) = true
isLowerTriangular(::FE2dGraph) = false
isUpperTriangular(::FE2dGraph) = false

isGloballyIndexed(::FE2dGraph) = true
isLocallyIndexed(::FE2dGraph) = true

function getNumEntriesInGlobalRow(graph::FE2dGraph{GID}, row::GID) where GID
    x = ((row-1) % graph.problemWidth) + 1
    y = fld(row-1, graph.problemWidth) + 1

    count = 3
    if x != 1 && x != graph.problemWidth
        count += 1
    end
    if y != 1 && y != graph.problemHeight
        count += 1
    end
    count
end

function getNumEntriesInLocalRow(graph::FE2dGraph{GID, PID, LID}, row::LID
        ) where {GID, PID, LID}
    getNumEntriesInGlobalRow(graph, gid(getRowMap(graph), row))
end

function getGlobalNumEntries(graph::FE2dGraph{GID}) where GID
    width = graph.problemWidth
    height = graph.problemHeight

    if width == 1
        if height == 1
            return 1
        else
            return 3*height-2
        end
    end
    if height == 1
        return 3*width-2
    end

    GID(5*(width-2)*(height-2) + 4*2*(height-2) + 4*2*(width-2) + 3*4)
end

function getLocalNumEntries(graph::FE2dGraph{GID, PID, LID}) where {GID, PID, LID}
    globalWidth = graph.problemWidth
    globalHeight = graph.problemHeight

    firstRow = minMyGIDs(getRowMap(graph))
    numRows = numMyElements(getRowMap(graph))

    firstX = mod(firstRow-1, globalWidth) + 1
    firstY = fld(firstRow-1, globalWidth) + 1

    lastX = mod(firstRow+numRows-1, globalWidth) + 1
    lastY = fld(firstRow+numRows-1, globalWidth) + 1

    firstRankYElements = firstY == 1 ? 1 : 2
    firstRank = (globalWidth-firstX+1)*(3+firstRankYElements) + (firstX == 1 ? -2 : -1)
    middleRanks = max(0, (lastY-firstY-1)*(5*globalWidth-2))
    lastRankYElements = lastY == globalHeight ? 2 : 1
    lastRank = (globalWidth-firstX+1)*(3+lastRankYElements) + (lastX == globalWidth ? -2 : -1)

    return firstRank + middleRanks+lastRank
end

function getGlobalMaxNumRowEntries(graph::FE2dGraph)
    globalWidth = graph.problemWidth
    globalHeight = graph.problemHeight

    if globalWidth >= 3
        max(3, globalHeight) + 2
    elseif globalWidth == 2
        max(3, globalHeight) + 1
    else #if globalWidth == 1
        max(3, globalHeight)
    end
end

function getLocalMaxNumRowEntries(graph::FE2dGraph)
    globalWidth = graph.problemWidth
    globalHeight = graph.problemHeight

    firstRow = minMyGIDs(getRowMap(graph))
    numRows = numMyElements(getRowMap(graph))

    if globalWidth >= 3
        if 1 + globalWidth <= firstRow <= globalHeight*globalWidth - globalWidth
            5
        elseif 1 + globalWidth <= firstRow + numRows <= globalHeight*globalWidth - globalWidth
            5
        elseif globalHeight >= 3
            if numRows == 1
                3
            else
                4
            end
        else #globalHeight <= 2
            min(2, globalHeight) + min(2, numRows)
        end
    elseif globalWidth == 2
        if 3 <= firstRow <= globalHeight - 2 || numRows > 6
            4
        elseif numRows > 2
            min(3, globalHeight) + 1
        else
            min(2, globalHeight) + 1
        end
    else #if globalWidth == 1
        if 1 < firstRow < globalHeight || numRows > 3
            3
        else
            min(numRows+1, globalHeight)
        end
    end
end

function getGlobalRowCopy(graph::FE2dGraph{GID}, row::GID) where GID
    globalWidth = graph.problemWidth
    globalHeight = graph.problemHeight

    globalX = mod(row-1, globalWidth) + 1
    globalY = fld(row-1, globalWidth) + 1

    inds = GID[]

    if globalY > 1
        push!(row-globalWidth, inds)
    end
    if globalX > 1
        push!(row-1, inds)
    end
    push!(row, inds)
    if globalX < globalWidth
        push!(row+1, inds)
    end
    if globalY < globalHeight
        push!(row+globalWidth, inds)
    end
    inds
end

function getLocalRowCopy(graph::FE2dGraph{GID, PID, LID}, row::LID) where {GID, PID, LID}
    rowMap = getRowMap(graph)
    map(x->lid(rowMap, x), getGlobalRowCopy(graph, gid(rowMap, row)))
end
