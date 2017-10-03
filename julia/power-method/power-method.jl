using julia_petra

#program settings
const useMPI = true

const commData = Float64
const commGID = UInt64
const commPID = UInt8
const commLID = UInt32


function powerMethod(A::RowMatrix{Data, GID, PID, LID}, niters::Integer,
        tolerance, verbose::Bool)::Tuple{Data, Bool} where {Data, GID, PID, LID}
    powerMethod(A, niters, Data(tolerance), verbose)
end

"""
Returns a tuple of the computed λ and if the λ is within tolerance
"""
function powerMethod(A::RowMatrix{Data, GID, PID, LID}, niters::Integer,
        tolerance::Data, verbose::Bool)::Tuple{Data, Bool} where {Data, GID, PID, LID}
    q = MultiVector{Data, GID, PID, LID}(getRowMap(A), 1)
    z = MultiVector{Data, GID, PID, LID}(getRowMap(A), 1)
    resid = MultiVector{Data, GID, PID, LID}(getRowMap(A), 1)

    #skipping flop counting

    randn!(z.data)

    valid = false

    for iter = 1:niters
        normz = norm2(z)[1]
        q = scale(z, 1.0/normz)
        apply!(z, A, q, Data(1), Data(0))
        λ = dot(q, z)[1]
        if iter%100 != 0 || iter+1 == niters
            #TODO improve - currently works, but is a little bit of a hack around MultiVector's lack of math operators
            @. resid.data = z.data - λ*q.data
            residual = norm2(resid)[1]

            if residual < tolerance
                return (λ, true)
            end
        end
    end
    (λ, false)
end

macro log(values...)
    esc(:(verbose && println($(values...))))
end

function main(comm::Comm{GID, PID, LID}, numGlobalElements, verbose, Data::Type) where{GID, PID, LID}

    const pid = myPid(comm)
    const nProc = numProc(comm)

    map = BlockMap(numGlobalElements, comm)

    numMyElements = julia_petra.numMyElements(map)
    myGlobalElements = julia_petra.myGlobalElements(map)
    numNz = Array{GID, 1}(numMyElements)
    for i = 1:numMyElements
        if myGlobalElements[i] == 1 || myGlobalElements[i] == numGlobalElements
            numNz[i] = 2
        else
            numNz[i] = 3
        end
    end

    const A = CSRMatrix{Data}(map, numNz, STATIC_PROFILE)

    values = Data[-1, -1]
    indices = Array{GID, 1}(2)#TODO does this need to be GID?
    two = Data(2)

    for i = 1:numMyElements
        if myGlobalElements[i] == 1
            indices = LID[1]
        elseif myGlobalElements[i] == numGlobalElements
            indices = LID[numGlobalElements-2]
        else
            indices = LID[myGlobalElements[i]-1, myGlobalElements[i]+1]
        end

        insertGlobalValues(A, myGlobalElements[i], indices, values)
        insertGlobalValues(A, myGlobalElements[i], LID[myGlobalElements[i]], Data[two])
    end

    fillComplete(A)

    const niters = numGlobalElements*10
    const tolerance = 1.0e-2

    #compile all the nessacery methods before timing
    powerMethod(A, niters, tolerance, verbose)

    tic()
    λ, success = powerMethod(A, niters, tolerance, verbose)
    #TODO look into storing FLOPS
    elapsedTime = toc()


    @log "λ = $λ; is within tolerance? $success"
    @log "\n\n total time for first solve = $elapsedTime sec\n\n"
    @log "increasing magnitude of first diagonal term, solving again\n\n"

    if myGID(map, 1)
        rowInds, rowVals = getLocalRowView(A, 0)
        for i = 1:getNumEntriesInGlobalRow(A, 0)
            if rowInds[i] == 0
                #using a view, so values update the original
                rowVals[i] *= 10
            end
        end
    end

    tic()
    λ, success = powerMethod(A, niters, tolerance, verbose)
    elapsedTime = toc()

    @log ""

    @log "λ = $λ; is within tolerance? $success"
    @log "\n\n total time for second solve = $elapsedTime sec\n\n"
end



if useMPI
    comm = MPIComm(commGID, commPID, commLID)
else
    comm = SerialComm{commGID, commPID, commLID}()
end

const pid = myPid(comm)
const nProc = numProc(comm)
const verbose = pid == 1

@log comm

if length(ARGS) != 1
    @log "Usage: $(ARGS[1]) number_of_equations"
#    @log "Exactly 2 arguments are required"
    exit(1)
end

numGlobalElements = parse(commLID, ARGS[1])
if numGlobalElements < nProc
    @log "numGlobalBlocks = $numGlobalElements cannot be < number of processors = $nProc"
    exit(1)
end

main(comm, numGlobalElements, verbose, commData)

print("\n\n\nSecond pass\n\n\n")

main(comm, numGlobalElements, verbose, commData)
