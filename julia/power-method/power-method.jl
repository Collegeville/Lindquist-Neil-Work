using JuliaPetra

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
    q = DenseMultiVector{Data}(getRowMap(A), 1)
    z = DenseMultiVector{Data}(getRowMap(A), 1)
    resid = DenseMultiVector{Data}(getRowMap(A), 1)

    #skipping flop counting

    λ::Data = 0

    randn!(z.data)

    const ONE = one(Data)
    const ZERO = zero(Data)

    for iter = 1:niters
        normz = norm(z, 2)[1]

        @. q = z/normz


        apply!(z, A, q, ONE, ZERO)
        λ = dot(q, z)[1]
        if iter%100 == 0 || iter+1 == niters
            @. resid = z - λ*q
            residual = norm(resid, 2)[1]

            if residual < tolerance
                return (λ, true)
            end
        end
    end
    (λ, false)
end

function log(values...)
    verbose && println(values...)
end

function main(comm::Comm{GID, PID, LID}, numGlobalElements, verbose, Data::Type) where{GID, PID, LID}

    const pid = myPid(comm)
    const nProc = numProc(comm)

    map = BlockMap(numGlobalElements, comm)

    numMyElements = JuliaPetra.numMyElements(map)
    myGlobalElements = JuliaPetra.myGlobalElements(map)
    numNz = Array{GID, 1}(numMyElements)
    for i = 1:numMyElements
        if myGlobalElements[i] == 1 || myGlobalElements[i] == numGlobalElements
            numNz[i] = 2
        else
            numNz[i] = 3
        end
    end

    const A = CSRMatrix{Data}(map, numNz, STATIC_PROFILE)

    const values = Data[-1, -1]
    const two = Data[Data(2)]

    for i = 1:numMyElements
        if myGlobalElements[i] == 1
            indices = LID[2]
        elseif myGlobalElements[i] == numGlobalElements
            indices = LID[numGlobalElements-2]
        else
            indices = LID[myGlobalElements[i]-1, myGlobalElements[i]+1]
        end

        insertGlobalValues(A, myGlobalElements[i], indices, values)
        insertGlobalValues(A, myGlobalElements[i], LID[myGlobalElements[i]], two)
    end

    fillComplete(A, map, map)#use map for domain and range

    const niters = numGlobalElements*10
    const tolerance = 1.0e-2

    log("starting tests")

    #compile all the nessacery methods before timing
    tic()
    λ, success = powerMethod(A, niters, tolerance, verbose)
    elapsedTime = toq()

    log("λ = $λ; is within tolerance? $success")
    log("total time for first solve (plus compile) = $elapsedTime sec\n\n")

    tic()
    #Profile.clear_malloc_data()
    #@time
    λ, success = powerMethod(A, niters, tolerance, verbose)
    #exit(0)
    #TODO look into storing FLOPS
    elapsedTime = toq()


    log("λ = $λ; is within tolerance? $success")
    log("total time for first solve (pre-compiled) = $elapsedTime sec\n\n")
    log("increasing magnitude of first diagonal term, solving again\n")

    if myGID(map, 1)
        rowInds, rowVals = getLocalRowView(A, 1)
        for i = 1:getNumEntriesInGlobalRow(A, 1)
            if rowInds[i] == 1
                #using a view, so values update the original
                rowVals[i] *= 10
            end
        end
    end

    tic()
    λ, success = powerMethod(A, niters, tolerance, verbose)
    elapsedTime = toq()

    log("")

    log("λ = $λ; is within tolerance? $success")
    log("total time for second solve = $elapsedTime sec\n\n")
end



if useMPI
    comm = MPIComm(commGID, commPID, commLID)
else
    comm = SerialComm{commGID, commPID, commLID}()
end

const pid = myPid(comm)
const nProc = numProc(comm)
const verbose = pid == 1

log(comm)

if length(ARGS) != 1
    log("Usage: 1 argument for the number_of_equations")
else
    targetNumGlobalElements = parse(commLID, ARGS[1])
    if targetNumGlobalElements < nProc
        log("numGlobalBlocks = $targetNumGlobalElements cannot be < number of processors = $nProc")
    else
        main(comm, targetNumGlobalElements, verbose, commData)
    end
end
