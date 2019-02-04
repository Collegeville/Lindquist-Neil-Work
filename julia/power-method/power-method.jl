using JuliaPetra
using Random
import Profile

include("MMToJP.jl")

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
        tolerance::Data, verbose::Bool) where {Data, GID, PID, LID}
    q = DenseMultiVector{Data}(getRowMap(A), 1)
    z = DenseMultiVector{Data}(getRowMap(A), 1)
    resid = DenseMultiVector{Data}(getRowMap(A), 1)

    #skipping flop counting

    λ::Data = 0

    randn!(z.data)

    ONE = one(Data)
    ZERO = zero(Data)

    iter = 1
    while iter <= niters
        normz = norm(z, 2)[1]

        @. q = z/normz


        apply!(z, A, q, ONE, ZERO)
        λ = dot(q, z)[1]
        if iter%100 == 0 || iter+1 == niters
            @. resid = z - λ*q
            residual = norm(resid, 2)[1]

            if residual < tolerance
                return (λ, true, iter)
            end
        end
        iter += 1
    end
    (λ, false, iter)
end

function log(values...)
    verbose && println(values...)
end

function main(comm::Comm{GID, PID, LID}, verbose, Data::Type) where{GID, PID, LID}

    #try
        A = readMM("power-method/matrix.mm", comm)
    #catch ex
    #    println("caught error: $ex")

    #    return
    #end

    niters = getGlobalNumRows(A)*10
    tolerance = 1.0e-2

    log("starting tests")

    #compile all the nessacery methods before timing

    elapsedTime = @elapsed λ, success, iter = powerMethod(A, niters, tolerance, verbose)

    log("λ = $λ; is within tolerance? $success")
    log("Took $iter iterations")
    log("total time for first solve (plus compile) = $elapsedTime sec\n\n")

    elapsedTime = @elapsed λ, success, iter = powerMethod(A, niters, tolerance, verbose)
    #elapsedTime = 0
    #Profile.init(;n=10000000)
    #Profile.clear()
    #Profile.@profile begin
    #    elapsedTime = @elapsed λ, success, iter = powerMethod(A, niters, tolerance, verbose)
    #end
    #verbose && Profile.print(;mincount=5)
    #TODO look into storing FLOPS


    log("λ = $λ; is within tolerance? $success")
    log("Took $iter iterations")
    log("total time for first solve (pre-compiled) = $elapsedTime sec\n\n")
    log("increasing magnitude of first diagonal term, solving again\n")

    if myGID(getRowMap(A), 1)
        rowInds, rowVals = getLocalRowView(A, 1)
        for i = 1:getNumEntriesInGlobalRow(A, 1)
            if rowInds[i] == 1
                #using a view, so values update the original
                rowVals[i] *= 10
            end
        end
    end

    elapsedTime = @elapsed λ, success, iter = powerMethod(A, niters, tolerance, verbose)

    log("")

    log("λ = $λ; is within tolerance? $success")
    log("Took $iter iterations")
    log("total time for second solve = $elapsedTime sec\n\n")
end



if useMPI
    comm = MPIComm(commGID, commPID, commLID)
else
    comm = SerialComm{commGID, commPID, commLID}()
end

const pid = myPid(comm)
const verbose = pid == 1
log(comm)
main(comm, verbose, commData)
