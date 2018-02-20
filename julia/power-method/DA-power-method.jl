
#Power method implementation using DistributedArrays.jl

@everywhere using DistributedArrays

"""
Returns a tuple of the computed λ and if the λ is within tolerance
"""
function powerMethod(A::DArray{Data, 2}, niters::Integer,
        tolerance::Data)::Tuple{Data, Bool} where Data
    vectDims = (size(A, 1),)
    vectDist = (length(procs(A)),)
    z = drand(Data, vectDims, procs(A), vectDist)::DArray{Data, 1, Array{Data, 1}}
    q = dfill(zero(Data), vectDims, procs(A), vectDist)::DArray{Data, 1, Array{Data, 1}}
    resid = dfill(zero(Data), vectDims, procs(A), vectDist)::DArray{Data, 1, Array{Data, 1}}


    #skipping flop counting

    λ::Data = 0

    for iter = 1:niters
        normz::Data = norm(z, 2)

        #map!(z->z/normz, q, z)
        z, q = q, z
        scale!(q, 1/normz)

        A_mul_B!(Data(1), A, q, Data(0), z)
        λ = dot(q, z)
        if iter%100 == 0 || iter+1 == niters
            #REVIEW mapping and broadcasting 2 dArrays -> 1 dArray not yet supported
            #map!((z,q)->z-λ*q, resid, z, q)
            copy!(resid, z)
            Base.axpy!(-λ, q, resid)

            residual::Data = norm(resid, 2)

            if residual < tolerance
                return (λ, true)
            end
        end
    end
    (λ, false)
end

function main(numGlobalElements)

    Data = Float64

    A = DArray((numGlobalElements, numGlobalElements), workers(), [length(workers()), 1]) do I
        rowIndices = I[1]
        colIndices = I[2]

        if rowIndices.start == 1
            lColIndices = (rowIndices.start:rowIndices.stop-1)
            lRowIndices = 2:length(rowIndices)
        else
            lColIndices = rowIndices-1
            lRowIndices = 1:length(lColIndices)
        end

        if rowIndices.stop == colIndices.stop
            rColIndices = (rowIndices.start+1:rowIndices.stop)
        else
            rColIndices = rowIndices+1
        end
        rRowIndices = 1:length(rColIndices)

        rows = vcat(1:length(rowIndices), lRowIndices, rRowIndices)
        cols = vcat(rowIndices, lColIndices, rColIndices)
        vals = vcat(fill(Data(2), length(rowIndices)), fill(Data(-1), length(lRowIndices)+length(rRowIndices)))
        sparse(rows, cols, vals, length(rowIndices), length(colIndices))
    end

    const niters = numGlobalElements*10
    const tolerance = Data(1.0e-2)

    println("starting tests")

    #compile all the nessacery methods before timing
    tic()
    λ, success = powerMethod(A, niters, tolerance)
    elapsedTime = toq()

    println("λ = $λ; is within tolerance? $success")
    println("total time for first solve (plus compile) = $elapsedTime sec\n\n")

    tic()
    #Profile.clear_malloc_data()
    #@time
    λ, success = powerMethod(A, niters, tolerance)
    #exit(0)
    #TODO look into storing FLOPS
    elapsedTime = toq()


    println("λ = $λ; is within tolerance? $success")
    println("total time for first solve (pre-compiled) = $elapsedTime sec\n\n")
    println("increasing magnitude of first diagonal term, solving again\n")

    #REVIEW this feels really low-level, is there a higher level way that has access to global indices?
    @sync for p in procs(A)
        @async remotecall_fetch(p, A) do A
            if (1 in localindexes(A)[1]) && (1 in localindexes(A)[2])
                DistributedArrays.localpart(A)[1, 1] *= 10
            end
          nothing
      end
    end

    tic()
    λ, success = powerMethod(A, niters, tolerance)
    elapsedTime = toq()

    println("")

    println("λ = $λ; is within tolerance? $success")
    println("total time for second solve = $elapsedTime sec\n\n")
end


if length(ARGS) != 1
    println("Usage: 1 argument for the number_of_equations")
else
    println("Running with $(length(workers())) workers")
    targetNumGlobalElements = parse(Int64, ARGS[1])
    main(targetNumGlobalElements)
end
