
#Power method implementation using Shared Arrays

@everywhere using ParallelSparseMatMul

"""
Computes the power method on the transpose of A
Returns a tuple of the computed λ and if the λ is within tolerance
"""
function powerMethod(A::AbstractSparseMatrix{Data}, niters::Integer,
        tolerance::Data)::Tuple{Data, Bool} where Data
    const numElements = size(A, 1)

    q = SharedArray{Data}(numElements; init = S-> S[Base.localindexes(S)]=0, pids=workers())
    z = SharedArray{Data}(numElements; init = S-> S[Base.localindexes(S)]=randn(length(Base.localindexes(S))), pids=workers())
    resid = SharedArray{Data}(numElements; init = S-> S[Base.localindexes(S)]=0, pids=workers())

    #skipping flop counting

    λ::Data = 0

    for iter = 1:niters
        normz::Data = norm(z, 2)

        #map!(z->z/normz, q, z)
        @. q = z/normz
        #@sync @parallel for i in 1:numElements
        #    q[i] = z[i]/normz
        #end

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

    rows = vcat(1:numGlobalElements, 2:numGlobalElements, 1:numGlobalElements-1)
    cols = vcat(1:numGlobalElements, 1:numGlobalElements-1, 2:numGlobalElements)
    vals = vcat(fill(Data(2), numGlobalElements), fill(Data(-1), numGlobalElements*2-2))
    A = share(sparse(rows, cols, vals))

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

    temp = sdata(A)
    temp[1, 1] *= 10
    A = share(temp)

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
