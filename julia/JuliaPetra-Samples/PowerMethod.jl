
# The power method for finding the largest eigenvalue
# This allows it to be used as a kernal for testing things

using JuliaPetra

function powerMethod(A::RowMatrix{Data, GID, PID, LID}, niters::Integer,
        tolerance) where {Data, GID, PID, LID}
    powerMethod(A, niters, Data(tolerance))
end

"""
Returns a tuple of the computed λ and if the λ is within tolerance
"""
function powerMethod(A::RowMatrix{Data, GID, PID, LID}, niters::Integer,
        tolerance::Data) where {Data, GID, PID, LID}
    q = MultiVector{Data}(getRowMap(A), 1)
    z = MultiVector{Data}(getRowMap(A), 1)
    resid = MultiVector{Data}(getRowMap(A), 1)

    #skipping flop counting

    λ::Data = 0

    randn!(z.data)

    const ONE = one(Data)
    const ZERO = zero(Data)

    iter = 1
    while true
        normz = norm2(z)

        @. q = z/normz

        apply!(z, A, q, ONE, ZERO)
        λ = dot(q, z)[1]
        if true#iter%100 == 0 || iter+1 == niters
            @. resid = z - λ*q
            residual = norm2(resid)[1]

            if residual < tolerance
                return λ
            end
        end
        iter += 1
    end
    @. resid = z - λ*q
    error("could not converge after $niters iteratations.  λ = $λ, residual = $(norm2(resid)[1])")
end
