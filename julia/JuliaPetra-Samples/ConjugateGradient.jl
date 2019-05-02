
using JuliaPetra

function conjugate_gradient(A, b::MultiVector{Data, GID, PID, LID}
        ; niters = 100, tol = 1e-6) where {Data, GID, PID, LID}
    conjugate_gradient(b, A, b; niters=niters, tol=tol)
end

"""
Unpreconditioned conjugate gradient
"""
function conjugate_gradient(x::MultiVector{Data, GID, PID, LID}, A, b::MultiVector{Data, GID, PID, LID},
        ; niters = 100, tol = 1e-6) where {Data, GID, PID, LID}
    x = copy(x)
    conjugate_gradient!(x, A, b; niters=niters, tol=tol)
    x
end

"""
Unpreconditioned conjugate gradient
x is modified in place
"""
function conjugate_gradient!(x::MultiVector{Data, GID, PID, LID}, A, b::MultiVector{Data, GID, PID, LID},
        ; niters = 100, tol = 1e-6) where {Data, GID, PID, LID}
    Ax = apply(b, A, x)
    r = apply(b, A, x, Data(-1))
    p = copy(r)
    rtr = dot(r, r)

    Ap = MultiVector{Data}(getMap(b), numVectors(b), false)

    for k in 1:niters
        apply!(Ap, A, p)
        pAp = dot(p, Ap)
        α = rtr ./ pAp

        @. x = x + α * p
        @. r = r - α * Ap

        rtr_old = rtr
        rtr = dot(r, r)
        if !any(sqrt.(rtr) .> tol)
            break
        end
        β = rtr ./ rtr_old
        @. p = r + β * p
    end
    x
end

function conjugate_gradient(A, b::MultiVector{Data, GID, PID, LID}, Minv
        ; niters = 100, tol = 1e-6) where {Data, GID, PID, LID}
    conjugate_gradient(b, A, b, Minv; niters=niters, tol=tol)
end

"""
Preconditioned conjugate gradient
"""
function conjugate_gradient(x::MultiVector{Data, GID, PID, LID}, A, b::MultiVector{Data, GID, PID, LID},
        Minv; niters = 100, tol = 1e-6) where {Data, GID, PID, LID}
    x = copy(x)
    conjugate_gradient!(x, A, b, Minv; niters=niters, tol=tol)
    x
end

"""
Preconditioned conjugate gradient
x is modified in place
"""
function conjugate_gradient!(x::MultiVector{Data, GID, PID, LID}, A, b::MultiVector{Data, GID, PID, LID},
        Minv; niters = 100, tol = 1e-6) where {Data, GID, PID, LID}
    r = apply(b, A, x, Data(-1), Data(1))
    z = MultiVector{Data}(getMap(x), numVectors(x), false)
    apply!(z, Minv, r)
    p = copy(z)

    Ap = MultiVector{Data}(getMap(b), numVectors(b), false)

    rtz = dot(r, z)

    for k in 1:niters
        apply!(Ap, A, p)
        pAp = dot(p, Ap)
        α = rtz ./ pAp

        @. x = x + α * p
        @. r = r - α * Ap

        apply!(z, Minv, r)

        rtz_old = rtz
        rtz = dot(r, z)
        if !any(sqrt.(rtz) .> tol)
            break
        end
        β = rtz/rtz_old
        @. p = z + β * p
    end
    x
end
