
struct DiagonalInversePreconditioner{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    mat::RowMatrix{Data, GID, PID, LID}
    importer::Nullable{Import{GID, PID, LID}}
end

function DiagonalInversePreconditioner(mat::RowMatrix{Data, GID, PID, LID}) where {Data, GID, PID, LID}
    #currently assuming range and domain are the same
    @assert sameAs(getRangeMap(mat), getDomainMap(mat))

    DiagonalInversePreconditioner(mat, Nullable{Import{GID, PID, LID}}())
end

getRangeMap(prec::DiagonalInversePreconditioner) = getDomainMap(prec.mat)
getDomainMap(prec::DiagonalInversePreconditioner) = getRangeMap(prec.mat)

function apply!(Y::MultiVector{Data, GID, PID, LID},
        operator::DiagonalInversePreconditioner{Data, GID, PID, LID},
        X::MultiVector{Data, GID, PID, LID},
        mode::TransposeMode, alpha::Data, beta::Data) where {Data, GID, PID, LID}

    if alpha == Data(0)
        if beta == Data(0)
            fill!(Y, Data(0))
        elseif beta != Data(1)
            scale!(Y, beta)
        end
        return Y
    end

    YIsOverwritten = (beta == ZERO)
    YIsReplicated = !distributedGlobal(Y) && numProc(getComm(Y)) != 0

    #part of special case for replicated MV output
    if YIsReplicated && myPid(getComm(Y)) != 1
        beta = Data(0)
    end

    diag = getLocalDiagCopy(operator.mat)

    rawX = X.data
    rawY = Y.data

    numRows = getLocalNumRows(operator.mat)
    for vect = LID(1):numVectors(Y)
        for row = LID(1):numRows
            @inbounds rawY[row, vect] = rawY[row, vect]*beta + applyConjugation(mode, rawX[row, vect]/diag[row]*alpha)
        end
    end

    if YIsReplicated
        commReduce(Y)
    end
    Y
end
