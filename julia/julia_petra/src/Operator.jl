
export apply!, apply, @operatorFunctions

"""
Operator is a description of all types that have a specific set of methods.  ``@operatorFunctions typ`` must be called
for Operator type ``typ``.  Operators must have 4 parametric types:
    Data - the type of the data
    GID  - the type of the global indexes
    PID  - the type of the processor ranks
    LID  - the type of the local indexes

All Operator types must implement the following methods (with Op standing in for the Operator):

apply!(operator::Op{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, mode::TransposeMode, alpha::Data, beta::Data)
    Computes ``Y = α\cdot A^{mode}\cdot X + β\cdot Y``, with the following exceptions
        If beta == 0, apply MUST overwrite Y, so that any values in Y (including NaNs) are ignored.
        If alpha == 0, apply MAY short-circuit the operator, so that any values in X (including NaNs) are ignored


domainMap(operator::Op{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    Returns the BlockMap associated with the domain of this operation

rangeMap(operator::Op{Data, GID, PID, LID})::BlockMap{GID, PID, LID}
    Returns the BlockMap associated with the range of this operation
The field operators contains all types that have had ``@operatorFunctions`` called on them
"""
const operators = []


"""
    apply!(operator, X::MultiVector, Y::MultiVector, mode::TransposeMode=NO_TRANS, alpha=1, beta=0)
    apply!(operator, X::MultiVector, Y::MultiVector, alpha=1, beta=0)

Computes ``Y = α\cdot A^{mode}\cdot X + β\cdot Y``, with the following exceptions:
* If beta == 0, apply MUST overwrite Y, so that any values in Y (including NaNs) are ignored.
* If alpha == 0, apply MAY short-circuit the operator, so that any values in X (including NaNs) are ignored
"""
function apply! end



@createFunctionGenerator operatorFunctions begin
    function apply!(matrix::$typ{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, mode::TransposeMode=TransposeMode.NO_TRANS, alpha::Data=1) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
        apply!(matrix, X, Y, mode, alpha, 0)
    end

    function apply!(matrix::$typ{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, alpha::Data, beta::Data=0) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
        apply!(matrix, X, Y, Transpose_Mode.NO_TRANS, alpha, beta)
    end

    """
    Computes ``Y' = α\cdot A^{mode}\cdot X + β\cdot Y`` and returns ``Y'``, with the following exceptions
        If beta == 0, apply MUST overwrite Y, so that any values in Y (including NaNs) are ignored.
        If alpha == 0, apply MAY short-circuit the operator, so that any values in X (including NaNs) are ignored
    """
    function apply(matrix::$typ{Data, GID, PID, LID}, X::MultiVector{Data, GID, PID, LID}, Y::MultiVector{Data, GID, PID, LID}, mode::TransposeMode=NO_TRANS, alpha::Data=1, beta=0) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
        Y = copy(Y)
        apply!(matrix, X, Y, mode, alpha, beta)
    end
    
    push!(operators, $typ)
end
