export CombineMode, ADD, INSERT, REPLACE, ABSMAX, ZERO
export TransposeMode, NO_TRANS, TRANS, CONJ_TRANS

"""
Tells petra how to combine data received from other processes with existing data on the calling process for specific import or export options.

Here is the list of combine modes:
  * ADD: Sum new values into existing values
  * INSERT: Insert new values that don't currently exist
  * REPLACE: REplace existing values with new values
  * ABSMAX: If ``x_{old}`` is the old value and ``x_{new}`` the incoming new value, replace ``x_{old}`` with ``\\max(x_{old}, x_{new})``
  * ZERO: Replace old values with zero
"""
@enum CombineMode ADD=1 INSERT=2 REPLACE=3 ABSMAX=4 ZERO=5



"""
Tells petra whether to use the transpose or conjugate transpose of the matrix
"""
@enum TransposeMode NO_TRANS=1 TRANS=2 CONJ_TRANS=3

"""
    isTransposed(mode::TransposeMode)::Bool

Checks whether the given TransposeMode is transposed
"""
@inline function isTransposed(mode::TransposeMode)::Bool
    mode != TransposeMode.NO_TRANS
end