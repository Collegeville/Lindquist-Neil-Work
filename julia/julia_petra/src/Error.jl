
export InvalidArgumentError

"""
    ArgumentValueError(msg)

The values passed as arguments are not valid.  Argument `msg`
is a descriptive error string.
"""
immutable InvalidArgumentError <: Exception
    msg::AbstractString
end