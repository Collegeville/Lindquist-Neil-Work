
export InvalidArgumentError

"""
    InvalidArgumentError(msg)

The values passed as arguments are not valid.  Argument `msg`
is a descriptive error string.
"""
struct InvalidArgumentError <: Exception
    msg::AbstractString
end