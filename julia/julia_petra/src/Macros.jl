#contains random utility macros

globalDebug = false
setGlobalDebug(val::Bool) = (globalDebug = val)

"""
Returns the global debug value

Has an optional ignored argument
"""
macro debug()
    globalDebug
end
