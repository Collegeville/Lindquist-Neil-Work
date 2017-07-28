
#DECISION how to deal with the subtypes having the same implementations for methods, but requiring fields
    #mandate certain fields be present
    #mandate certain getters and setters (or other "low level" operators) be present
    #have each subtype implement (violate Don't Repeat Yourself)

#TODO figure out args for import and export
"""
A base type for constructing and using dense multi-vectors, vectors and matrices in parallel.

All subtypes must have the following methods, with Impl standing in for the subtype:

"""
abstract type DistObject{GID <:Integer, PID <: Integer, LID <: Integer} <: SrcDistObject{GID, PID, LID}
end
