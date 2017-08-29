export MultiVector
export localLength, globalLength, numVectors, map
export scale!, scale
export getVectorView, getVectorCopy
export commReduce


#TODO implement dot(::MultiVector{...}, ::MultiVector{...})::Array{Data, 1}
#TODO implement norm2(::MultiVector{...})::Array{Data, 1}

"""
MultiVector represents a dense multi-vector.  Note that all the vectors in a single MultiVector are the same size
"""
type MultiVector{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistObject{GID, PID, LID}
    data::Array{Data, 2} # data[1, 2] is the first element of the second matrix
    localLength::LID
    globalLength::GID
    numVectors::LID
    
    map::BlockMap{GID, PID, LID}
end

## Constructors ##

"""
    MultiVector{Data, GID, PID, LID}(::BlockMap{GID, PID, LID}, numVecs::Integer, zeroOut=True)

Creates a new MultiVector based on the given map
"""
function MultiVector{Data, GID, PID, LID}(map::BlockMap{GID, PID, LID}, numVecs::Integer, zeroOut=True) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    localLength = numMyElements(map)
    if zeroOut
        data = zeros(Data, (localLength, numVecs))
    else
        data = Array{Data, 2}(localLength, numVecs)
    end
    MultiVector{Data, GID, PID, LID}(data, localLength, numGlobalElements(map), numVecs, map)
end

"""
    MultiVector{Data, GID, PID, LID}(map::BlockMap{GID, PID, LID}, data::Array{Data, 2})

Creates a new MultiVector wrapping the given array.  Changes to the MultiVector or Array will affect the other
"""
function MultiVector{Data, GID, PID, LID}(map::BlockMap{GID, PID, LID}, data::Array{Data, 2}) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    localLength = numMyElements(map)
    if size(data, 1) != localLength
        throw(InvalidArgumentError("Length of vectors does not match local length indicated by map"))
    end
    MultiVector{Data, GID, PID, LID}(data, localLength, numGlobalElements(map), size(data, 2), map)
end

## External methods ##

"""
    copy(::MutliVector{Data, GID, PID, LID})::MultiVector{Data, GID, PID, LID}
Returns a copy of the multivector
"""
function Base.copy(vect::MultiVector{Data, GID, PID, LID})::MultiVector{Data, GID, PID, LID} where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    MultiVector{Data, GID, PID, LID}(copy(vect.data), vect.localLength, vect.globalLength, vect.numVectors, vect.map)
end

function Base.copy!(dest::MultiVector{Data, GID, PID, LID}, src::MultiVector{Data, GID, PID, LID})::MultiVector{Data, GID, PID, LID} where {Data, GID, PID, LID}
    copy!(dest.data, src.data)
    dest.localLength = src.localLength
    dest.globalLength = src.globalLength
    dest.numVectors = src.numVectors
    dest.map = src.map
    
    dest
end

"""
    localLength(::MutliVector{Data, GID, PID, LID})::LID

Returns the local length of the vectors in the multivector
"""
function localLength(vect::MultiVector{Data, GID, PID, LID})::LID where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    vect.localLength
end

"""
    globalLength(::MultiVector{Data, GID, PID, LID})::GID

Returns the global length of the vectors in the mutlivector
"""
function globalLength(vect::MultiVector{Data, GID})::GID where {Data <: Number, GID <: Integer}
    vect.globalLength
end

"""
    numVectors(::MultiVector{Data, GID, PID, LID})::LID

Returns the number of vectors in this multivector
"""
function numVectors(vect::MultiVector{Data, GID, PID, LID})::LID where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    vect.numVectors
end

"""
    numVectors(::MultiVector{Data, GID, PID, LID})::BlockMap{GID, PID, LID}

Returns the BlockMap used by this multivector
"""
function map(vect::MultiVector{Data, GID, PID, LID})::BlockMap{GID, PID, LID} where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    vect.map
end


# have to use Base.scale! to avoid requiring module qualification everywhere
"""
    scale!(::MultiVector{Data, GID, PID, LID}, ::Data})::MultiVector{Data, GID, PID, LID}

Scales the mulitvector in place and returns it
"""
function Base.scale!(vect::MultiVector{Data, GID, PID, LID}, alpha::Data)::MultiVector{Data, GID, PID, LID} where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    @. vect.data = vect.data*alpha
    vect
end

"""
    scale!(::MultiVector{Data, GID, PID, LID}, ::Data)::MultiVector{Data, GID, PID, LID}

Scales a copy of the mulitvector and returns the copy
"""
function scale(vect::MultiVector{Data, GID, PID, LID}, alpha::Data)::MultiVector{Data, GID, PID, LID} where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    scale!(copy(vect), alpha)
end

"""
    scale!(::MultiVector{Data, GID, PID, LID}, ::Array{Data, 1})::MultiVector{Data, GID, PID, LID}

Scales each column of the mulitvector in place and returns it
"""
function Base.scale!(vect::MultiVector{Data, GID, PID, LID}, alpha::Array{Data, 1})::MultiVector{Data, GID, PID, LID} where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    for v = 1:vect.numVectors
        vect.data[:, v] *= alpha[v]
    end
    vect
end

"""
    scale(::MultiVector{Data, GID, PID, LID}, ::Array{Data, 1})::MultiVector{Data, GID, PID, LID}

Scales each column of a copy of the mulitvector and returns the copy
"""
function scale(vect::MultiVector{Data, GID, PID, LID}, alpha::Array{Data, 1})::MultiVector{Data, GID, PID, LID} where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    scale!(copy(vect), alpha)
end


"""
    getVectorView(::MultiVector{Data}, columns)::AbstractArray{Data}

Gets a view of the requested column vector(s) in this multivector
"""
function getVectorView(mVect::MultiVector{Data}, column)::AbstractArray{Data}  where {Data}
    view(mVect.data, :, column)
end

"""
    getVectorCopy(::MultiVector{Data}, columns)::Array{Data}

Gets a copy of the requested column vector(s) in this multivector
"""
function getVectorCopy(mVect::MultiVector{Data}, column)::Array{Data} where {Data}
    mVect.data[:, column]
end

function Base.fill!(mVect::MultiVector, values)
    fill!(mVect.data, values)
    mVect
end

"""
    commReduce(::MultiVector)

Reduces the content of the MultiVector across all processes.  Note that the MultiVector cannot be distributed globally.
"""
function commReduce(mVect::MultiVector)
    if distributedGlobal(mVect)
        throw(InvalidArgumentError("Cannot reduce distributed MultiVector"))
    end
    
    mVect.data = sumAll(comm(mVect), mVect.data)
end


## DistObject interface ##

function checkSizes(source::MultiVector{Data, GID, PID, LID},
        target::MultiVector{Data, GID, PID, LID})::Bool where {
            Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    (source.numVectors == target.numVectors 
        && source.globalLength == target.globalLength 
        && source.localLength == target.localLength)
end


function copyAndPermute(source::MultiVector{Data, GID, PID, LID},
        target::MultiVector{Data, GID, PID, LID}, numSameIDs::LID,
        permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1}
        ) where {Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    target.data[1:numSameIDs, :] = source.data[1:numSameIDs, :]
    #don't need to sort permute[To/From]LIDs, since the orders match
    target.data[permuteToLIDs, :] = source.data[permuteFromLIDs, :]
end

function packAndPrepare(source::MultiVector{Data, GID, PID, LID},
        target::MultiVector{Data, GID, PID, LID}, exportLIDs::Array{LID, 1},
        distor::Distributor{GID, PID, LID})::Array where {
            Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    exports = Array{Array{Data, 1}}(length(exportLIDs))
    for i = 1:length(exports)
        exports[i] = source.data[exportLIDs[i], :]
    end
    exports
end

function unpackAndCombine(target::MultiVector{Data, GID, PID, LID},
        importLIDs::Array{LID, 1}, imports::Array,
        distor::Distributor{GID, PID, LID},cm::CombineMode) where {
            Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer}
    for i = 1:length(importLIDs)
        target.data[importLIDs[i], :] = imports[i]
    end
end
