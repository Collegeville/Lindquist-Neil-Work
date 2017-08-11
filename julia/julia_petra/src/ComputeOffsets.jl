
#used in implementation of CRSGraph, CSRMatrix and, if added, FixedHashTable

function computeOffsets(rowPtrs::Array{<: Integer, 1}, numEnts::Integer)
    numOffsets = length(rowPtrs)
    @simd for i = 1:numOffsets
        @inbounds rowPtrs[i] = numEnts*(i-1)+1
    end
    rowPtrs
end

function computeOffsets(rowPtrs::Array{<: Integer, 1}, numEnts::Array{<: Integer, 1})
    numOffsets = length(rowPtrs)
    numCounts = length(numEnts)
    if numCounts >= numOffsets
        throw(InvalidArgumentError("length(numEnts) = $numCounts "
                * ">= length(rowPtrs) = $numOffsets"))
    end
    sum = 1
    for i = 1:numCounts
        @inbounds rowPtrs[i] = sum
        @inbounds sum += numEnts[i]
    end
    @inbounds rowPtrs[numCounts+1:numOffsets] = sum
    rowPtrs
end