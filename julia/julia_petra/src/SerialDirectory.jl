
export SerialDirectory

type SerialDirectory <: Directory
    map::BlockMap
end

function getDirectoryEntries(directory::SerialDirectory, map::BlockMap,
        globalEntries::Array{GID}, high_rank_sharing_procs::Bool
        )::Tuple{Array{Integer}, Array{Integer}} where GID <: Integer
        #where PID <: Integer where LID <:Integer
    numEntries = length(globalEntries)
    procs = Array{Integer}(numEntries)
    localEntries = Array{Integer}(numEntries)
    
    for i = 1:numEntries
        lidVal = lid(map, globalEntries[i])
        
        if lidVal == 0
            procs[i] = 0
            warn("GID $(globalEntries[i]) is not part of this map")
        else
            procs[i] = 1
        end
        localEntries[i] = lidVal
    end
    
    (procs, localEntries)
end
    

function gidsAllUniquelyOwned(directory::SerialDirectory)
    true
end