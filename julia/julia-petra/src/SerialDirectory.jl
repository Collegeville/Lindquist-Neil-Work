

type SerialDirectory <: Directory
    map::Map
end

function getDirectoryEntries(directory::SerialDirectory, map::Map,
        globalEntries::Array{GID}, high_rank_sharing_procs::Bool
        )::Tuple{Array{Integer}, Array{Integer}} where GID <: Integer
        #where PID <: Integer where LID <:Integer
    numEntries = length(globalEntries)
    procs = Array{Integer}(numEntries)
    localEntries = Array{Integer}(numEntries)
    
    for i = 1:numEntries
        lid = map.lid(globalEntries[i])
        
        if lid == 0
            procs[i] = 0
            warn("GID $globalEntries[i] is not part of this map")
        else
            procs[i] = 1
        end
        localEntries[i] = lid
    end
    
    procs, localEntries
end
    

function gidsAllUniquelyOwned(directory::SerialDirectory)
    true
end