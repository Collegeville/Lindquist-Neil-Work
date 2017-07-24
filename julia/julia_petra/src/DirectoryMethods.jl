
export getDirectoryEntries, gidsAllUniquelyOwned

function getDirectoryEntries(directory::Directory, map::BlockMap,
        globalEntries::Array{GID})::Tuple{Array, Array} where GID <: Integer
        # where PID <: Integer where LID <: Integer
    getDirectoryEntries(directory, map, globalEntries, false)
end
