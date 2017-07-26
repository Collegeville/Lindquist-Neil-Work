
export getDirectoryEntries, gidsAllUniquelyOwned
export createDirectory

function getDirectoryEntries(directory::Directory, map::BlockMap,
        globalEntries::Array{GID})::Tuple{Array, Array} where GID <: Integer
        # where PID <: Integer where LID <: Integer
    getDirectoryEntries(directory, map, globalEntries, false)
end


"""
    createDirectory(comm::Comm, map::BlockMap)
Create a directory object for the given Map
"""
function createDirectory(comm::Comm, map::BlockMap)::BasicDirectory
    BasicDirectory(map)
end