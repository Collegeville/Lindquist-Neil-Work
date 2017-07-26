
export getDirectoryEntries, gidsAllUniquelyOwned
export createDirectory

function getDirectoryEntries(directory::Directory{GID, PID, LID}, map::BlockMap{GID, PID, LID},
        globalEntries::Array{GID})::Tuple{Array{PID}, Array{LID}} where GID <: Integer where PID <: Integer where LID <: Integer
    getDirectoryEntries(directory, map, globalEntries, false)
end


"""
    createDirectory(comm::Comm, map::BlockMap)
Create a directory object for the given Map
"""
function createDirectory(comm::Comm{GID, PID, LID}, map::BlockMap{GID, PID, LID})::BasicDirectory{GID, PID, LID} where GID <: Integer where PID <: Integer where LID <: Integer
    BasicDirectory{GID, PID, LID}(map)
end
