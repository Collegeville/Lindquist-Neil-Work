
# methods and docs based straight off Epetra_Directory to match Comm


# DECISION: figure out types for global and local IDs
# DECISION: figure out types for entry sizes
"""
A base type as an interface to allow Map and BlockMap objects to reerence non-local
elements.

All subtypes must have the following methods, with DirectoryImpl standing in for
the subtype:

getDirectoryEntries(directory::DirectoryImpl, Map::BlockMap, globalEntries::Array{GID},
        high_rank_sharing_procs::Bool) where GID <: Integer
    - Returns processor and local id infor for non-local map entries.  Returns a tuple
        containing
            1 - an Array of processors owning the global ID's in question
            2 - an Array of local IDs of the global on the owning processor
            3 - an Array containing the size of the objects associated with each global id

gidsAllUniquelyOwned()
    - Returns true if all GIDs appear on just one processor
"""
abstract type Directory
end


getDirectoryEntries(directory::Directory, map::BlockMap, globalEntries::Array{GID})
        where GID <: Integer
    getDirectoryEntries(directory, map, globalEntries, false)
end