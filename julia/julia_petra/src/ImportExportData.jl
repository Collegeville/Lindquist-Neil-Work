

type ImportExportData{GID <: Integer, PID <: Integer, LID <: Integer}
    source::BlockMap{GID, PID, LID}
    target::BlockMap{GID, PID, LID}
    
    permuteToLIDs::Array{LID, 1}
    permuteFromLIDs::Array{LID, 1}
    remoteLIDs::Array{LID, 1}
    
    exportLIDs::Array{LID, 1}
    exportPIDs::Array{PID, 1}
    
    numSameIDs::GID
    distributor::Distributor
    
    isLocallyComplete::Bool
end

## Constructors ##
function ImportExportData(source::BlockMap{GID, PID, LID}, target::BlockMap{GID, PID, LID})::ImportExportData{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    ImportExportData{GID, PID, LID}(source, target, [], [], [], [], [], 0, getDistributor(comm(source)), true)
end
    
## Getters ##
"""
Get the source map for the given ImportExportData
"""
function sourceMap(data::ImportExportData{GID, PID, LID})::BlockMap{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    data.source
end

"""
Get the target map for the given ImportExportData
"""
function targetMap(data::ImportExportData{GID, PID, LID})::BlockMap{GID, PID, LID} where GID <: Integer where PID <:Integer where LID <: Integer
    data.target
end

"""
List of elements in the target map that are permuted.
"""
function permuteToLIDs(data::ImportExportData{GID, PID, LID})::Array{LID} where GID <: Integer where PID <:Integer where LID <: Integer
    data.permuteToLIDs
end

"""
List of elements in the source map that are permuted.
"""
function permuteFromLIDs(data::ImportExportData{GID, PID, LID})::Array{LID} where GID <: Integer where PID <:Integer where LID <: Integer
    data.permuteFromLIDs
end

"""
List of elements in the target map that are coming from other processors
"""
function remoteLIDs(data::ImportExportData{GID, PID, LID})::Array{LID} where GID <: Integer where PID <: Integer where LID <: Integer
    data.remoteLIDs
end

"""
Sets the list of elements in the target map that are coming from other processors
"""
function remoteLIDs(data::ImportExportData{GID, PID, LID}, remoteLIDs::Array{LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    data.remoteLIDs = remoteLIDs
end

"""
List of elements that will be sent to other processors
"""
function exportLIDs(data::ImportExportData{GID, PID, LID})::Array{LID} where GID <: Integer where PID <: Integer where LID <: Integer
    data.exportLIDs
end

"""
Sets the list of elements that will be sent to other processors
"""
function exportLIDs(data::ImportExportData{GID, PID, LID}, exportLIDs::Array{LID}) where GID <: Integer where PID <: Integer where LID <: Integer
    data.exportLIDs = exportLIDs
end

"""
List of processors to which elements will be sent `exportLID[i]` will be sent to processor `exportPIDs[i]`
"""
function exportPIDs(data::ImportExportData{GID, PID, LID})::Array{PID} where GID <: Integer where PID <: Integer where LID <: Integer
    data.exportPIDs
end

"""
Sets the list of processors to which elements will be sent `exportLID[i]` will be sent to processor `exportPIDs[i]`
"""
function exportPIDs(data::ImportExportData{GID, PID, LID}, exportPIDs::Array{PID}) where GID <: Integer where PID <: Integer where LID <: Integer
    data.exportPIDs = exportPIDs
end

"""
Returns the number of elements that are identical between the source and target maps, up to the first different ID
"""
function numSameIDs(data::ImportExportData{GID, PID, LID})::LID where GID <: Integer where PID <: Integer where LID <: Integer
    data.numSameIDs
end

"""
Sets the number of elements that are identical between the source and target maps, up to the first different ID
"""
function numSameIDs(data::ImportExportData{GID, PID, LID}, numSame::GID)::LID where GID <: Integer where PID <: Integer where LID <: Integer
    data.numSameIDs = numSame
end


"""
Returns the distributor being used
"""
function distributor(data::ImportExportData{GID, PID, LID})::Distributor{GID, PID, LID} where GID <: Integer where PID <: Integer where LID <: Integer
    data.distributor
end

"""
Returns whether the import or export is locally complete
"""
function isLocallyComplete(data::ImportExportData{GID, PID, LID})::Bool where GID <: Integer where PID <: Integer where LID <: Integer
    data.isLocallyComplete
end

"""
Sets whether the import or export is locally complete
"""
function isLocallyComplete(data::ImportExportData{GID, PID, LID}, isLocallyComplete::Bool) where GID <: Integer where PID <: Integer where LID <: Integer
    data.isLocallyComplete = isLocallyComplete
end