
#TODO figure out Packet
"""
A base type for constructing and using dense multi-vectors, vectors and matrices in parallel.

All subtypes must have the following methods, with Impl standing in for the subtype, in addition to the methods required by SrcDistObject:

checkSizes(source::SrcDistObject{GID, PID, LID}, target::Impl{GID, PID, LID})::Bool
    Compare the source and target objects for compatiblity

copyAndPermute(source::SrcDistObject{GID, PID, LID}, target::Impl{GID, PID, LID},
        numSameIDs::LID, permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1})
    Perform copies and permutations that are local to this process.

packAndPrepare(source::SrcDistObject{GID, PID, LID}, target::Impl{GID, PID, LID},
        exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID}
        )::Tuple{Array{Packet}, Union{Array{Integer}, Integer}}
    Perform any packing or preparation required for communications

unpackAndCombine(target::Impl{GID, PID, LID}, importLIDs::Array{LID, 1}, imports::Array{Packet},
        numPacketsPerLID::Union{Array{Integer}, Integer}}, distor::Distributor{GID, PID, LID},
        cm::CombineMode)
    Perform any unpacking and combining after communication
"""
abstract type DistObject{GID <:Integer, PID <: Integer, LID <: Integer} <: SrcDistObject{GID, PID, LID}
end


## import/export interface ##

"""
    doImport(target::Impl{GID, PID, LID}, source::SrcDistObject{GID, PID, LID}, importer::Import{GID, PID, LID}, cm::CombineMode)
Import data into this object using an Import object ("forward mode")
"""
function doImport(source::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID},
        importer::Import{GID, PID, LID}, cm::CombineMode) where {
        GID <:Integer, PID <: Integer, LID <: Integer}
    #TODO add checks for map equality when debuging
    
    doTransfer(source, target, cm, numSameIDs(importer), permuteToLIDs(importer),
        permuteFromLIDs(importer), remoteLIDs(importer), exportLIDs(importer),
        distributor(importer), false)
end

"""
    doExport(target::Impl{GID, PID, LID}, source::SrcDistObject{GID, PID, LID}, exporter::Export{GID, PID, LID}, cm::CombineMode)
Export data into this object using an Export object ("forward mode")
"""
function doExport(source::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID},
        exporter::Export{GID, PID, LID}, cm::CombineMode) where {
        GID <:Integer, PID <: Integer, LID <: Integer}
    #TODO add checks for map equality when debuging
    
    doTransfer(source, target, cm, numSameIDs(exporter), permuteToLIDs(exporter),
        permuteFromLIDs(exporter), remoteLIDs(exporter), exportLIDs(exporter),
        distributor(exorter), false)
end

"""
    doImport(source::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID}, exporter::Export{GID, PID, LID}, cm::CombineMode)
Import data into this object using an Export object ("reverse mode")
"""
function doImport(source::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID},
        exporter::Export{GID, PID, LID}, cm::CombineMode) where {
            GID <:Integer, PID <: Integer, LID <: Integer}
    #TODO add checks for map equality when debuging
    
    doTransfer(source, target, cm, numSameIDs(exporter), permuteToLIDs(exporter),
        permuteFromLIDs(exporter), remoteLIDs(exporter), exportLIDs(exporter),
        distributor(exporter), true)
end

"""
    doExport(source::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID}, importer::Import{GID, PID, LID}, cm::CombineMode)
Export data into this object using an Import object ("reverse mode")
"""
function doExport(source::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID},
        importer::Import{GID, PID, LID}, cm::CombineMode) where {
            GID <:Integer, PID <: Integer, LID <: Integer}
    #TODO add checks for map equality when debugin
    
    doTransfer(source, target, cm, numSameIDs(exporter), permuteToLIDs(exporter),
        permuteFromLIDs(exporter), remoteLIDs(exporter), exportLIDs(exporter),
        distributor(exporter), true)
end


## import/export functionality ##

"""
    constantNumberOfPackets(::Impl)::Integer
Whether the implementation's instance promises always to have a constant
number of packets per LID (local index), and if so, how many packets
per LID there are.  ``0`` indicates the number of packets may not be consistant.

The default implementation returns 0, but may be overrided by subtypes to reduce
allocation
"""
function constantNumberOfPackets(obj::DistObject)::Integer
    0
end

"""
    doTransfer(src::SrcDistObject{GID, PID, LID}, target::Impl{GID, PID, LID}, cm::CombineMode, numSameIDs::LID, permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1}, remoteLIDs::Array{LID, 1}, exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID}, reversed::Bool)
Perform actual redistribution of data across memory images
"""
function doTransfer(src::SrcDistObject{GID, PID, LID}, target::DistObject{GID, PID, LID},
        cm::CombineMode, numSameIDs::LID, permuteToLIDs::Array{LID, 1},
        permuteFromLIDs::Array{LID, 1}, remoteLIDs::Array{LID, 1}, exportLIDs::Array{LID, 1},
        distor::Distributor{GID, PID, LID}, reversed::Bool) where {
            GID <: Integer, PID <: Integer, LID <: Integer}
    
    #TODO implement 
end