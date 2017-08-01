
#DECISION how to deal with the subtypes having the same implementations for methods, but requiring fields
    #mandate certain fields be present
    #mandate certain getters and setters (or other "low level" operators) be present
    #have each subtype implement (violate Don't Repeat Yourself)

"""
A base type for constructing and using dense multi-vectors, vectors and matrices in parallel.

All subtypes must have the following methods, with Impl standing in for the subtype, in addition to the methods required by SrcDistObject:
doTransfer(src::SrcDistObject{GID, PID, LID}, target::Impl{GID, PID, LID}, cm::CombineMode,
        numSameIDs::LID, permuteToLIDs::Array{LID, 1}, permuteFromLIDs::Array{LID, 1},
        remoteLIDs::Array{LID, 1}, exportLIDs::Array{LID, 1}, distor::Distributor{GID, PID, LID},
        reversed::Bool) where {GID <: Integer, PID <: Integer, LID <: Integer}
    Perform actual redistribution of data across memory images
"""
abstract type DistObject{GID <:Integer, PID <: Integer, LID <: Integer} <: SrcDistObject{GID, PID, LID}
end

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