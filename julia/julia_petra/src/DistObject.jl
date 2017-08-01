
#TODO add getter and setter methods for fields required by transfer in requirements

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
        )::Tuple{Array{Packet}, Union{Array{<:Integer}, Integer}}
    Perform any packing or preparation required for communications.  The returned tuple contains
        1 the allocated export array
        2 if constantPacketsPerLID != 0: constantPacketsPerLID else: array of number of Packets per LID

unpackAndCombine(target::Impl{GID, PID, LID}, importLIDs::Array{LID, 1}, imports::Array{Packet},
        numPacketsPerLID::Union{Array{<:Integer}, Integer}}, distor::Distributor{GID, PID, LID},
        cm::CombineMode)
    Perform any unpacking and combining after communication

numExportPacketsPerLID(target::Impl, numExportPackets::Array{<:Integer})
numExportPacketsPerLID(target::Impl)::Integer
    getter and setter for Integer field numExportPacketsPerLID

numImportPacketsPerLID(target::Impl, numImportPackets::Array{<:Integer})
numImportPacketsPerLID(target::Impl)::Integer
    getter and setter for Integer field numImportPacketsPerLID


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
    # used doTransferOld since the implementation seemed the same, except for
    # extra complications converting to the new data types in doTransfer
    
    debug = false #DECISION add plist for debug setting? get debug from (im/ex)porter
    
    if !checkSizes(source, target)
        throw(InvalidArgumentError("checkSize() indicates that the destination object " *
                "is not a legal target for redistribution from the source object.  This " *
                "probably means that they do not have the same dimensions.  For example, " *
                "MultiVectors must have the same number of rows and columns"))
    end
    
    readAlso = true #from TPetras rwo
    if cm == CombineMode.INSERT || cm == CombineMode.REPLACE
        numIDsToWrite = numSameIDs + length(permuteToLIDs) + length(remoteIDs)
        if numIDsToWrite == numMyElements(map(target))
            # overwriting all local data in the destination, so write-only suffices
            
            #TODO look at FIXME on line 503
            readAlso = false
        end
    end
    
    #TODO look at FIXME on line 514
    createViews(source)
    
    #tell target to create a view of its data
    #TODO look at FIXME on line 531
    createViewNonConst(targetreadAlso)
    
    
    if numSameIDs + length(permuteToLIDs) != 0
        copyAndPermute(source, target, numSameIDs, permuteToLIDs, permuteFromLIDs)
    end
    
    # if there is a known, constant number of packets per LID
    constantNumPackets = constantNumberOfPackets(target)
    
    
    # only need to pack comm buffers if combine mode is not ZERO
    # ZERO combine mode indicates results are the same as if all zeros were recieved
    if cm != CombineMode.ZERO
        if constantNumPackets == 0
            #TODO figure out if allocation nessacery
            numImportPacketsPerLID(target, Array{Integer, 1}(length(remoteLIDs)))
        end
        
        #begin

        exports, numPackets = packAndPrepare(source, target, exportLIDs, distor)
        if isa(numPackets, Array)
            numExportPacketsPerLID(target, numPackets)
        else
            #TODO check if need allocation for numExportPacketsPerLID if constantNumPackets != 0
            constantNumPackets = numPackets
        end
        
        exportsLen = length(exports)
        
        #Line 596
        #TODO finish
    end
end

"""
    createViews(obj)

doTransfer calls this on the source object, by default it does nothing, but the source object can use this as a hint to fetch data from a compute buffer on an off-CPU decice (such as GPU) into host memory
"""
function createViews(obj)
end

"""
    createViewsNonConst(obj, readAlso::Bool)

doTransfer calls this on the target object, by default it does nothing, but the source object can use this as a hint to fetch data from a compute buffer on an off-CPU decice (such as GPU) into host memory
readAlso indicates whether the doTransfer might read from the original buffer
"""
function createViewsNonConst(obj, readAlso::Bool)
end
