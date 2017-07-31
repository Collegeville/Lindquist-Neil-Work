
#DECISION how to deal with the subtypes having the same implementations for methods, but requiring fields
    #mandate certain fields be present
    #mandate certain getters and setters (or other "low level" operators) be present
    #have each subtype implement (violate Don't Repeat Yourself)

"""
A base type for constructing and using dense multi-vectors, vectors and matrices in parallel.

All subtypes must have the following methods, with Impl standing in for the subtype, in addition to the methods required by SrcDistObject:
doImport(target::Impl{GID, PID, LID}, source::SrcDistObject{GID, PID, LID},
        importer::Import{GID, PID, LID}, cm::CombineMode)
    - Import data into this object using an Import object ("forward mode")

doExport(target::Impl{GID, PID, LID}, source::SrcDistObject{GID, PID, LID},
        exporter::Export{GID, PID, LID}, cm::CombineMode)
    - Export data into this object using an Export object ("forward mode")

doImport(target::Impl{GID, PID, LID}, source::SrcDistObject{GID, PID, LID},
        exporter::Export{GID, PID, LID}, cm::CombineMode)
    - Import data into this object using an Export object ("reverse mode")

doExport(target::Impl{GID, PID, LID}, source::SrcDistObject{GID, PID, LID},
        importer::Import{GID, PID, LID}, cm::CombineMode)
    - Export data into this object using an Import object ("reverse mode")
"""
abstract type DistObject{GID <:Integer, PID <: Integer, LID <: Integer} <: SrcDistObject{GID, PID, LID}
end
