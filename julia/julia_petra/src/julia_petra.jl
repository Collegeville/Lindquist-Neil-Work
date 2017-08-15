module julia_petra

#TODO refactor to accept AbstractArray's instead of Array's

# Internal Utilities
include("Enums.jl")
include("Error.jl")
include("Macros.jl")

include("ComputeOffsets.jl")


# Communication interface
include("Distributor.jl")
include("Directory.jl")
include("Comm.jl")

include("BlockMapData.jl")
include("BlockMap.jl")

include("DirectoryMethods.jl")
include("BasicDirectory.jl")


# Serial Communication
include("SerialDistributor.jl")
#include("SerialDirectory.jl")
include("SerialComm.jl")


# MPI Communication
include("MPIUtil.jl")
include("MPIComm.jl")
include("MPIDistributor.jl")


# Data interface
include("ImportExportData.jl")
include("Import.jl")
include("Export.jl")

include("SrcDistObject.jl")
include("DistObject.jl")


# Dense Data types
include("MultiVector.jl")


# Sparse Data types
include("RowMatrix.jl")

include("LocalCRSGraph.jl")

include("CRSGraphConstructors.jl")
include("CRSGraphInternalMethods.jl")
include("CRSGraphExternalMethods.jl")

include("CSRMatrix.jl")

include("Operator.jl")

end # module
