export CSRMatrix

mutable struct CSRMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistRowMatrix{Data, GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::Nullable{BlockMap{GID, PID, LID}}
    
    #TODO look into how nessacery these are
    importMV::Nullable{MultiVector{Data, GID, PID, LID}}
    exportMV::Nullable{MultiVector{Data, GID, PID, LID}}

    myGraph::CRSGraph{GID, PID, LID}
    
    #TODO figure out
    localMatrix::LocalCSRMatrix{Data, LID}

    #TODO look at FIXME on line 293 - try to delete values1D?
    values1D::Array{Data, 1}
    values2D::Array{Array{Data, 1}, 1}

    #pull storageStatus and fillComplete  from graph

    #Dict keys are row indices
    #first element of each tuple is a column index
    #second element of each tuple is the matching entry
    nonlocals::Dict{GID, Array{Tuple{GID, Data}, 1}}

    function CSRMatrix{Data, GID, PID, LID}(rowMap::BlockMap{GID, PID, LID}, colMap::Nullable{BlockMap{GID, PID, LID}}, myGraph::CRSGraph{GID, PID, LID}, values1D::Array{Data, 1}, values2D::Array{Array{Data, 1}, 1}) where {Data, GID, PID, LID}
        new(rowMap,
            colMap,
            Nullable{MultiVector{Data, GID, PID, LID}}(),
            Nullable{MultiVector{Data, GID, PID, LID}}(),
            myGraph,
            values1D,
            values2D,
            Dict{GID, Array{Tuple{GID, Data}, 1}}())
    end
        
end

#### Constructors ####
#TODO implement Constructors
#TODO make sure to add constructor to convert from Julia's data types



#### Internal methods ####

#### External methods ####
#TODO implement RowMatrix methods
#TODO implement DistObject methods
#TODO implement SrcDistObject methods