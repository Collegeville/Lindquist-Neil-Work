
#TODO get a good understanding of whats going on for CSRMatrix and CRSGraph
#TODO make sure to add constructor to convert from Julia's sparce matrix

mutable struct CSRMatrix{Data <: Number, GID <: Integer, PID <: Integer, LID <: Integer} <: DistRowMatrix{Data, GID, PID, LID}
    rowMap::BlockMap{GID, PID, LID}
    colMap::BlockMap{GID, PID, LID}
    
end

#TODO implement RowMatrix methods
#TODO implement DistObject methods
#TODO implement SrcDistObject methods