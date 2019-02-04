
using JuliaPetra

function readMM(filename, comm::Comm{GID, PID, LID}) where {GID, PID, LID}

    open(filename) do file
        headerline = readline(file)
        m = match(r"%%MatrixMarket matrix coordinate ([a-z]+) ([a-z]+)", headerline)
        if m == nothing
            error("Could not process first line")
        end
        field, symmetry = m.captures
        if field == "complex" || field == "pattern"
            error("Unsupported field type $field")
        end
        if symmetry != "general"
            error("Unsupported symmetry type $symmetry")
        end

        line = readline(file)
        while line[1] == '%'
            line = readline(file)
        end
        rows, cols, nnz = map(x->parse(GID,x),split(line))

        rowmap = BlockMap(rows, comm)
        A = CSRMatrix{Float64}(rowmap, min(cols, nnz, 8), DYNAMIC_PROFILE)

        col = GID[0]
        val = Float64[0]

        while !eof(file)
            i, j, value = split(readline(file))
            gRow = parse(GID, i)

            lRow = lid(rowmap, gRow)
            if lRow != 0
                #local value

                val[1] = parse(Float64, value)
                if val[1] != 0
                    col[1] = parse(GID, j)
                    insertGlobalValues(A, gRow, col, val)
                end
            end
        end

        fillComplete(A, rowmap, rowmap; optimizeStorage=true)

        A
    end
end
