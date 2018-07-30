# this file has basic setup for running a small test of JuliaPetra code

using JuliaPetra

#program settings
const useMPI = true

const commData = Float64
const commGID = UInt64
const commPID = UInt8
const commLID = UInt32


function log(values...)
    pid == 1 && println(values...)
    #print(string(myPid(comm)) * ":" * string(values[1]) * "\n")
end


if useMPI
    comm = MPIComm(commGID, commPID, commLID)
else
    comm = SerialComm{commGID, commPID, commLID}()
end

const pid = myPid(comm)
const nProc = numProc(comm)

log(comm)

#### user code here

try
    include("FiniteElementGraph.jl")
    include("FiniteElementMatrix.jl")
    #include("PowerMethod.jl")
    include("ConjugateGradient.jl")
    include("DiagonalInversePreconditioner.jl")

    problemWidth = parse(commGID, ARGS[1])
    problemHeight = parse(commGID, ARGS[2])

    map = BlockMap(problemWidth*problemHeight, comm)

    #A = make3PointSencilCSR(map, commData)
    A = FE2dMatrix{Float64}(problemWidth, problemHeight, map)
    bData = zeros(numMyElements(map), 1)
    for i in 1:length(bData)
        gRow = gid(map, i)
        x = mod(gRow-1, problemWidth) + 1
        y = fld(gRow-1, problemWidth) + 1
        if (x == 1 || x == problemWidth) && (y == 1 || y == problemWidth)
            bData[i] = 2
        elseif (x == 1 || x == problemWidth) || (y == 1 || y == problemWidth)
            bData[i] = 1
        end
    end
    preconditioner = DiagonalInversePreconditioner(A)
    b = MultiVector(map, bData)

    log("running warmup")
    tic()
    conjugate_gradient(A, b; niters = 500)
    toq()

    log("Starting conjugate gradient")
    x = conjugate_gradient(A, b; niters = 500)

    log("Results")
    barrier(getComm(x))
    print("$(myPid(getComm(x))): x = $x\n")
    barrier(getComm(x))
    log("Residual: $(norm2(apply(b, A, x, Float64(-1), Float64(1)))[1])")
catch err
    #if pid == 1
        rethrow(err)
    #end
end
