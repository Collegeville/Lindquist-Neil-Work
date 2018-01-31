

function stable_test_method(foo::UInt32)
    foo+foo
end

#Used to ensure check_method gets right method
function stable_test_method(foo::Number)
    bar = "unstable"
    if foo > 10
        bar = foo
    end
    bar
end

function unstable_variables(foo::UInt8)
    bar = 1
    if foo > 8
        bar = foo
    end
    "result:$bar"
end

function unstable_return(foo::UInt32)
    if foo > 1
        1
    else
        1.0
    end
end

function unstable_combo(foo::UInt32)
    bar = 1
    if foo > 1
        bar = foo
    end
    bar
end

@testset "Stability Analysis" begin
    #TODO test check_function
    @test Tuple{Any, StabilityReport}[((UInt32,), StabilityReport())] == check_function(stable_test_method, ((UInt32,),))

    @test Tuple{Any, StabilityReport}[((UInt32,), StabilityReport()), ((Float64,), StabilityReport(Tuple{Symbol, Type}[(:bar, Union{String, Float64})], Nullable(Union{String, Float64})))] == check_function(stable_test_method, ((UInt32,),(Float64,)))

    #check_method
    @test StabilityReport() == check_method(stable_test_method, (UInt32,))
    @test StabilityReport() == check_method(unstable_variables, (UInt8,); unstable_vars=Dict(:bar=>Integer))
    @test StabilityReport() == check_method(unstable_return, (UInt32,); unstable_return=Number)

    @test StabilityReport(Tuple{Symbol, Type}[(:bar, Union{UInt8, Int})]) == check_method(unstable_variables, (UInt8,))
    @test StabilityReport(Nullable(Union{Int, Float64})) == check_method(unstable_return, (UInt32,))

    @test StabilityReport(Tuple{Symbol, Type}[(:bar, Union{UInt32, Int})]) == check_method(unstable_combo, (UInt32,); unstable_return=Number)
    @test StabilityReport(Nullable(Union{Int, UInt32})) == check_method(unstable_combo, (UInt32,); unstable_vars=Dict(:bar=>Integer))
    @test StabilityReport(Tuple{Symbol, Type}[(:bar, Union{UInt32, Int})], Nullable(Union{Int, UInt32})) == check_method(unstable_combo, (UInt32,))

    #is_stable(::StabilityReport)
    @test is_stable(StabilityReport())
    @test is_stable(StabilityReport(Nullable{Type}()))
    @test is_stable(StabilityReport(Vector{Tuple{Symbol, Type}}(0), Nullable{Type}()))
    @test !is_stable(StabilityReport(Tuple{Symbol, Type}[(:x, Number)]))
    @test !is_stable(StabilityReport(Tuple{Symbol, Type}[(:x, Number),
                                                         (:y, Array),
                                                         (:z, AbstractArray{UInt8, 1})]))
    @test !is_stable(StabilityReport(Nullable(Number)))
    @test !is_stable(StabilityReport(Tuple{Symbol, Type}[(:x, Number)], Nullable(Number)))

    #is_stable(::Vector{StabilityReport})
    @test is_stable([StabilityReport()])
    @test is_stable([StabilityReport(), StabilityReport(), StabilityReport()])
    @test !is_stable([StabilityReport(Vector{Tuple{Symbol, Type}}(0), Nullable(Number))])
    @test !is_stable([StabilityReport(), StabilityReport(Vector{Tuple{Symbol, Type}}(0), Nullable(Number))])
    @test !is_stable([StabilityReport(Vector{Tuple{Symbol, Type}}(0), Nullable(Number)), StabilityReport()])
end
