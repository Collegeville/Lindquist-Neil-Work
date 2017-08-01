
@createFunctionGenerator(generator,
begin
    function f1(arg1::$typ, arg2::Int)
        println("f1 Success with type ", $typ)
    end

    function f2(arg1::$typ, arg2::$typ)
        println("f2 Success with type ", $typ)
    end
end)

# have to create methods with the name's f1 and f2 to allow method_exists to do the pre-generate fail
# note the argument lists are different from the generated methods
function f1()
end
function f2()
end

@test !method_exists(f1, (String, Int))
@test !method_exists(f2, (String, String))

@test !method_exists(f1, (AbstractString, Int))
@test !method_exists(f2, (AbstractString, AbstractString))

@test !method_exists(f1, (Int, Int))
@test !method_exists(f2, (Int, Int))

@generator String

@test method_exists(f1, (String, Int64))
@test method_exists(f2, (String, String))

@test !method_exists(f1, (AbstractString, Int))
@test !method_exists(f2, (AbstractString, AbstractString))

@test !method_exists(f1, (Int, Int))
@test !method_exists(f2, (Int, Int))