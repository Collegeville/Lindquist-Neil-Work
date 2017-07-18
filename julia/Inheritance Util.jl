#For adding inheritance, see
#https://github.com/abeschneider/AbstractFields.jl/blob/master/src/AbstractFields.jl


doc"""
    create_function_generator(name, functions)

Creates a macro to generate instances of the given functions.
The resulting macro takes a single type, which can be accessed
with `$typ`.

This macro can be used to simulate multiple-inheritance, by using
the resulting generator on each child type.

# Examples
```jldoctest
julia> @create_function_generator generate_combine begin
           function combine(a::$typ, b::$typ)
               return a+2*b
           end
       end
@generate_combine (macro with 1 method)

julia> @generate_combine Int
combine (generic function with 1 method)

julia> @generate_combine Array{Int64, 1}
combine (generic function with 2 methods)

julia> combine(1, 3)
7

julia> combine([1, 2], [3, 4])
2-element Array{Int64,1}:
  7
 10

julia> combine('a', 'b')
ERROR: MethodError: no method matching combine(::Char, ::Char)
```
"""
macro create_function_generator(name, functions)
    return esc(quote
        macro $name(typ)
            return $(Expr(:quote, esc(functions)))
        end
    end)
end


