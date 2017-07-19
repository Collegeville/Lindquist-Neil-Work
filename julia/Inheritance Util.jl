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




#The following methods are based off https://github.com/abeschneider/AbstractFields.jl/blob/master/src/AbstractFields.jl
#The AbstractFields.jl package is licensed under the MIT "Expat" License:
#
#Copyright (c) 2014: Abraham Schneider.
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


# global registry of all abstract classes defined
abstract_declarations = Dict{Symbol, Array}()

macro _abstract(sym, parents, block)
    abstract_body(sym, parents, block)
end

macro _abstract(sym, block)
    abstract_body(sym, (), block)
end

function abstract_body(sym, parents, block)
    declarations = combine_declarations(parents, block)
    
    sym_name = get_symbol_name(sym)
    
    abstract_declarations[sym_name] = declarations
    
    Expr(:abstract, esc(sym_name), abstract_declarations[sym_name])
end


macro _type(sym, parents, block)
    declarations = combine_declarations(parents, block)

    sym_name = get_symbol_name(sym)
    
    if isa(parents, Symbol)
        parent = parents
    else
        parent = eval(parents)[1]
    end
    
    Expr(:type, true, Expr(:<:, esc(sym_name), esc(parent)), Expr(:block, declarations...))
end


macro _immutable(sym, parents, block)
    declarations = combine_declarations(parents, block)

    sym_name = get_symbol_name(sym)
    
    Expr(:type, false, Expr(:<:, esc(sym_name), esc(parents)), Expr(:block, declarations...))
end


function combine_declarations(parents, block)
    if isa(parents, Symbol)
        parent_declarations = abstract_declarations[parents]
    else
        #is tuple or something
        parents = eval(parents)
        parent_declarations = reduce(vcat, [], [abstract_declarations[get_symbol_name(p)] for p in parents])
    end
    child_declarations = [Expr(:(::), get_symbol_name(var.args[1]), var.args[2])
                                for var in block.args[2:2:end]]
    return unique(vcat(parent_declarations, child_declarations))
end


"""
If the given value is a symbol, returns the symbol.
If the given value is a globalref Expr, returns the local name
"""
function get_symbol_name(sym)::Symbol
    if isa(sym, Type)
        Symbol(sym)
    elseif isa(sym, Expr) && sym.head == :globalref
        sym.args[2]
    else
        sym
    end
end