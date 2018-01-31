
export check_function, check_method
export StabilityReport, is_stable

"""
    check_function(func, signatures; unstable_vars=Dict(), unstable_return::Type)

Check that the function is stable under each of the given signatures.

Return an array of method signature-`StabilityReport` pairs from
[`check_method`](@ref).
"""
function check_function(func, signatures; unstable_vars=Dict{Symbol, Type}(),
        unstable_return::Type=Bool)
    result = Tuple{Any, StabilityReport}[]
    for params in signatures
        push!(result, (params, check_method(func, params; unstable_vars=unstable_vars, unstable_return=unstable_return)))
    end
    result
end

"""
    check_method(func, signature; unstable_vars=Dict(), unstable_return::Type)

Create a `StabilityReport` object describing the type stability of the method.

Compute non-concrete types of variables and return value, returning them in
a [`StabilityReport`](@ref) Object

#Arguments
-`unstable_vars=Dict()`: A mapping of variables that are allowed be non-concrete
types.  `get` is called with the mapping, the variable's symbol and `Bool` to
get the variable's allowed type.
-`unstable_return::Type`: A supertype of allowed, non-concrete return types.
"""
#Based off julia's code_warntype
function check_method(func, signature; unstable_vars=Dict{Symbol, Type}(), unstable_return::Type=Bool)
    function slots_used(ci, slotnames)
        used = falses(length(slotnames))
        scan_exprs!(used, ci.code)
        return used
    end

    function scan_exprs!(used, exprs)
        for ex in exprs
            if isa(ex, Slot)
                used[ex.id] = true
            elseif isa(ex, Expr)
                scan_exprs!(used, ex.args)
            end
        end
    end

    #loop over possible methods for the given argument types
    code = code_typed(func, signature)
    if length(code) == 0
        error("No methods found for $func matching $signature")
    elseif length(code) != 1
        warn("Mutliple methods for $func matching $signature")
    end

    unstable_vars_list = Array{Tuple{Symbol, Type}, 1}(0)
    unstable_ret = Nullable{Type}()

    for (src, rettyp) in code
        #check variables
        slotnames = Base.sourceinfo_slotnames(src)
        used_slotids = slots_used(src, slotnames)

        if isa(src.slottypes, Array)
            for i = 1:length(slotnames)
                if used_slotids[i]
                    name = Symbol(slotnames[i])
                    typ = src.slottypes[i]
                    if (!isleaftype(typ) || typ == Core.Box) && !(typ <: get(unstable_vars, name, Bool))
                        push!(unstable_vars_list, (name, typ))
                    end

                    #else likely optmized out
                end
            end
        else
            warn("Can't access slot types of CodeInfo")
        end

        if (!isleaftype(rettyp) || rettyp == Core.Box) && !(rettyp <: unstable_return)
            unstable_ret = Nullable(rettyp)
        end

        #TODO check body
    end

    return StabilityReport(unstable_vars_list, unstable_ret)
end

"""
    StabilityReport()
    StabilityReport(unstable_variables::Vector{Tuple{Symbol, Type}}, unstable_return::Nullable{Type})

Holds information about the stability of a method.

If `unstable_vars` and `unstable_return` are present, sets the respective
fields.  Otherwise, creates an empty list and null value respectively.

See [`is_stable`](@ref)
"""
struct StabilityReport
    "A list of unstable variables and their values"
    unstable_variables::Array{Tuple{Symbol, Type}, 1}
    "The return type, if not concrete.  Otherwise `Nullable()`"
    unstable_return::Nullable{<:Type}
end

StabilityReport() = StabilityReport(Array{Tuple{Symbol, Type}, 1}(0), Nullable{Type}())
StabilityReport(vars::Vector{Tuple{Symbol, Type}}) = StabilityReport(vars, Nullable{Type}())
StabilityReport(ret::Nullable{<:Type}) = StabilityReport(Vector{Tuple{Symbol, Type}}(0), ret)

function Base.:(==)(x::StabilityReport, y::StabilityReport)
    (x.unstable_variables == y.unstable_variables
        && if isnull(x.unstable_return)
            isnull(y.unstable_return)
        else
            !isnull(y.unstable_return) && get(x.unstable_return) == get(y.unstable_return)
        end)
end

"""
    is_stable(report::StabilityReport)::Bool
    is_stable(reports::AbstractArray{Tuple{Any, StabilityReport}})::Bool

Check if the given [`StabilityReport`](@ref)s don't have any unstable types.
"""
is_stable(report::StabilityReport)::Bool = length(report.unstable_variables) == 0 && isnull(report.unstable_return)
is_stable(reports::AbstractArray{StabilityReport})::Bool = all(@. is_stable(reports))
