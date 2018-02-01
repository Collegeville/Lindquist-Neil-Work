
export enable_inline_stability_checks, inline_stability_checks_enabled
export @stable_function, stability_warn


run_inline_stability_checks = false

"""
    enable_inline_stability_checks(::Bool)

Sets whether to run inline stability checks from [@stable_function](@ref).

If it is set to `false` (the default value), @stable_function does not perform
any type stability checks.

The value is checked when @stable_function is evaluated, so this should useually
be set at the begining of a package definition.

See [inline_stability_checks_enabled](@ref).
"""
function enable_inline_stability_checks(enabled::Bool)
    global run_inline_stability_checks
    run_inline_stability_checks = enabled
end

"""
    inline_stability_checks_enabled()::Bool

Returns whether inline stability checks are enabled.

See [enable_inline_stability_checks](@ref).
"""
function inline_stability_checks_enabled()::Bool
    run_inline_stability_checks
end


"""
    @stable_function arg_lists function_definition(s)
    @stable_function arg_lists acceptable_instability function_definitions(s)
    @stable_function arg_lists function_name
    @stable_function arg_lists acceptable_instability function_name

Checks the type stability of the function under the given argument lists.

If the second value is a function definition, that is evaluated before checking
type stability.
"""
macro stable_function(arg_lists, unstable, func)
    if unstable isa Void || unstable == :nothing
        unstable = Dict{Symbol, Type}()
    end
    if run_inline_stability_checks
        (func_names, body) = parsebody(func)
        esc(quote
            $body
            $((:(TypeStability.stability_warn($name, TypeStability.check_function($name, $arg_lists, $unstable)))
                for name in func_names)...)
        end)
    else
        esc(func)
    end
end

macro stable_function(arg_lists, func)
    esc(:(@stable_function $arg_lists nothing $func))
end

"""
    parsebody(func)

Internal method to parse the last argument of @stable_function
"""
function parsebody(func::Expr)
    if func.head == :function
        func_names = [func.args[1].args[1]]
    elseif func.head == :block
        func_names = Symbol[]
        for expr in func.args
            if expr.head == :function
                if expr.args[1].head == :where
                    push!(func_names, expr.args[1].args[1].args[1])
                else
                    push!(func_names, expr.args[1].args[1])
                end
            end
        end
        func_names = unique(func_names)
        if length(func_names) == 0
            error("Cannot find function name in $func")
        end
    else
        error("Cannot find function name in $func")
    end
    (func_names, func)
end

function parsebody(func::Symbol)
    ([func], begin end)
end

"""
    stability_warn(func_name, report)

Displays warnings about the function if any of the reports are not stable
"""
function stability_warn(func_name, reports)
    for (args, report) in reports
        if !is_stable(report)
            println(STDERR, "$func_name($(join(args, ", "))) is not stable")
            for (var, typ) in report.unstable_variables
                println(STDERR, "  $var is of type $typ")
            end
        end
    end
end
