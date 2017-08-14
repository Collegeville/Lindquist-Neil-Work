

@test :(Base.get(graph.plist, :debug, false)) == @macroexpand julia_petra.@debug graph
@test :(Base.get(obj.plist, :debug, false))   == @macroexpand julia_petra.@debug obj