#contains random utility macros

"""
Expands into get(\$graph.plist, :debug, false)
"""
macro debug(graph)
    :($(esc(:(Base.get)))($(esc(graph)).plist, :debug, false))
end