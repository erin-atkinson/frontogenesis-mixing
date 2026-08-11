@inline function b_forcing_func(x, y, z, t, v, sp)
    return -sp.M² * v * spindown_func(t, sp)
end
# ---------------------------------------

# ---------------------------------------
b_forcing = Forcing(b_forcing_func; field_dependencies=(:v, ), parameters=sp)
# ---------------------------------------

forcing = (;)
