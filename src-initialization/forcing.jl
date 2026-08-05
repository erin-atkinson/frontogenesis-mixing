@inline function v_forcing_func(x, y, z, t, w, sp)
    return -sp.S * w * spindown_func(t, sp)
end

@inline function b_forcing_func(x, y, z, t, u, sp)
    return -sp.M² * u * spindown_func(t, sp)
end
# ---------------------------------------

# ---------------------------------------
v_forcing = Forcing(v_forcing_func; field_dependencies=(:w, ), parameters=sp)
b_forcing = Forcing(b_forcing_func; field_dependencies=(:u, ), parameters=sp)
# ---------------------------------------

forcing = (; v=v_forcing, b=b_forcing)
