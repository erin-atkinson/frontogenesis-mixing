@inline function v_forcing_func(x, y, z, t, u, sp)
    return -sp.S * w * background_spinup(t, sp)
end

@inline function b_forcing_func(x, y, z, t, u, sp)
    return -sp.M² * u * background_spinup(t, sp)
end
# ---------------------------------------


# ---------------------------------------
v_forcing = Forcing(v_forcing_func; field_dependencies=(:w, ))
b_forcing = Forcing(b_forcing_func; field_dependencies=(:u, ))
# ---------------------------------------

forcing = (; v=v_forcing, b=b_forcing)
