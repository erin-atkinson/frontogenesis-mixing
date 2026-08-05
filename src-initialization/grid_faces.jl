# grid_faces.jl

# ---------------------------------------
# All grid faces
@inline function get_grid_faces(simulation_parameters)
    sp = simulation_parameters

    zs = (-sp.Lz, -sp.H_ml)
    # Other dimensions are uniform spacing
    xs = (-sp.Lx/2, sp.Lx/2)
    ys = (-sp.Ly/2, sp.Ly/2)
    
    (; xs, ys, zs)
end
# ---------------------------------------
