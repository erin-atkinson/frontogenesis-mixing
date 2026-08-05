# grid_faces.jl

# ---------------------------------------
# Calculating the variable spaced grid
@inline function z_faces_func(n, a, N)
    n <= N && return n
    n > N && return n + a * (n - N)^3 / 3
end

@inline function z_faces(sp)
    N_ml = Integer(13sp.Nz / 16)
    N_th = sp.Nz - N_ml
    Δz = sp.H_ml / N_ml
    
    a = 3 * (sp.Lz / Δz - sp.Nz) / N_th^3
    
    zs = map(n->-Δz * z_faces_func(-n, a, N_ml), -sp.Nz:0)

    return zs
end
# ---------------------------------------

# ---------------------------------------
# All grid faces
@inline function get_grid_faces(simulation_parameters)
    sp = simulation_parameters
    
    # z spacing varies
    zs = z_faces(sp)
    
    # Other dimensions are uniform spacing
    xs = (-sp.Lx/2, sp.Lx/2)
    ys = (-sp.Ly/2, sp.Ly/2)
    
    (; xs, ys, zs)
end
# ---------------------------------------
