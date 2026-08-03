# grid_cells.jl

# ---------------------------------------
# Calculating the variable spaced grid
@inline function z_faces_func(n, α, N)
    n <= N && return n
    n > N && return N + (α^(n - N) - 1)
end

@inline function z_faces(sp)
    N_ml = Integer(3sp.Nz / 4) # Want this to fail if not possible
    N_th = sp.Nz - N_ml

    Δz = sp.H_ml / N_ml
    α = (1 + H / Δz - N_ml)^(1 / N_th)
    zs = map(n->-Δz * z_faces_func(-n, α, N_ml), -sp.Nx:0)

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
