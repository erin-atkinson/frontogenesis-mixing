# base_state.jl
# Functions describing the initial state of the simulations

@inline u₀(x, y, z) = u_initial(y, z, sp) + 1e-8 * sp.U * randn()# + 1e-3 * sp.U * sin(6π * x / sp.Lx) * sin(6π * y / sp.Ly)
@inline v₀(x, y, z) = 1e-8 * sp.U * randn()# + 1e-3 * sp.U * cos(6π * x / sp.Lx) * cos(6π * y / sp.Ly) * sp.Lx / sp.Ly
@inline w₀(x, y, z) = 0
@inline b₀(x, y, z) = b_initial(y, z, sp)
