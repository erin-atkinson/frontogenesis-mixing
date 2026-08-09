# base_state.jl
# Functions describing the initial state of the simulations

@inline u₀(x, y, z) = 1e-3 * sp.U * sin(10π * x / sp.Lx) * sin(10π * y / sp.Lx)
@inline v₀(x, y, z) = 1e-3 * sp.U * cos(10π * x / sp.Lx) * cos(10π * y / sp.Lx)
@inline w₀(x, y, z) = 0
@inline b₀(x, y, z) = b_initial(z, sp)
