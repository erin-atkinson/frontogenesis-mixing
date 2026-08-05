# base_state.jl
# Functions describing the initial state of the simulations

@inline u₀(x, y, z) = 1e-8 * sp.U * randn()
@inline v₀(x, y, z) = 1e-8 * sp.U * randn()
@inline w₀(x, y, z) = 0
@inline b₀(x, y, z) = b_initial(z, sp)
