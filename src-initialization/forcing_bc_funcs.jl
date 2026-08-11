# forcing_bc_funcs.jl

# Background stratification
@inline b_initial(y, z, sp) = sp.N² * z + sp.Ly * sp.M² * cos(π * y / sp.Ly)^2 / π
@inline u_initial(y, z, sp) = sp.S * 2sin(π * y / sp.Ly) * cos(π * y / sp.Ly) * (z + (sp.H + sp.H_ml) / 2)
#-sp.S * (z + (sp.H + sp.H_ml) / 2) 
@inline function spindown_func(t, sp)
    t_spindown = sp.stop_time - sp.spindown_time
    t < t_spindown && return one(t)
    return zero(t)
end
