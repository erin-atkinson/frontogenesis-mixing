# forcing_bc_funcs.jl

# Background stratification
@inline b_initial(z, sp) = sp.N² * z

@inline function spindown_func(t, sp)
    t_spindown = sp.stop_time - sp.spindown_time
    t < t_spindown && return one(t)
    return zero(t)
end
