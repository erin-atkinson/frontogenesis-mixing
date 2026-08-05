# forcing_bc_funcs.jl

# Background stratification
@inline b_initial(z, sp) = sp.N² * z

@inline function spindown_func(t, sp)
    t_spindown = sp.stop_time / 2
    t < t_spindown && return one(t)
    return zero(t)
    return one(t) - (t - t_spindown) / t_spindown
end
