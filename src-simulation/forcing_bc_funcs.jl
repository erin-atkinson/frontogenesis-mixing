# forcing_bc_funcs.jl

# Background stratification

@inline g(s) = log(1 + exp(s))
@inline g′(s) = 1 / (1 + exp(-s))
@inline g′′(s) = exp(-s) / (1 + exp(-s))^2

@inline b_initial(z, sp) = -sp.λ * sp.H * sp.N² * g(-(z + sp.H) / (sp.λ * sp.H))

# We turn on the background with coriolis timescale
@inline function background_spinup(t, sp)
    return 1 - exp(-t * sp.f)
end

@inline function b_flux_func(x, y, t, sp) 
    return sp.B
end

# Wind has to avoid too much oscillation
# θ: angle relative to a down-front wind
@inline function u_flux_func(x, y, t, sp) 
    min(1, t / sp.T_spinup + 1)
    return sp.τ * turnon * sin(sp.θτ)
end

@inline function v_flux_func(x, y, t, sp) 
    min(1, t / sp.T_spinup + 1)
    return sp.τ * turnon * cos(sp.θτ)
end
