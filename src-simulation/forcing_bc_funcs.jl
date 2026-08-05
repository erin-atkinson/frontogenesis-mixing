# forcing_bc_funcs.jl

# Background stratification

@inline g(s) = if s > 10
    s
elseif s < -10
    exp(s)
else
    log(1 + exp(s))
end

@inline g′(s) = 1 / (1 + exp(-s))
@inline b_initial(z, sp) = -sp.λ * sp.H_ml * sp.N² * g(-(z + sp.H_ml) / (sp.λ * sp.H_ml))

# We turn on the surface with Coriolis timescale
@inline function surface_spinup(t, sp)
    return max(1 - exp(-t * sp.f), zero(t))
end

@inline function b_flux_func(x, y, t, sp) 
    return sp.B
end

@inline function u_flux_func(x, y, t, sp) 
    return sp.τ * surface_spinup(t, sp) * sin(sp.θτ)
end

@inline function v_flux_func(x, y, t, sp) 
    return sp.τ * surface_spinup(t, sp) * cos(sp.θτ)
end

@inline function sponge_profile(z, sp)
    z > -sp.H_ml && return zero(z)
    return min((z + sp.H_ml)^2 / sp.H_ml^2, one(z))
end
