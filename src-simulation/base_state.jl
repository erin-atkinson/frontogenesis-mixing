# base_state.jl
# Functions describing the initial state of the simulations
# We read in a previously made simulation and assume thermal wind balance in the mixed layer...

function with_mixed_layer_u(x, y, z, u_bc, by, sp)
    S = -Oceananigans.interpolate((x, y), by) / sp.f
    u₀ = Oceananigans.interpolate((x, y), u_bc)
    return S * (z - sp.H_ml) + u₀
end

function with_mixed_layer_v(x, y, z, v_bc, bx, sp)
    S = Oceananigans.interpolate((x, y), bx) / sp.f
    v₀ = Oceananigans.interpolate((x, y), v_bc)
    return S * (z - sp.H_ml) + v₀
end

function with_mixed_layer_b(x, y, z, b_bc, sp)
    b₀ = Oceananigans.interpolate((x, y), b_bc)
    return -sp.N²_ml * z^2 / 2sp.H_ml + b₀
end

begin 
    bx = Field(∂x(b_bc))
    by = Field(∂y(b_bc))

    u₀(x, y, z) = with_mixed_layer_u(x, y, z, u_bc, by, sp)
    v₀(x, y, z) = with_mixed_layer_v(x, y, z, v_bc, bx, sp)
    b₀(x, y, z) = with_mixed_layer_b(x, y, z, b_bc, sp)
    c₀(x, y, z) = with_mixed_layer_b(x, y, z, b_bc, sp)
    
    set!(model; u=u₀, v=v₀, b=b₀, c=c₀)
end
