# base_state.jl
# Functions describing the initial state of the simulations
# We read in a previously made simulation and assume thermal wind balance in the mixed layer...

function with_mixed_layer_w(x, y, z, w_init, z₀, sp)
    w_value = Oceananigans.interpolate((x, y, min(z, z₀)), w_init)
    return w_value
end

function with_mixed_layer_b(x, y, z, b_init, z₀, sp)
    b_value = Oceananigans.interpolate((x, y, min(z, z₀)), b_init)
    z < z₀ && return b_value
    return b_value + sp.N²_ml * (z - z₀)
end

function with_mixed_layer_u(x, y, z, u_init, u_thermal, z₀, sp)
    u_value = Oceananigans.interpolate((x, y, min(z, z₀)), u_init)
    z < z₀ && return u_value
    
    u_thermal = Oceananigans.interpolate((x, y, z), u_init) - Oceananigans.interpolate((x, y, z₀), u_init)
    return u_value + u_thermal
end

initfilename = joinpath(replace(output_folder, "frontogenesis-mixing" => "frontogenesis-mixing-initialization"), "INS.jld2")
@info "Reading in an initial condition from $initfilename"
let initfds = FieldDataset(initfilename; backend=OnDisk()),
    u_init = initfds.u[end],
    v_init = initfds.v[end],
    w_init = initfds.w[end],
    b_init = initfds.b[end]
    
    u₀ = XFaceField(grid)
    v₀ = YFaceField(grid)
    w₀ = ZFaceField(grid)
    b₀ = CenterField(grid)

    set!(b₀, (x, y, z)->with_mixed_layer_b(x, y, z, b_init, znodes(b_init)[end], sp))
    u_thermal = Field(@at (Center, Face, Center) CumulativeIntegral(-∂y(b₀) / sp.f; dims=3))
    v_thermal = Field(@at (Center, Face, Center) CumulativeIntegral(∂x(b₀) / sp.f; dims=3))
    
    set!(u₀, (x, y, z)->with_mixed_layer_u(x, y, z, u_init, u_thermal, znodes(u_init)[end], sp))
    set!(v₀, (x, y, z)->with_mixed_layer_u(x, y, z, v_init, v_thermal, znodes(v_init)[end], sp))
    set!(w₀, (x, y, z)->with_mixed_layer_w(x, y, z, w_init, znodes(w_init)[end], sp))

    set!(model; u=u₀, v=v₀, w=w₀, b=b₀)
end
