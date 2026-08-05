# base_state.jl
# Functions describing the initial state of the simulations
# We read in a previously made simulation
function with_mixed_layer(x, y, z, c_init, sp)
    return Oceananigans.interpolate((x, y, min(z, -sp.H_ml)), c_init)
end

initfilename = joinpath(replace(output_folder, "frontogenesis-mixing" => "frontogenesis-mixing-initialization"), "INS.jld2")
@info "Reading in an initial condition from $initfilename"
(u₀, v₀, w₀, b₀) = let initfds = FieldDataset(initfilename; backend=OnDisk()),
    u_init = initfds.u[end],
    v_init = initfds.v[end],
    w_init = initfds.w[end],
    b_init = initfds.b[end]
    
    u₀ = XFaceField(grid)
    v₀ = YFaceField(grid)
    w₀ = ZFaceField(grid)
    b₀ = CenterField(grid)

    set!(u₀, (x, y, z)->with_mixed_layer(x, y, z, u_init, sp))
    set!(v₀, (x, y, z)->with_mixed_layer(x, y, z, v_init, sp))
    set!(w₀, (x, y, z)->with_mixed_layer(x, y, z, w_init, sp))
    set!(b₀, (x, y, z)->with_mixed_layer(x, y, z, b_init, sp))

    (u₀, v₀, w₀, b₀)
end
