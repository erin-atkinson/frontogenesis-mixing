using Oceananigans.Operators
using Oceananigans.Advection: advective_momentum_flux_Uu,
                              advective_momentum_flux_Vu,
                              advective_momentum_flux_Wu,
                              advective_momentum_flux_Uv,
                              advective_momentum_flux_Vv,
                              advective_momentum_flux_Wv,
                              advective_momentum_flux_Uw,
                              advective_momentum_flux_Vw,
                              advective_momentum_flux_Ww,
                              advective_tracer_flux_x,
                              advective_tracer_flux_y,
                              advective_tracer_flux_z

using Oceananigans.Utils: SumOfArrays
using Oceananigans.OutputWriters: AveragedSpecifiedTimes

@info "Constructing time averaged outputs..."

u, v, w = model.velocities
p = PressureField(model)
b = model.tracers.b

advection = model.advection

@inline u_sq(i, j, k, grid, u) = @inbounds u[i, j, k]^2
function ke_func(i, j, k, grid, u, v, w)
    u² = ℑxᶜᵃᵃ(i, j, k, grid, u_sq, u)
    v² = ℑyᵃᶜᵃ(i, j, k, grid, u_sq, v)
    w² = ℑzᵃᵃᶜ(i, j, k, grid, u_sq, w)

    return (u² + v² + w²) / 2
end

ke = KernelFunctionOperation{Center, Center, Center}(ke_func, grid, u, v, w)

output_fields = (; u, v, w, b, ke)

writing_times = 0:sp.save_time:sp.stop_time

simulation.output_writers[:averages] = JLD2Writer(model, output_fields; 
    filename = "$output_folder/AVG.jld2", 
    schedule = AveragedSpecifiedTimes(writing_times[2:end], sp.save_time),
    overwrite_existing = true,
    with_halos = true,
    init = init_jld2!
)

simulation.output_writers[:fields] = JLD2Writer(model, (; u, v, w, b, p); 
    filename = "$output_folder/INS.jld2", 
    schedule = SpecifiedTimes(writing_times),
    overwrite_existing = true,
    with_halos = true,
    init = init_jld2!
)
