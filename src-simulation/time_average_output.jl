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
c = model.tracers.c

u_surface = KernelFunctionOperation{Face, Center, Nothing}(top_slice_func, grid, u)
v_surface = KernelFunctionOperation{Center, Face, Nothing}(top_slice_func, grid, v)
p_surface = KernelFunctionOperation{Center, Center, Nothing}(top_slice_func, grid, p)
b_surface = KernelFunctionOperation{Center, Center, Nothing}(top_slice_func, grid, b)
c_surface = KernelFunctionOperation{Center, Center, Nothing}(top_slice_func, grid, c)

output_fields = (; u, v, w, b, p, c)
surface_fields = (; 
    u = u_surface,
    v = v_surface,
    p = p_surface,
    b = b_surface,
    c = c_surface,
)

writing_times_pos = filter(t-> t > prev_time, 0:sp.save_time:sp.stop_time)
writing_times_neg = filter(t-> t > prev_time, (sp.start_time:sp.save_time:0)[1:end-1])
writing_times = [writing_times_neg; writing_times_pos]

surface_writing_times_pos = filter(t-> t > prev_time, 0:(sp.save_time / 10):sp.stop_time)
surface_writing_times_neg = filter(t-> t > prev_time, (sp.start_time:(sp.save_time / 10):0)[1:end-1])
surface_writing_times = [surface_writing_times_neg; surface_writing_times_pos]

surface_symbol = Symbol(:surface, prev_iteration)
output_symbol = Symbol(:ins, prev_iteration)
checkpointer_symbol = Symbol(:checkpointer, prev_iteration)

simulation.output_writers[surface_symbol] = JLD2Writer(model, surface_fields; 
    filename = "$output_folder/SURFACE.jld2", 
    schedule = SpecifiedTimes(surface_writing_times),
    overwrite_existing = false,
    with_halos = true,
    init = init_jld2!
)

simulation.output_writers[output_symbol] = JLD2Writer(model, output_fields; 
    filename = "$output_folder/INS.jld2", 
    schedule = SpecifiedTimes(writing_times),
    overwrite_existing = false,
    with_halos = true,
    init = init_jld2!
)

simulation.output_writers[checkpointer_symbol] = Checkpointer(model;
    schedule = SpecifiedTimes(writing_times),
    dir = output_folder, 
    overwrite_existing = true,
    verbose = true,
    cleanup = true
)
