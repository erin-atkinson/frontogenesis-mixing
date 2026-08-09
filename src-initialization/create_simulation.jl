using Oceananigans
using CUDA
using JLD2
using Printf

include("base_state.jl")
include("grid_faces.jl")
include("parameters.jl")
include("forcing_bc_funcs.jl")

const sp = create_simulation_parameters(simulation_parameters)

function init_jld2!(file, model)
    file["metadata/author"] = "Erin Atkinson"
    file["metadata/comment"] = simulation_comment
    file["metadata/parameters"] = sp
    return nothing
end

!isdir(output_folder) && mkdir(output_folder)

# Get the grid
grid_faces = get_grid_faces(sp)
@info "Created grid faces"
(xs, ys, zs) = grid_faces

grid = RectilinearGrid(GPU();
    x=xs,
    y=ys,
    z=zs,
    size=(sp.Nx, sp.Ny, sp.Nz),
    topology=(Periodic, Periodic, Bounded),
    halo=(3, 3, 3)
)
@info grid

# Forcing functions
include("forcing.jl")

# Boundary conditions
include("boundary_conditions.jl")

# Closure
include("closure.jl")

model = NonhydrostaticModel(grid;
    advection = WENO(; order=5),
    coriolis = FPlane(; sp.f),
    tracers = (:b, ),
    closure,
    buoyancy = BuoyancyTracer(),
    forcing,
    boundary_conditions,
    hydrostatic_pressure_anomaly = CenterField(grid)
)

@info model
set!(model; u=u₀, v=v₀, w=w₀, b=b₀)

# Some initial timestep...
Δt = 1e-3 / sp.f

simulation = Simulation(model; Δt, stop_time=sp.stop_time, wall_time_limit=3 * 3600)

include("time_average_output.jl")

# Variable time step
wizard = TimeStepWizard(; cfl=0.5, max_Δt=0.1/sp.f)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# Output progress
const prev_t = [0.0]
const prev_wall_time = [0.0]
function progress(sim)
    i = iteration(sim)

    t = time(sim)
    Δ_time = t - prev_t[1]
    prev_t[1] = t
    t_str = @sprintf " -- Time: %.3e" t

    wall_time = sim.run_wall_time
    Δ_wall_time = wall_time - prev_wall_time[1]
    prev_wall_time[1] = wall_time

    t_per_hour = Δ_time / (Δ_wall_time / 3600)
    tph_str = @sprintf " -- Time / wall hour: %.3e" t_per_hour

    remaining_time = wall_time / 3600 + (sim.stop_time - t) / t_per_hour
    remaining_str = @sprintf "%.1f hr / %.1f hr" (wall_time / 3600) remaining_time

    str = string("Iteration: ", i, t_str, tph_str, " -- Progress: ", remaining_str)
    
    print(rpad("\r$str", 100))

    return nothing
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# ---------------------------------------
# We set the large-scale flow to zero
u, v, w = model.velocities
b = model.tracers.b

b_ref = Field{Nothing, Nothing, Center}(grid)
set!(b_ref, z->b_initial(z, sp))

u_avg = Field(Average(u; dims=(1, 2)))
v_avg = Field(Average(v; dims=(1, 2)))
b_avg = Field(Average(b; dims=(1, 2)))

u_target = Field(u - u_avg)
v_target = Field(v - v_avg)
b_target = Field(b_ref + b - b_avg)

function halt_cascade!(sim, p)
    compute!(p.u_target)
    compute!(p.v_target)
    compute!(p.b_target)
    
    set!(p.u, p.u_target)
    set!(p.v, p.v_target)
    set!(p.b, p.b_target)

    return nothing
end

simulation.callbacks[:halt_cascade] = Callback(halt_cascade!, IterationInterval(1); 
    parameters = (; u, v, b, u_target, v_target, b_target)
)

# ---------------------------------------

@info simulation
run!(simulation)
