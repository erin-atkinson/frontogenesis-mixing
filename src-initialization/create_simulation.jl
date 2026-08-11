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
    advection = Centered(; order=2),
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
wizard = TimeStepWizard(; cfl=0.5, max_Δt=0.1/sp.f, diffusive_cfl=0.5)
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

@info simulation
run!(simulation)
