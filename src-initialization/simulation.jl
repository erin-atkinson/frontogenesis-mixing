#= simulation.jl
    Initialize a baroclinic instability simulation

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder

        [02]: Simulation stop time (s)
        [03]: Save interval (s)

        [04]: Coriolis parameter (s⁻¹)
        [05]: Thermocline deformation radius (m)
        [06]: Domain height (m)

        [07]: x grid size
        [08]: y grid size
        [09]: z grid size

        [10]: Relative mixed layer depth
        [11]: Relative baroclinic growth rate

        [12]: Relative domain extent x
        [13]: Relative domain extent y

        [14]: Comment
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

simulation_parameters = begin
    stop_time = parse(Float64, ARGS[2])
    save_time = parse(Float64, ARGS[3])

    f = parse(Float64, ARGS[4])
    L = parse(Float64, ARGS[5])
    H = parse(Float64, ARGS[6])

    Nx = parse(Int64, ARGS[7])
    Ny = parse(Int64, ARGS[8])
    Nz = parse(Int64, ARGS[9])

    βH_ml = parse(Float64, ARGS[10])
    βσ = parse(Float64, ARGS[11])
    
    βx = parse(Float64, ARGS[12])
    βy = parse(Float64, ARGS[13])
    
    (;
        stop_time, save_time,
        f, L, H,
        Nx, Ny, Nz,
        βH_ml,
        βσ,
        βx, βy,
    )
end

simulation_comment = join(ARGS[14:end], " ")
simulation_comment = join(ARGS[14:end], " ")

include("create_simulation.jl")
