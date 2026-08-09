#= simulation.jl
    Initialize a baroclinic instability simulation

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder

        [02]: Simulation stop time (s)
        [03]: Simulation spindown time (s)
        [04]: Save interval (s)

        [05]: Coriolis parameter (s⁻¹)
        [06]: Thermocline deformation radius (m)
        [07]: Domain height (m)

        [08]: x grid size
        [09]: y grid size
        [10]: z grid size

        [11]: Relative mixed layer depth
        [12]: Relative baroclinic growth rate

        [13]: Relative domain extent x
        [14]: Relative domain extent y

        [15]: Comment
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

simulation_parameters = begin
    stop_time = parse(Float64, ARGS[2])
    spindown_time = parse(Float64, ARGS[3])
    save_time = parse(Float64, ARGS[4])

    f = parse(Float64, ARGS[5])
    L = parse(Float64, ARGS[6])
    H = parse(Float64, ARGS[7])

    Nx = parse(Int64, ARGS[8])
    Ny = parse(Int64, ARGS[9])
    Nz = parse(Int64, ARGS[10])

    βH_ml = parse(Float64, ARGS[11])
    βσ = parse(Float64, ARGS[12])
    
    βx = parse(Float64, ARGS[13])
    βy = parse(Float64, ARGS[14])
    
    (;
        stop_time, spindown_time, save_time,
        f, L, H,
        Nx, Ny, Nz,
        βH_ml,
        βσ,
        βx, βy,
    )
end

simulation_comment = join(ARGS[15:end], " ")

include("create_simulation.jl")
