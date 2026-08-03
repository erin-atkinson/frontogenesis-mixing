#= simulation.jl
    Create a simulation of a front forced by strain flow and surface cooling

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder

        [02]: Simulation stop time (s)
        [03]: Simulation spinup time (s)
        [04]: Save interval (s)

        [05: Coriolis parameter (s⁻¹)
        [06]: Thermocline deformation radius (m)
        [07]: Domain height (m)

        [08]: x grid size
        [09]: y grid size
        [10]: z grid size

        [11]: Relative mixed layer deformation radius 
        [12]: Relative mixed layer height

        [13]: Relative domain extent x
        [14]: Relative domain extent y

        [15]: Cooling parameter
        [16]: Wind parameter
        [17]: Wind angle relative to background flow

        [18]: Comment
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

simulation_parameters = begin
    stop_time = parse(Float64, ARGS[2])
    spinup_time = parse(Float64, ARGS[3])
    save_time = parse(Float64, ARGS[4])

    f = parse(Float64, ARGS[5])
    L = parse(Float64, ARGS[6])
    H = parse(Float64, ARGS[7])

    Nx = parse(Int64, ARGS[8])
    Ny = parse(Int64, ARGS[9])
    Nz = parse(Int64, ARGS[10])

    βL_ml = parse(Float64, ARGS[11])
    βH_ml = parse(Float64, ARGS[12])
    
    βx = parse(Float64, ARGS[13])
    βy = parse(Float64, ARGS[14])

    # Background
    βB = parse(Float64, ARGS[15])
    βτ = parse(Float64, ARGS[16])
    θτ = parse(Float64, ARGS[17])

    comment = join(ARGS[18:end], " ")
    
    (;
        stop_time, spinup_time, save_time,
        f, L, H,
        Nx, Ny, Nz,
        βL_ml, βH_ml,
        βx, βy,
        βB, βτ, θτ,
        comment
    )
end

include("create_simulation.jl")
