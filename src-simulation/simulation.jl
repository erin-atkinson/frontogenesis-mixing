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

        [11]: Relative mixed layer height
        [12]: Relative mixed layer deformation radius
        [13]: Relative baroclinic growth rate

        [14]: Relative domain extent x
        [15]: Relative domain extent y

        [16]: Cooling parameter
        [17]: Wind parameter
        [18]: Wind angle relative to background flow

        [19]: Comment
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

    βH_ml = parse(Float64, ARGS[11])
    βL_ml = parse(Float64, ARGS[12])
    βσ = parse(Float64, ARGS[13])
    
    βx = parse(Float64, ARGS[14])
    βy = parse(Float64, ARGS[15])

    # Background
    βB = parse(Float64, ARGS[16])
    βτ = parse(Float64, ARGS[17])
    θτ = parse(Float64, ARGS[18])
    
    (;
        stop_time, spinup_time, save_time,
        f, L, H,
        Nx, Ny, Nz,
        βH_ml, βL_ml,
        βσ,
        βx, βy,
        βB, βτ, θτ,
    )
end

simulation_comment = join(ARGS[19:end], " ")

include("create_simulation.jl")
