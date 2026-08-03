using JLD2
using Oceananigans.TimeSteppers: Clock

function check_completion(simname)
    foldername = "../scratch/turbulence-at-many-fronts/$simname"
    
    !isdir(foldername) && return "$simname: Uninitialised (no folder)"
    filenames = readdir(foldername)
    !mapreduce(a->startswith(a, "checkpoint"), |, filenames; init=false) && return "$simname: Uninitialised (no checkpoint)"
    
    checkpointnames = filter(a->startswith(a, "checkpoint"), filenames)
    prev_time = mapreduce(max, checkpointnames; init=-Inf) do checkpoint_file
        str = "simulation/model/clock"
        checkpoint_path = joinpath(foldername, checkpoint_file)
        
        jldopen(file->file[str].time, checkpoint_path)
    end
    prev_time < 0 && return "$simname: Uninitialised (no checkpoint at t=0)"
    
    !mapreduce(a->startswith(a, "INS"), |, filenames) && return "$simname: Initialised only (no output)"
    
    iterations = jldopen(joinpath(foldername, "INS.jld2")) do file
        keys(file["timeseries/t"])
    end
    ts = jldopen(joinpath(foldername, "INS.jld2")) do file
        [file["timeseries/t/$iteration"] for iteration in keys(file["timeseries/t"])]
    end
    sp = jldopen(joinpath(foldername, "INS.jld2")) do file
        file["metadata/parameters"]
    end
    "$simname: Run until ft = $(sp.f * ts[end]), αt = $(sp.α * ts[end]), $(iterations[end]) $checkpointnames"
end

function make_filename(sp, ext=nothing, pre=""; βH=sp.βH, βα=sp.βα, βB=sp.βB, βτ=sp.βτ, θτ=sp.θτ)
    strs = map([βH, βα, βB, βτ, θτ]) do β
        replace(string(β), "."=>"_")
    end
    ext = isnothing(ext) ? "" : ".$ext"
    return joinpath(pre, join(strs, "-") * ext)
end

function θτ_from_str(ip)
    ip.βτ == 0 && return "C"
    ip.θτ == 0 && return "N"
    ip.θτ == π/2 && return "E"
end

function make_preamble(jobname, T)
    """
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --gpus-per-node=1
    #SBATCH --time=$T
    #SBATCH --job-name=$jobname
    #SBATCH --output=../scratch/logs/$jobname.txt
    
    module load julia/1.12.5
    
    # Launch from scratch
    export JULIA_DEPOT_PATH=\$SCRATCH/julia-trig
    # export JULIA_CUDA_SOFT_MEMORY_LIMIT=10%
    cd ~/turbulence-at-many-fronts
    """
end

function make_body(ip, filename)
    """
    
    output_folder=\$SCRATCH/frontogenesis-mixing/$filename

    stop_time=$(ip.stop_time)
    spinup_time=$(ip.T_spinup)
    save_time=$(ip.save_time)
    
    f=$(ip.f)
    L=$(ip.L)
    H=$(ip.H)
    
    Nx=$(ip.Nx)
    Ny=$(ip.Ny)
    Nz=$(ip.Nz)

    betaLml=$(ip.βL_ml)
    betaHml=$(ip.βH_ml)

    betax=$(ip.βx)
    betay=$(ip.βy)
    
    betaB=$(ip.βB)
    betatau=$(ip.βτ)
    thetatau=$(ip.θτ)
    
    comment="$(ip.comment)"

    julia -t 8 -- src-simulation/simulation.jl \$output_folder \$stop_time \$spinup_time \$save_time \$f \$L \$H \$Nx \$Ny \$Nz \$betaLml \$betaHml \$betax \$betay \$betaB \$betatau \$thetatau \$comment
    """
end

function make_copy(filename::String)
    """
    
    mkdir \$output_folder/../$filename
    cp \$output_folder/checkpoint* \$output_folder/../$filename
    """
end

function make_copy(filenames)
    mapreduce(make_copy, *, filenames; init="")
end


function make_script(jobname, T, ip, filename, destfilenames=[])
    preamble = make_preamble(jobname, T)
    body = make_body(ip, filename)
    return preamble * body * make_copy(destfilenames)
end

function save_script(jobname, T, ip, filename; loc="", copy_to=[])
    scriptpath = joinpath(loc, jobname * ".sh")
    println(check_completion(filename))
    write(scriptpath, make_script(jobname, T, ip, filename, copy_to))
    return nothing
end

include("ensemble.jl")

# All the simulations run until wall-time is 3 hours, then checkpoint.
# Run this file again to check for completion, then submit others
T = "4:00:00"

for (k, v) in pairs(ensemble)
    println()
    println("TEST: $(v.name)")
    path = "jobs-simulation/$k"
    mkpath(path)

    for (ip, filename) in zip(v.ips, v.filenames)
        save_script(filename, T, ip, filename; loc=path)
    end

    println()
end

println("DONE!")