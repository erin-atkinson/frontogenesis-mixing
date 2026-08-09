using JLD2
using Oceananigans.TimeSteppers: Clock
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
    cd ~/frontogenesis-mixing
    """
end

function make_body(ip, filename)
    """
    
    output_folder=\$SCRATCH/frontogenesis-mixing-initialization/$filename

    stop_time=$(ip.stop_time)
    spindown_time=$(ip.spindown_time)
    save_time=$(ip.save_time)
    
    f=$(ip.f)
    L=$(ip.L)
    H=$(ip.H)
    
    Nx=$(ip.Nx)
    Ny=$(ip.Ny)
    Nz=$(ip.Nz)

    betaH_ml=$(ip.βH_ml)
    betasigma=$(ip.βσ)

    betax=$(ip.βx)
    betay=$(ip.βy)
    
    comment="$(ip.comment)"

    julia -t 8 -- src-initialization/simulation.jl \$output_folder \$stop_time \$spindown_time \$save_time \$f \$L \$H \$Nx \$Ny \$Nz \$betaH_ml \$betasigma \$betax \$betay \$comment
    """
end


function make_script(jobname, T, ip, filename)
    preamble = make_preamble(jobname, T)
    body = make_body(ip, filename)
    return preamble * body
end

function save_script(jobname, T, ip, filename; loc="")
    scriptpath = joinpath(loc, jobname * ".sh")
    write(scriptpath, make_script(jobname, T, ip, filename))
    return nothing
end

include("ensemble.jl")

# All the simulations run until wall-time is 3 hours, then checkpoint.
# Run this file again to check for completion, then submit others
T = "1:00:00"

for (k, v) in pairs(ensemble)
    println()
    println("TEST: $(v.name)")
    path = "jobs-initialization/$k"
    mkpath(path)

    for (ip, filename) in zip(v.ips, v.filenames)
        save_script(filename, T, merge(ip, init_dims), filename; loc=path)
    end

    println()
end

println("DONE!")
