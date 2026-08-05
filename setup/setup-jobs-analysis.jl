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

function make_preamble(jobname, scriptname, T)
    """
    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=192
    #SBATCH --time=$T
    #SBATCH --job-name=$jobname
    #SBATCH --output=../scratch/logs/$jobname.txt
    
    module load julia/1.12.5
    
    # Launch from scratch
    export JULIA_DEPOT_PATH=\$SCRATCH/julia-tri
    cd ~/frontogenesis-mixing
    
    PPFILE=$scriptname
    RAM=/dev/shm/$jobname
    mkdir \$RAM
    """
end

function make_body(foldername, scriptname; filename="AVG.jld2", outputfilename="$scriptname.jld2")
    """
    julia -t 24 -- src-analysis/postprocess/postprocess.jl \$SCRATCH/frontogenesis-mixing/$foldername/$filename \$PPFILE $outputfilename \$RAM/$foldername &
    """
end

function make_cleanup()
    """
    wait
    
    rmdir \$RAM
    """
end

function make_script(jobname, foldernames, scriptname, T; filename="AVG.jld2", outputfilename="$scriptname.jld2")
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, scriptname; filename, outputfilename)
    end
    return make_preamble(jobname, scriptname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, scriptname, T; loc="", filename="AVG.jld2", outputfilename="$scriptname.jld2")
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, scriptname, T; filename, outputfilename))
    return nothing
end

include("ensemble.jl")

for (k, v) in pairs(ensemble)
    println()
    println("TEST: $(v.name)")
    path = "jobs-simulation/$k"
    mkpath(path)

    save_script("$setname-MEAN", v.filenames, "MEAN", "0:30:00"; loc="jobs-analysis", outputfilename="BAR.jld2")

    println()
end
