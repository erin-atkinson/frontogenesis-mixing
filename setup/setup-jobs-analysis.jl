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
    """
end

function make_body(foldername, scriptname; filename="INS.jld2", outputfilename="$scriptname.jld2")
    """
    julia -t 24 -- src-analysis/postprocess/postprocess.jl \$SCRATCH/frontogenesis-mixing/$foldername/$filename \$PPFILE $outputfilename
    """
end

function make_cleanup()
    """
    
    """
end

function make_script(jobname, foldernames, scriptname, T; filename="INS.jld2", outputfilename="$scriptname.jld2")
    body = mapreduce(*, foldernames) do foldername
        make_body(foldername, scriptname; filename, outputfilename)
    end
    return make_preamble(jobname, scriptname, T) * body * make_cleanup()
end

function save_script(jobname, foldernames, scriptname, T; loc="", filename="INS.jld2", outputfilename="$scriptname.jld2")
    write(joinpath(loc, jobname * ".sh"), make_script(jobname, foldernames, scriptname, T; filename, outputfilename))
    return nothing
end

include("ensemble.jl")

for (k, v) in pairs(ensemble)
    println()
    println(v.name)

    loc = "jobs-analysis"

    save_script("$(v.name)-PROFILE", v.filenames, "PROFILE", "0:30:00"; loc)
    save_script("$(v.name)-FRONTOGENESIS", v.filenames, "FRONTOGENESIS", "0:30:00"; loc)

    save_script("$(v.name)-STRAINTENSOR", v.filenames, "STRAINTENSOR", "0:30:00"; loc, filename="SURFACE.jld2")
    save_script("$(v.name)-ENERGY", v.filenames, "ENERGY", "0:30:00"; loc)

    println()
end
