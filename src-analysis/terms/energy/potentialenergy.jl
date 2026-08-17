@inline function potentialenergydensity(i, j, k, grid, b)
    (z, ) = Oceananigans.Grids.node(i, j, k, grid, nothing, nothing, Center())
    
    return @inbounds -z * b[i, j, k]
end

function PotentialEnergyDensity(b)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(potentialenergydensity, grid, b)
end
