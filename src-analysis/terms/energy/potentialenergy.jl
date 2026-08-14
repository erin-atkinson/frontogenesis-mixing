@inline function potentialenergydensity(i, j, k, grid, b)
    (z, ) = Oceananigans.Grids.node(i, j, k, grid, Nothing, Nothing, Center)
    
    return @inbounds -z * b[i, j, k]
end

function PotentialEnergyDensity(b)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(potentialenergydensity, grid, b)
end
