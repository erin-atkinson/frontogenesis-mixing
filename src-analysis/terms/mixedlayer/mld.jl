using Oceananigans.Operators
using Oceananigans: location

# Integrate a field down to the mixed layer depth
@inline function ml_average_func(i, j, _, grid, ℓx, ℓy, ℓz, mld, field)
    zs = znodes(grid, ℓx, ℓy, ℓz)
    Δzs = zspacings(grid, ℓx, ℓy, ℓz)
    
    res = zero(eltype(field))
    h = @inbounds mld[i, j, 1]

    for k in axes(zs, 1)
        z = @inbounds zs[k]
        Δz = @inbounds Δzs[k]
        res += @inbounds field[i, j, k] * exp(z / h) / h * Δz
    end
    
    return res
end

function ML_Average(field, mld)
    (ℓx, ℓy, ℓz) = location(field)
    grid = field.grid
    return KernelFunctionOperation{ℓx, ℓy, Nothing}(ml_average_func, grid, ℓx(), ℓy(), ℓz(), mld, field)
end
