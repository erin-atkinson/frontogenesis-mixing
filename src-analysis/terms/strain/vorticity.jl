@inline function vorticity_x(i, j, k, grid, v, w)
    return ∂yᶜᶠᶠ(i, j, k, grid, w) - ∂zᶜᶠᶠ(i, j, k, grid, v)
end

@inline function vorticity_y(i, j, k, grid, u, w)
    return ∂zᶠᶜᶠ(i, j, k, grid, u) - ∂xᶠᶜᶠ(i, j, k, grid, w)
end

@inline function vorticity_z(i, j, k, grid, u, v)
    return ∂xᶠᶠᶜ(i, j, k, grid, v) - ∂yᶠᶠᶜ(i, j, k, grid, u)
end

function VorticityX(v, w)
    grid = v.grid
    loc = locationornothing((Center, Face, Face), v)
    return KernelFunctionOperation{loc...}(vorticity_x, grid, v, w)
end

function VorticityY(u, w)
    grid = u.grid
    loc = locationornothing((Face, Center, Face), u)
    return KernelFunctionOperation{loc...}(vorticity_y, grid, u, w)
end

function VorticityZ(u, v)
    grid = u.grid
    loc = locationornothing((Face, Face, Center), u)
    return KernelFunctionOperation{loc...}(vorticity_z, grid, u, v)
end
