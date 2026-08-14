@inline function horizontaldivergence(i, j, k, grid, u, v)
    return ∂xᶜᶜᶜ(i, j, k, grid, u) + ∂yᶜᶜᶜ(i, j, k, grid, v)
end

function HorizontalDivergence(u, v)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(horizontaldivergence, grid, u, v)
end
