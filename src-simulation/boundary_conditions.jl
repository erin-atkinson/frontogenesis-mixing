# boundary_conditions.jl

initfilename = joinpath(replace(output_folder, "frontogenesis-mixing" => "frontogenesis-mixing-initialization"), "INS.jld2")
@info "Reading in boundary conditions from $initfilename"

@inline function top_slice_func(i, j, k, grid, field)
    return Oceananigans.Operators.ℑzᵃᵃᶠ(i, j, grid.Nz+1, grid, field)
end

(u_bc, v_bc, b_bc) = begin
    initfds = FieldDataset(initfilename; backend=OnDisk())
    
    u_init = initfds.u[end]
    v_init = initfds.v[end]
    b_init = initfds.b[end]
    grid_init = u_init.grid
    
    u_bc = Field(KernelFunctionOperation{Face, Center, Nothing}(top_slice_func, grid_init, u_init))
    v_bc = Field(KernelFunctionOperation{Center, Face, Nothing}(top_slice_func, grid_init, v_init))
    b_bc = Field(KernelFunctionOperation{Center, Center, Nothing}(top_slice_func, grid_init, b_init))
    
    (u_bc, v_bc, b_bc)
end

# ---------------------------------------
# Boundary conditions have no flow outside of the domain
b_bcs = FieldBoundaryConditions(;
    top = GradientBoundaryCondition(0.0),
    bottom = GradientBoundaryCondition(sp.N²_ml)
)
c_bcs = FieldBoundaryConditions(;
    top = GradientBoundaryCondition(0.0),
    bottom = GradientBoundaryCondition(sp.N²_ml)
)
u_bcs = FieldBoundaryConditions(;
    top = FluxBoundaryCondition(u_flux_func; parameters=sp),
)
v_bcs = FieldBoundaryConditions(;
    top=FluxBoundaryCondition(v_flux_func; parameters=sp),
)
# ---------------------------------------


boundary_conditions = (; u=u_bcs, v=v_bcs, b=b_bcs, c=c_bcs)
@info "Created boundary conditions"
