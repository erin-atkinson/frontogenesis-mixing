# boundary_condtions.jl

# ---------------------------------------
# Boundary conditions have no flow outside of the domain
b_bcs = FieldBoundaryConditions(;
    bottom=GradientBoundaryCondition(sp.N₀²),
    top=FluxBoundaryCondition(b_flux_func; parameters=sp)
)
u_bcs=FieldBoundaryConditions(;
    top=FluxBoundaryCondition(u_flux_func; parameters=sp)
)
v_bcs=FieldBoundaryConditions(;
    top=FluxBoundaryCondition(v_flux_func; parameters=sp)
)
# ---------------------------------------


boundary_conditions = (; u=u_bcs, v=v_bcs, b=b_bcs)
@info "Created boundary conditions"
