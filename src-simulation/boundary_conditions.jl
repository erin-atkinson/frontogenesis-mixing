# boundary_conditions.jl

# ---------------------------------------
# Boundary conditions have no flow outside of the domain
b_bcs = FieldBoundaryConditions(;
    bottom = GradientBoundaryCondition(sp.N²),
    top = GradientBoundaryCondition(0.0)
)
u_bcs=FieldBoundaryConditions(;
    top = FluxBoundaryCondition(u_flux_func; parameters=sp),
    bottom = ValueBoundaryCondition(0.0),
)
v_bcs=FieldBoundaryConditions(;
    top=FluxBoundaryCondition(v_flux_func; parameters=sp),
    bottom = ValueBoundaryCondition(0.0),
)
# ---------------------------------------


boundary_conditions = (; u=u_bcs, v=v_bcs, b=b_bcs)
@info "Created boundary conditions"
