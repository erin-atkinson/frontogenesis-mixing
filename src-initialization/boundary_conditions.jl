# boundary_conditions.jl

# ---------------------------------------
b_bcs = FieldBoundaryConditions(;
    bottom = GradientBoundaryCondition(sp.N²),
    top = GradientBoundaryCondition(sp.N²),
)
u_bcs = FieldBoundaryConditions(;
    bottom = ValueBoundaryCondition(0.0),
)
v_bcs = FieldBoundaryConditions(;
    bottom = ValueBoundaryCondition(0.0),
)
# ---------------------------------------


boundary_conditions = (; u=u_bcs, v=v_bcs, b=b_bcs)
@info "Created boundary conditions"
