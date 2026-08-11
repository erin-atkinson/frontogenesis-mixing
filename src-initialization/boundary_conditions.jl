# boundary_conditions.jl

# ---------------------------------------
b_bcs = FieldBoundaryConditions(;
    bottom = GradientBoundaryCondition(sp.N²),
    top = GradientBoundaryCondition(sp.N²),
)
# ---------------------------------------

boundary_conditions = (; b=b_bcs)
@info "Created boundary conditions"
