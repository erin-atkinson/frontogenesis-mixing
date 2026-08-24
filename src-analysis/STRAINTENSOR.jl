# Terms in the strain tensor for the surface velocity
include("terms/terms.jl")

include("terms/strain/vorticity.jl")
include("terms/strain/divergence.jl")
include("terms/strain/strain.jl")

fields = (:u, :v, :b, :c)
frames = frames[1:5:end]

u = input_fields.u
v = input_fields.v
b = input_fields.b
c = input_fields.c

ζ = Field(VorticityZ(u, v))
δ = Field(HorizontalDivergence(u, v))
σₙ = Field(NormalStrain(u, v))
σₛ = Field(ShearStrain(u, v))
σ = Field(sqrt(σₙ^2 + σₛ^2))

γb = Field(StrainEfficiency(u, v, b))
γc = Field(StrainEfficiency(u, v, c))

bx = Field(@at (Center, Center, Nothing) ∂x(b))
by = Field(@at (Center, Center, Nothing) ∂y(b))
Fb = Field((bx^2 + by^2) / 2)

cx = Field(@at (Center, Center, Nothing) ∂x(c))
cy = Field(@at (Center, Center, Nothing) ∂y(c))
Fc = Field((cx^2 + cy^2) / 2)

dependency_fields = (; ζ, δ, σₙ, σₛ, σ, γb, γc, bx, by, Fb, cx, cy, Fc)
output_fields = (; ζ, δ, σ, γb, γc, Fb, Fc)
skip_update = filter(a->a ∉ fields, keys(input_fields))
