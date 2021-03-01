using KernelDensityEstimate

## in 1D
p1 = kde!(randn(1,100))

# evaluate poinst
p1([-2:0.1:2;])
# 41-element Vector{Float64}:
#  0.05066429779912228
#  0.06350081520156599
#  ...


## in 3D
p3 = kde!(randn(3,75))

v = hcat([0;0;0.0], [1;0;0.0])
# 3×2 Matrix{Float64}:
#  0.0  1.0
#  0.0  0.0
#  0.0  0.0

p3(v)
# 2-element Vector{Float64}:
#  0.06352422688815289
#  0.06171254847809395

# NOTE see #75 for progress on better dispatch API
