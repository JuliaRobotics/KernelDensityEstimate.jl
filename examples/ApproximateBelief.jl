# approximating belief

using KernelDensityEstimate, DataFrames, Gadfly, Distributions


# moved into KernelDensityEstimate
# phic(x::Float64; μ::Float64=0.0, f::Float64=1.0) = cos((x-μ)*2*pi*f)
# phis(x::Float64; μ::Float64=0.0, f::Float64=1.0) = sin((x-μ)*2*pi*f)


ran = -1.0:0.001:5
N = 200
κ  = 4

b1 = Beta(1.5,3.0)
u1 = Uniform(2.0,4.0)

truepdf = (x) -> (pdf(b1,x) + pdf(u1,x))*0.5


pTru = kde!([rand(b1,5000);rand(u1,5000)]);
μ1 = Base.mean(getPoints(pTru),2)[1]
σ1 = Base.std(getPoints(pTru),2)[1]
fr = 1.0/(κ*σ1)

phic1 = (x) -> phic(x, μ=μ1, f=fr)
phis1 = (x) -> phis(x, μ=μ1, f=fr)


fn = DataFrame(
  x = ran,
  y = truepdf.(ran),
  Function = "True"
)


nearest = DataFrame(
  x=ran,
  y=evaluateDualTree(pTru, collect(ran)),
  Function = "Closest"
)


iter = 40
K = Vector(iter)
phi1Points = zeros(iter, 2)

for i in 1:iter
  p = kde!([rand(b1,N);rand(u1,N)]);

  K[i] = DataFrame(
    x=ran,
    y=evaluateDualTree(p, collect(ran)),
    Function = "p$(i)"
  )

  x = approxHilbertInnerProd(p, phic1, N=2000)
  y = approxHilbertInnerProd(p, phis1, N=2000)
  phi1Points[i,:] = [x;y]
end



si = DataFrame(
  x = ran,
  y = phis.(ran, μ=μ1, f=fr),
  Function = "ϕs"
)

co = DataFrame(
  x = ran,
  y = phic.(ran, μ=μ1, f=fr),
  Function = "ϕc"
)



PLL = vcat(fn, K[1],K[2], co, si)


pl = plot(PLL,
  Geom.line,
  Guide.yticks(ticks=[-0.5,0.0,0.5,1.0,1.25]),
  Guide.ylabel("density, and Φ-basis"),
  Guide.xlabel("Θ"),
  Guide.title("μ=$(round(μ1,3)), frequency=$(round(fr,3)) "),
  Coord.cartesian(ymin=-0.5),
  x=:x,
  y=:y,
  color=:Function
)

draw(PDF("approxbeliefillustration.pdf",12cm,7cm),pl)

# @async run(`evince approxbeliefillustration.pdf`)


phipts = DataFrame(
  x = phi1Points[:,1],
  y = phi1Points[:,2],
  Projections="Estimates"
)

trupts = DataFrame(
  x = [approxHilbertInnerProd(pTru, phic1, N=2000)],
  y = [approxHilbertInnerProd(pTru, phis1, N=2000)],
  Projections="Closest"
)

fpts = vcat(phipts,trupts)


pl2 = plot(fpts,
  Geom.point,
  # Guide.yticks(ticks=collect(linspace(0.34,0.38,5))),
  # Guide.xticks(ticks=collect(linspace(1.035,1.105,5))),
  x=:x,
  y=:y,
  color=:Projections,
  Guide.ylabel("Imaginary, ϕs"),
  Guide.xlabel("Real, ϕc"),
  Guide.title("μ=$(round(μ1,3)), frequency=$(round(fr,3)) ")
)


draw(PDF("approxbelieffreqproj.pdf",10cm,7cm),pl2)




#
