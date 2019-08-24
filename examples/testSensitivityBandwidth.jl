using KernelDensityEstimate
using KernelDensityEstimatePlotting
using Gadfly, Cairo, Fontconfig
Gadfly.set_default_plot_size(35cm,25cm)



theta = pi/4
R = [cos(theta) sin(theta); -sin(theta) cos(theta)]
offset = 3.0
leng = 10.0

x = [offset.+leng*rand(30)' 0.1*randn(30)' 10.0.+3.0*randn(10)'; 0.1.*randn(30)' offset.+leng*rand(30)' 3.0.+3.0*randn(10)']
Rx = R*x

# setOverride(0.03)
p = kde!(Rx)

plotKDE(p) |> PDF("/tmp/test.pdf")
@async run(`evince /tmp/test.pdf`)






plot(x=Rx)

#
