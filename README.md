# KernelDensityEstimate.jl

[![codecov.io](https://codecov.io/github/dehann/KernelDensityEstimate.jl/coverage.svg?branch=master)](https://codecov.io/github/dehann/KernelDensityEstimate.jl?branch=master)
[![KernelDensityEstimate](http://pkg.julialang.org/badges/KernelDensityEstimate_0.4.svg)](http://pkg.julialang.org/?pkg=KernelDensityEstimate&ver=0.4)
[![KernelDensityEstimate](http://pkg.julialang.org/badges/KernelDensityEstimate_0.5.svg)](http://pkg.julialang.org/?pkg=KernelDensityEstimate&ver=0.5)

Kernel Density Estimate with product approximation using multiscale Gibbs sampling.

All code is implemented in native Julia, including plotting which uses Gadfly. The main focus of this module is the ability to take the product between multiple KDEs, and makes this module unique from other KDE implementations. This package also supports n-dimensional KDEs. Please see examples below for details. The implementation is already fairly optimized from a symbolic standpoint and is based on work by:

    Sudderth, Erik B.; Ihler, Alexander, et al. "Nonparametric belief propagation." Communications of the ACM 53.10 (2010): 95-103.

The package has built in plotting functionality, using [Gadfly](https://github.com/dcjones/Gadfly.jl). This package is heavily used by [IncrementalInference](https://github.com/dehann/IncrementalInference.jl). Comments welcome.

Installation
============

    Pkg.add("KernelDensityEstimate")

Newest stuff:

    Pkg.checkout("KernelDensityEstimate")

Example
=======

Bring the module into the workspace

    using KernelDensityEstimate

Basic one dimensional examples

    # using leave-one-out likelihood cross validation for bandwidth estimation
    p100 = kde!([randn(50);10.0+2*randn(50)])
    p2 = kde!([0.0;10.0],[1.0]) # multibandwidth still to be added
    p75 = resample(p2,75)
    plotKDE([p100;p2;p75],c=["red";"green";"blue"]) # using Gadfly under the hood

![alt tag](https://raw.githubusercontent.com/dehann/KernelDensityEstimate.jl/master/test/FirstExamplePlot.png)

Multidimensional example

    pd2 = kde!(randn(3,100));
    @time pd2 = kde!(randn(3,100)); # defaults to loocv
    pm12 = marginal(pd2,[1;2]);
    pm2 = marginal(pm12,[2]);
    plotKDE(pm2);

Multiscale Gibbs product approximation example

    p = kde!(randn(2,100))
    q = kde!(2.0+randn(2,100))
    dummy = kde!(rand(2,100),[1.0]);
    mcmciters = 5
    pGM, = prodAppxMSGibbsS(dummy, [p;q], Union{}, Union{}, mcmciters)
    pq = kde!(pGM)
    pq1 = marginal(pq,[1])
    Pl1 = plotKDE([marginal(p,[1]);marginal(q,[1]);marginal(pq,[1])],c=["red";"green";"black"])

Direct histogram of points from the product

    using Gadfly
    Pl2 = Gadfly.plot(x=pGM[1,:],y=pGM[2,:],Geom.histogram2d);
    draw(PDF("product.pdf",15cm,8cm),hstack(Pl1,Pl2))

![alt tag](https://raw.githubusercontent.com/dehann/KernelDensityEstimate.jl/master/test/product.png)

KDE product between non-gaussian distributions

    using Distributions
    p = kde!(rand(Beta(1.0,0.45),300));
    q = kde!(rand(Rayleigh(0.5),100)-0.5);
    dummy = kde!(rand(1,100),[1.0]);
    pGM, = prodAppxMSGibbsS(dummy, [p;q], Union{}, Union{}, 5)
    pq = kde!(pGM)
    plotKDE([p;q;pq],c=["red";"green";"black"])

![alt tag](https://raw.githubusercontent.com/dehann/KernelDensityEstimate.jl/master/test/RayleighBetaProduct.png)

Draw multidimensional distributions as marginalized 2D contour plots

    using KernelDensityEstimate, Gadfly
    axis=[[-5.0;5]';[-2.0;2.0]';[-10.0;10]';[-5.0;5]']
    draw(PDF("test.pdf",30cm,20cm),
     plotKDE( kde!(randn(4,200)) ) )
    draw(PNG("MultidimPlot.png",30cm,20cm),
     plotKDE( kde!(randn(4,200)), axis=axis, dims=2:4, dimLbls=["w";"x";"y";"z"]) )

![alt tag](https://raw.githubusercontent.com/dehann/KernelDensityEstimate.jl/master/test/MultidimPlot.png)
