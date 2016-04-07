# KDE.jl
Kernel Density Estimate with product approximation using multiscale Gibbs sampling.

This implementation is based on, but with further optimizations:

    Sudderth, Erik B., et al. "Nonparametric belief propagation." Communications of the ACM 53.10 (2010): 95-103.

Installation
============

    Pkg.clone("https://github.com/dehann/KDE.jl.git")

Example
=======

Bring the module into the workspace

    using KDE

Basic one dimensional examples

    # using leave-one-out likelihood cross validation for bandwidth estimation
    p100 = kde!([randn(50);10.0+2*randn(50)])
    p2 = kde!([0.0;10.0],[1.0]) # multibandwidth still to be added
    p75 = resample(p2,75)
    plotKDE([p100;p2;p75],c=["red";"green";"blue"]) # using Gadfly under the hood

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
