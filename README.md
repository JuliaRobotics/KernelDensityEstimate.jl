# KDE.jl
Kernel Density Estimate with product approximation using multiscale Gibbs sampling.

This implementation is based on, but with further optimizations:

    Sudderth, Erik B., et al. "Nonparametric belief propagation." Communications of the ACM 53.10 (2010): 95-103.

Installation
============

    Pkg.clone("https://github.com/dehann/KDE.jl.git")

Example
=======

    using KDE
    p = kde!([randn(50);10.0+2*randn()])
    plotKDE(p) # using Gadfly under the hood
