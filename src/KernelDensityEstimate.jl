__precompile__(true)

module KernelDensityEstimate

using Gadfly, Colors, Cairo, Fontconfig

export
    kde!,
    getPoints,
    getBW,
    root,
    Npts,
    Ndim,
    getWeights,
    marginal,
    sample,
    resample,
    evaluateDualTree,
    BallTree,
    BallTreeDensity,
    getKDERange,
    getKDEMax,
    getKDEMean,
    # approximate intersection volume
    intersIntgAppxIS,

    # product operation with Multiscale Gibbs sampling
    prodAppxMSGibbsS,

    # Gadfly plotting functions
    plotKDE,
    stackMarginals,
    vstackedPlots,
    drawHorDens


include("BallTree01.jl")
include("BallTreeDensity01.jl")
include("DualTree01.jl")
include("KDE01.jl")
include("KDEPlotting01.jl")
include("MSGibbs01.jl")

"""
A Julia package for Kernel Density Estimation and approximations of their products

"""


end
