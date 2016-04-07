__precompile__(true)

module KDE

using Gadfly, Color

export
    kde!,
    getPoints,
    getBW,
    getWeights,
    marginal,
    sample,
    resample,
    evaluateDualTree,

    # product operation with Multiscale Gibbs sampling
    prodAppxMSGibbsS,

    # Gadfly plotting functions
    plotKDE,
    getKDERange,
    stackMarginals,
    vstackedPlots


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
