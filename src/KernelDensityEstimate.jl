__precompile__(true)

module KernelDensityEstimate

# using Gadfly, Colors, Cairo, Fontconfig, Compose
using Distributions, Compat

import Distributions: sample
import Base: promote_rule, *, rand

export
    MixtureDensity,
    kde!,
    getPoints,
    getBW,
    root,
    Npts,
    Ndim,
    getWeights,
    marginal,
    sample,
    rand,
    resample,
    evaluateDualTree,
    BallTree,
    BallTreeDensity,
    getKDERange,
    getKDEMax,
    getKDEMean,
    getKDEfit,
    kld,
    evalAvgLogL,
    # approximate intersection volume
    intersIntgAppxIS,

    # product operation with Multiscale Gibbs sampling
    prodAppxMSGibbsS,

    # add * operator for kde product approximate
    *,
    VectorRange


# useful aliases
@compat VoidUnion{T}   = Union{Void, T}
@compat VectorRange{T} = Union{Vector{T},UnitRange{T}}

include("BallTree01.jl")
include("BallTreeDensity01.jl")
include("DualTree01.jl")
include("KDE01.jl")
# include("KDEPlotting01.jl")
include("MSGibbs01.jl")

"""
A Julia package for Kernel Density Estimation and approximations of their products

"""


end
