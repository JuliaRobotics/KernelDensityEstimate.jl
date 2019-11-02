__precompile__(true)

module KernelDensityEstimate

using LinearAlgebra
using Statistics
using Distributions
using DocStringExtensions

import Distributions: sample
import Base: promote_rule, *, rand, string, convert

export
    # override Base
    string,
    convert,

    # new stuff
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
    minkld,
    evalAvgLogL,
    # approximate intersection volume
    intersIntgAppxIS,

    # product operation with Multiscale Gibbs sampling
    prodAppxMSGibbsS,

    # add * operator for kde product approximate
    *,
    VectorRange,
    setForceEvalDirect!


# on-manifold quick fix for avoiding computational geometry speedups, see evaluate(...)
global FORCE_EVAL_DIRECT = true

function setForceEvalDirect!(flag=false)
  global FORCE_EVAL_DIRECT = flag
  @info "KernelDensityEstimate.FORCE_EVAL_DIRECT = $FORCE_EVAL_DIRECT"
  nothing
end


VectorRange{T} = Union{Vector{T},UnitRange{T}}

include("BallTree01.jl")
include("BallTreeDensity01.jl")
include("DualTree01.jl")
include("CrossValidation.jl")
include("KDE01.jl")

include("MSGibbs01.jl")
include("StringSerialization.jl")


"""
A Julia package for Kernel Density Estimation and approximations of their products
"""
end
