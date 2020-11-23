

# need to declare these here, for kernel
mutable struct pArrHdls
    pMin::Array{Float64,1}
    pMax::Array{Float64,1}
    pAdd::Array{Float64,2}
    pErr::Array{Float64,1}
    mini::Array{Float64,1}
    maxi::Array{Float64,1}
    restmp::Array{Float64,1}
end

function distGauss!(restmp::Array{Float64, 1},
                    bd::BallTreeDensity,
                    dRoot::Int,
                    atTree::BallTreeDensity,
                    aRoot::Int,
                    minmaxFnc::Function,
                    minmaxFncUni::Function,
                    mainop=(-,),
                    diffop=(-,),
                    saturate::Bool=false,
                    isX86Arch::Bool= (Base.Sys.ARCH in [:x86_64;])  )
  #
  # internal helper to for IR optimization to remove

  #
  @fastmath @inbounds begin
    restmp[1] = 0.0
    #tmp = 0.0
    for k=1:Ndim(atTree.bt)
      ## TODO upgrade for more general manifolds
      restmp[2] = abs( diffop[k]( center(atTree.bt, aRoot, k), center(bd.bt, dRoot, k)) )
      restmp[2] = mainop[k](restmp[2], rangeB(atTree.bt, aRoot, k) )
      restmp[2] = mainop[k](restmp[2], rangeB(bd.bt, dRoot, k) )
      # reasoning behind saturation for minDistGauss not clear (or well documented yet)
      restmp[2] = (saturate && (restmp[2] <= 0)) ? 0.0 : restmp[2]
      restmp[1] += (restmp[2]*restmp[2])/minmaxFnc(bd, dRoot, k)
      if  !bwUniform(bd)
        restmp[1] += log(minmaxFncUni(bd, dRoot, k))
      end
    end
    restmp[1] = exp(-0.5*restmp[1])
  end
  nothing
end


# function
maxDistGauss!(rettmp::Vector{Float64},
              bd::BallTreeDensity,
              dRoot::Int,
              atTree::BallTreeDensity,
              aRoot::Int,
              addop=(+,),
              diffop=(-,) ) = distGauss!(rettmp,bd,dRoot,atTree,aRoot,bwMin,bwMax,addop,diffop)
#


# function mainop := diffop in this case, and again diffop = diffop
minDistGauss!(restmp::Array{Float64, 1},
              bd::BallTreeDensity,
              dRoot::Int,
              atTree::BallTreeDensity,
              aRoot::Int,
              addop=(+,),
              diffop=(-,) ) = distGauss!(restmp,bd,dRoot,atTree,aRoot,bwMax,bwMin,diffop,diffop,true)
#



#   Bounds on kernel values between points in this subtree & another
maxDistKer!(rettmp,
            bd::BallTreeDensity,
            dRoot::Int,
            atTree::BallTreeDensity,
            aRoot::Int,
            addop=(+,),
            diffop=(-,) ) = maxDistGauss!(rettmp, bd, dRoot, atTree, aRoot, addop, diffop)



minDistKer!(rettmp,
            bd::BallTreeDensity,
            dRoot::Int,
            atTree::BallTreeDensity,
            aRoot::Int,
            addop=(+,),
            diffop=(-,) ) = minDistGauss!(rettmp, bd, dRoot, atTree, aRoot, addop, diffop)



function pushDownLocal(atTree::BallTreeDensity, aRoot::Int, hdl::pArrHdls)
    if !(isLeaf(atTree, aRoot))
      close = atTree.left(aRoot);
      if (close != NO_CHILD)
        hdl.pAdd[1,close] += hdl.pAdd[1,aRoot]
      end
      close = right(atTree, aRoot);
      if (close != NO_CHILD)
        hdl.pAdd[1,close] += hdl.pAdd[1,aRoot]
      end
      hdl.pAdd[1,aRoot] = 0
    end
end

function pushDownAll(locations::BallTreeDensity, hdl::pArrHdls)
  for j in root():(leafFirst(locations,root())-1)
      hdl.pAdd[1,  left(locations, j)] += hdl.pAdd[1,j]
      hdl.pAdd[1, right(locations, j)] += hdl.pAdd[1,j]
      hdl.pAdd[1,j] = 0
    end
    for j in leafFirst(locations, root()):(leafLast(locations, root())+1)
      hdl.pMin[j] += hdl.pAdd[1,j] - hdl.pErr[j]
      hdl.pMax[j] += hdl.pAdd[1,j] + hdl.pErr[j]
      hdl.pAdd[1,j] = 0; hdl.pErr[j] = 0;
    end
end

function recurseMinMax(atTree::BallTreeDensity, aRoot::Int, hdl::pArrHdls)
  l = left(atTree, aRoot); r = right(atTree, aRoot);
  if !(isLeaf(atTree, l)) recurseMinMax(atTree, l, hdl) end
  if !(isLeaf(atTree, r)) recurseMinMax(atTree, r, hdl) end
  hdl.pMin[aRoot] = hdl.pMin[l]; hdl.pMax[aRoot] = hdl.pMax[l];
  if (hdl.pMin[aRoot] > hdl.pMin[r]) hdl.pMin[aRoot] = hdl.pMin[r] end
  if (hdl.pMax[aRoot] < hdl.pMax[r]) hdl.pMax[aRoot] = hdl.pMax[r] end
end

function evalDirect(bd::BallTreeDensity,
                    dRoot::Int,
                    atTree::BallTreeDensity,
                    aRoot::Int,
                    hdl::pArrHdls,
                    addop=(+,),
                    diffop=(-,)  )
  #
  # firstFlag = true;
  minVal=2e22;
  maxVal=0.0;
  # restmp = Array{Float64,1}(undef, 2)
  # d = 0.0
  for j in leafFirst(atTree.bt, aRoot):leafLast(atTree.bt, aRoot)
    for i in leafFirst(bd.bt,dRoot):leafLast(bd.bt, dRoot)
      # TODO PoC sensitivity weighting here
      if (bd != atTree || i != j)                                  #  Check leave-one-out condition;
        #d = weight(bd.bt, i) * maxDistKer(bd, i, atTree, j)       #  Do direct N^2 kernel evaluation
        maxDistKer!(hdl.restmp, bd, i, atTree, j, addop, diffop)
        # d = bd.bt.weights[i] * hdl.restmp[1]                       #  Do direct N^2 kernel evaluation
        hdl.restmp[1] *= bd.bt.weights[i]                            #  Do direct N^2 kernel evaluation
        @inbounds hdl.pMin[j] += hdl.restmp[1] # hdl.pMin[j] + d
        @inbounds hdl.pMax[j] += hdl.restmp[1] # hdl.pMax[j] + d
      end
    end

    @inbounds if (hdl.pMin[j] < minVal) minVal = hdl.pMin[j]; end  # determine min & max value in this block
    @inbounds if (hdl.pMax[j] > maxVal) maxVal = hdl.pMax[j]; end
  end
  @inbounds hdl.pMin[aRoot] = minVal
  @inbounds hdl.pMax[aRoot] = maxVal
  nothing
end

function recurseOnSubtrees(bd::BallTreeDensity,
                           dRoot::Int,
                           atTree::BallTreeDensity,
                           aRoot::Int,
                           maxErr::Float64,
                           hdl::pArrHdls,
                           addop=(+,),
                           diffop=(-,)  )
  # RECURSE ON SUBTREES

  # Find the subtree in closest to the other tree's left child and do
  # that first so that the values are higher and there is a better
  # chance of being able to skip a recursion.
  close_ = closer!(atTree, left(atTree, aRoot), right(atTree, aRoot), bd, left(bd, dRoot), addop, diffop)
  if (left(bd, dRoot) != NO_CHILD && close_ != NO_CHILD)
    evaluate(bd, left(bd, dRoot), atTree, close_, maxErr, hdl, addop, diffop);
  end
  far   = (close_ == left(atTree, aRoot)) ? right(atTree, aRoot) : left(atTree, aRoot);
  if (left(bd, dRoot) != NO_CHILD && far != NO_CHILD)
    evaluate(bd, left(bd, dRoot), atTree, far, maxErr, hdl, addop, diffop);
  end

  # Now the same thing for the density's right child
  close_ = closer!(atTree, left(atTree, aRoot), right(atTree, aRoot), bd, right(bd, dRoot), addop, diffop)
  if (right(bd, dRoot) != NO_CHILD && close_ != NO_CHILD)
    evaluate(bd, right(bd, dRoot), atTree, close_, maxErr, hdl, addop, diffop)
  end
  far   = (close_ == left(atTree, aRoot)) ? right(atTree, aRoot) : left(atTree, aRoot);
  if (right(bd, dRoot) != NO_CHILD && far != NO_CHILD)
    evaluate(bd, right(bd, dRoot), atTree, far, maxErr, hdl, addop, diffop);
  end

  # Propogate additions in children's minimum value to this node
  if !(isLeaf(atTree, aRoot))
    hdl.pMin[aRoot] = hdl.pMin[ left(atTree, aRoot) ]
    hdl.pMax[aRoot] = hdl.pMax[ left(atTree, aRoot) ]
    if (right(atTree, aRoot) != NO_CHILD)
      if (hdl.pMin[aRoot] > hdl.pMin[ right(atTree, aRoot) ])
        hdl.pMin[aRoot] = hdl.pMin[ right(atTree, aRoot) ]
      end
      if (hdl.pMax[aRoot] < hdl.pMax[ right(atTree, aRoot) ])
        hdl.pMax[aRoot] = hdl.pMax[ right(atTree, aRoot) ]
      end
    end
  end
  nothing
end

function dontRecurseSubtrees(bd::BallTreeDensity,
                             dRoot::Int,
                             atTree::BallTreeDensity,
                             aRoot::Int,
                             hdl::pArrHdls,
                             Kmin::Float64,
                             Kmax::Float64,
                             addop=(+,),
                             diffop=(-,) )
  #
  if (bd == atTree && aRoot==dRoot)            # LEAVE-ONE-OUT (and same subtree)
    for k in leafFirst(atTree, aRoot):(leafLast(atTree, aRoot))
      hdl.pMin[k] += Kmin * (1.0 - weight(bd, k)/weight(bd, dRoot))   # leave our weight out of it
      hdl.pMax[k] += Kmax * (1.0 - weight(bd, k)/weight(bd, dRoot))
    end
    recurseMinMax(atTree, aRoot, hdl)
  else                                         #     NO L-O-O => just add away
    # hdl.pAdd[1,aRoot] += (Kmin + Kmax)/2.0; hdl.pErr[aRoot] = (Kmax-Kmin)/2.0;
    # !!! Should *not* do this -- instead add to local expansion (constant term)
    for k in leafFirst(atTree, aRoot):(leafLast(atTree, aRoot))
      hdl.pMin[k] += Kmin
      hdl.pMax[k] += Kmax
    end

    if !(isLeaf(atTree, aRoot))
      hdl.pMin[aRoot] += Kmin
      hdl.pMax[aRoot] += Kmax
    end
  end
  nothing
end

# Recursively evaluate the density implied by the samples of the
# subtree (rooted at dRoot) of densTree at the locations given by
# the subtree (rooted at aRoot) of *this, to within the error
# percentage "maxErr"
function evaluate(bd::BallTreeDensity,
                  dRoot::Int,
                  atTree::BallTreeDensity,
                  aRoot::Int,
                  maxErr::Float64,
                  hdl::pArrHdls,
                  addop=(+,),
                  diffop=(-,)  )
  #
  global DirectSize
  global FORCE_EVAL_DIRECT

  ndims = bd.bt.dims
  addopT = length(addop)!=ndims ? ([ (addop[1]) for i in 1:ndims]...,) : addop
  diffopT = length(diffop)!=ndims ? ([ (diffop[1]) for i in 1:ndims]...,) : diffop


  # TODO Nixies Issue -- proper way to compute on-manifold distances between balls
  if !FORCE_EVAL_DIRECT
    restmp = hdl.restmp # Array{Float64,1}(undef, 2)

    # find the minimum and maximum effect of these two balls on each other
    minDistKer!(restmp, bd, dRoot, atTree, aRoot, addopT, diffopT)
    Kmax = restmp[1]
    maxDistKer!(restmp, bd, dRoot, atTree, aRoot, addopT, diffopT)
    Kmin = restmp[1]

    total = hdl.pMin[ aRoot ];		   	         # take pmin of data below this level
    #total += hdl.pAdd[1,aRoot] - hdl.pErr[aRoot]; # add lower bound from local expansion
    total += weight(bd, dRoot)*Kmin;               # also add minimum for this block

    # if the weighted contribution of this multiply is below the
    #    threshold, no need to recurse; just treat as constant
    # CHECK change to on-manifold evalutation
    if ( Kmax - Kmin <= maxErr * total)            # APPROXIMATE: PERCENT
      Kmin *= weight(bd, dRoot);
      Kmax *= weight(bd, dRoot);
      dontRecurseSubtrees(bd, dRoot, atTree, aRoot, hdl, Kmin, Kmax, addopT, diffopT)
    else
      # TODO this if statement call consumes a lot of memory for some reason
      if (Npts(bd, dRoot)*Npts(atTree, aRoot)<=DirectSize)  # DIRECT EVALUATION
        evalDirect(bd, dRoot, atTree, aRoot, hdl, addopT, diffopT)
      else
        recurseOnSubtrees(bd, dRoot, atTree, aRoot, maxErr, hdl, addopT, diffopT)
      end
    end
  else
    evalDirect(bd, dRoot, atTree, aRoot, hdl, addopT, diffopT)
  end

  return nothing
end

# Dual Tree evaluation: estimate the values at this ball tree's
# points given the other tree as the samples from a distribution.
function evaluate(bd::BallTreeDensity,
                  locations::BallTreeDensity,
                  p::Array{Float64,1},
                  maxErr::Float64,
                  addop=(+,),
                  diffop=(-,) )
  #
  if bd.bt.dims != locations.bt.dims
    error("evaluate -- dimensions of two BallTreeDensities must match")
  end

  # required memory allocations TODO consider moving out of loocv loop
  hdl = pArrHdls(zeros(2*locations.bt.num_points),
                 zeros(2*locations.bt.num_points),
                 zeros(1,0), #zeros(1,2*locations.bt.dims),
                 zeros(0),   #zeros(2*locations.bt.num_points),
                 [0.], [0.],
                 zeros(2))
  #
  evaluate(bd, root(), locations, root(), 2.0*maxErr, hdl, addop, diffop)

  # Gaussian kernel (add other kernel types here)
  norm = (2.0*pi)^((bd.bt.dims)/2.0)
  if (bwUniform(bd))
    @inbounds @fastmath @simd for i in 1:bd.bt.dims
        norm *= sqrt(bd.bandwidthMax[i])
    end
  end

  lRoot = root()
  if (bd == locations)    # if we need to do leave-one-out
    for j in leafFirst(locations, lRoot):(leafLast(locations, lRoot))
      p[getIndexOf(locations, j)] = 0.5*(hdl.pMin[j]+hdl.pMax[j])/norm/(1.0-weight(bd, j))
    end
  else
    for j in leafFirst(locations, lRoot):(leafLast(locations, lRoot))
      p[getIndexOf(locations, j)] = 0.5*(hdl.pMin[j]+hdl.pMax[j])/norm;
    end
  end

  hdl.pMin = zeros(0)
  hdl.pMax = zeros(0)
  nothing
end


function makeDualTree(bd1::BallTreeDensity,
                      bd2::BallTreeDensity,
                      errTol::Float64,
                      addop=(+,),
                      diffop=(-,) )
    #
    pRes = zeros(Npts(bd2))
    # NOTE should be evaluate!
    evaluate(bd1, bd2, pRes, errTol, addop, diffop)
    return pRes
end

function makeDualTree(bd::BallTreeDensity, errTol::Float64, addop=(+,), diffop=(-,))
    #densTree = bd ## just a new pointer to the same data...
    pRes = zeros(Npts(bd))
    #println("makeDualTree -- leave one out, going to call evaluate(densTree")
    evaluate(bd, bd, pRes, errTol, addop, diffop) # this seems excessive

    return pRes
end

function evaluateDualTree(bd::BallTreeDensity,
                          pos::Array{Float64,2},
                          lvFlag::Bool=false,
                          errTol::Float64=1e-3,
                          addop=(+,),
                          diffop=(-,)  )
    #
    ndims = bd.bt.dims
    addopT = length(addop)!=ndims ? ([ (addop[1]) for i in 1:ndims]...,) : addop
    diffopT = length(diffop)!=ndims ? ([ (diffop[1]) for i in 1:ndims]...,) : diffop


    if (bd.bt.dims != size(pos,1)) error("bd and pos must have the same dimension") end
    if (lvFlag)
        p = makeDualTree(bd, errTol, addopT, diffopT)
    else
      posKDE = makeBallTreeDensity(pos, ones(size(pos,2))./size(pos,2), GaussianKer, addopT, diffopT);
      p = makeDualTree(bd,posKDE,errTol, addopT, diffopT)
    end
    return p
end

# should make this be the Union again TODO ??
# TODO: not using diffop yet
function evaluateDualTree(bd::BallTreeDensity,
                          pos::AA,
                          lvFlag::Bool=false,
                          errTol::Float64=1e-3,
                          addop=(+,),
                          diffop=(-,)  ) where {AA <: AbstractArray{Float64,1}}
    @warn "evaluateDualTree vector evaluation API is changing for single point evaluation across multiple dimensions rather than assuming multiple points on a univariate kde."
    pos2 = zeros(1,length(pos))
    pos2[1,:] = pos[:]
    return evaluateDualTree(bd, pos2, lvFlag, errTol, addop, diffop)
end


function evaluateDualTree(bd::BallTreeDensity,
                          pos::BallTreeDensity,
                          lvFlag::Bool=false,
                          errTol::Float64=1e-3,
                          addop=(+,),
                          diffop=(-,)  )
    #
    if (bd.bt.dims != pos.bt.dims) error("bd and pos must have the same dimension") end
    if (lvFlag)
        p = makeDualTree(bd, errTol, addop, diffop)
    else
        p = makeDualTree(bd, pos, errTol, addop, diffop)
    end
    return p
end


"""
    $SIGNATURES

Evaluate the KDE object at given points.

> **Note**, must use Array{Float64,2} when passing in evaluation points.
"""
function (bd::BallTreeDensity)(pos::Array{Float64,2},
                               lvFlag::Bool=false,
                               errTol::Float64=1e-3,
                               addop=(+,),
                               diffop=(-,) )
  evaluateDualTree(bd, pos, lvFlag, errTol, addop, diffop)
end
function (bd::BallTreeDensity)(pos::Array{Float64,1},
                               lvFlag::Bool=false,
                               errTol::Float64=1e-3,
                               addop=(+,),
                               diffop=(-,) )
  # TODO should it not be reshape(pos,1,:) instead?
  # @warn "whoa! check reshape inside this eval balltree function"
  evaluateDualTree(bd, reshape(pos,1,:), lvFlag, errTol, addop, diffop)
end



function evalAvgLogL(bd1::BallTreeDensity,
                     bd2::BallTreeDensity,
                     addop=(+,),
                     diffop=(-,) )
  #
  L = evaluateDualTree(bd1, bd2, false, 1e-3, addop, diffop) # true
  #printBallTree(bd1)
  W = getWeights(bd2)

  # TODO convert ind to a for loop to avoid memory allocation
  ind = findall(L.==0.0)
  ll = nothing
  if sum(findall(x->x!=0, W[ind])) > 0
    # println("evalAvgLogL -- in if")
    ll=-Inf
  else
    # println("evalAvgLogL -- in else")
    L[ind] .= 1.0
    ll = (log.(L)')*W
  end
  return ll
end

function evalAvgLogL(bd1::BallTreeDensity, at::Array{Float64,1})
    error("evalAvgLogL(bd1::BallTreeDensity, at::Array{Float64,1}) -- not implemented yet")
end

# estimate KL-divergence D_{KL}(p1 || p2)
function kld(p1::BallTreeDensity,
             p2::BallTreeDensity;
             method::Symbol=:direct,
             addop=(+,),
             diffop=(-,) )
  #
  D = Ndim(p1)

  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=D ? ([ (addop[1]) for i in 1:D]...,) : addop
  diffopT = length(diffop)!=D ? ([ (diffop[1]) for i in 1:D]...,) : diffop

  if method == :direct
    return evalAvgLogL(p1,p1, addopT, diffopT) - evalAvgLogL(p2,p1, addopT, diffopT)
  elseif method == :unscented
    N = Npts(p1)
    ptsE = getPoints(p1)
    ptsE = repeat(ptsE,1,2*D+1)
    bw = getBW(p1);
    for i in 1:D
      ptsE[i,(i-1)*N+(1:N)] = ptsE[i,(i-1)*N+(1:N)] + bw[i,:];
      ptsE[i,(2*i-1)*N+(1:N)] = ptsE[i,(2*i-1)*N+(1:N)] - bw[i,:];
    end;
    pE = kde!(ptsE);
    return evalAvgLogL(p1,pE, addopT, diffopT) - evalAvgLogL(p2,pE, addopT, diffopT);
  end
end

function entropy(bd::BallTreeDensity, addop=(+,), diffop=(-,))
    H = -evalAvgLogL(bd,bd, addop, diffop)
    return H[1]
end

minkld(p::BallTreeDensity,q::BallTreeDensity) = minimum(abs.([kld(p,q);kld(q,p)]))

function getKDERange(bd::BallTreeDensity; extend::Float64=0.1, addop=(+,), diffop=(-,) )
  rangeV = nothing
  ndims = bd.bt.dims

  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=ndims ? ([ (addop[1]) for i in 1:ndims]...,) : addop
  diffopT = length(diffop)!=ndims ? ([ (diffop[1]) for i in 1:ndims]...,) : diffop

  pts = getPoints(bd)
  if false && (bd.bt.dims == 1)
    rangeV = [minimum(pts),maximum(pts)]
    dr = extend*diffopT[1](rangeV[2], rangeV[1])
    rangeV[1] = diffopT[1](rangeV[1], dr);
    rangeV[2] = addopT[1](rangeV[2], dr);
  else
    rangeV = zeros(bd.bt.dims,2)
    for i in 1:bd.bt.dims
      rangeV[i,1] = minimum(pts[i,:])
      rangeV[i,2] = maximum(pts[i,:])
      dr = extend*diffopT[i](rangeV[i,2],rangeV[i,1])
      rangeV[i,1] = diffopT[i](rangeV[i,1], dr);
      rangeV[i,2] = addopT[i](rangeV[i,2], dr);
    end
  end

  return rangeV
end

function getKDERange(BD::Vector{BallTreeDensity}; extend::Float64=0.1, addop=(+,), diffop=(-,))
  rangeV = getKDERange(BD[1], extend=extend, addop=addop, diffop=diffop )
  for i in 2:length(BD)
    tmprv = getKDERange(BD[i], extend=extend, addop=addop, diffop=diffop )
    for j in 1:Ndim(BD[i])
      rangeV[j,1] = rangeV[j,1] < tmprv[j,1] ? rangeV[j,1] : tmprv[j,1]
      rangeV[j,2] = rangeV[j,2] > tmprv[j,2] ? rangeV[j,2] : tmprv[j,2]
    end
  end
  return rangeV
end

function getKDERangeLinspace(bd::BallTreeDensity; extend::Float64=0.1, N::Int=200, addop=(+,), diffop=(-,))
  v = getKDERange(bd,extend=extend, addop=addop, diffop=diffop )
  # TODO: udpate to use on-manifold
  return range(v[1], stop=v[2], length=N)
end

function getKDEMax(p::BallTreeDensity; N=200, addop=(+,), diffop=(-,))
  m = zeros(p.bt.dims)
  for i in 1:p.bt.dims
    mm = marginal(p,[i])
    rangeV = getKDERange(mm, addop=addop, diffop=diffop )
    X = range(rangeV[1],stop=rangeV[2],length=N)
    # yV = evaluateDualTree(mm,X)
    yV = mm(reshape(collect(X), 1, N)) # TODO should allow AbstractArray
    m[i] = X[something(findfirst(isequal(maximum(yV)), yV), 0)] # findfirst(yV,maximum(yV))
  end
  return m
end

function getKDEMean(p::BallTreeDensity)
  # TODO: update for on-manifold operations
  return vec(Statistics.mean(getPoints(p),dims=2))
end
function getKDEfit(p::BallTreeDensity; distribution=MvNormal)
  # TODO: update for on-manifold operations
  fit(distribution, getPoints(p))
end


function intersIntgAppxIS(p::BallTreeDensity,
                          q::BallTreeDensity;
                          N=201,
                          addop=(+,),
                          diffop=(-,)  )

  ndims = Ndim(p)
  xx = zeros(ndims, N)
  dx = zeros(ndims)
  LD = Array{LinRange,1}(undef, ndims)
  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=ndims ? ([ (addop[1]) for i in 1:ndims]...,) : addop
  diffopT = length(diffop)!=ndims ? ([ (diffop[1]) for i in 1:ndims]...,) : diffop
  for d in 1:ndims
    LD[d] = getKDERangeLinspace(marginal(p,[d]), N=N, extend=0.3, addop=addopT, diffop=diffopT)
    dx[d] = diffopT[d](LD[d][2], LD[d][1])
  end
  xx[1,:] = collect(getKDERangeLinspace(marginal(p,[1]), N=N, extend=0.3, addop=addopT, diffop=diffopT))
  acc = 0.0
  if ndims == 1
    yy=evaluateDualTree(p,xx[:], false, 1e-3, addopT, diffopT)
    yy= yy .* evaluateDualTree(q,xx[:], false, 1e-3, addopT, diffopT)
    acc = addopT[1](acc, sum(yy)*dx[1])
  elseif ndims == 2
    for i in 1:N
      for j in 1:N
        @inbounds xx[2,j] = LD[2][i]
      end
      yy=evaluateDualTree(p, xx, false, 1e-3, addopT, diffopT)
      yy = yy .* evaluateDualTree(q, xx, false, 1e-3, addopT, diffopT)
      @warn "two dimensional intersIntgAppxIS might not be calculating correctly -- consider using AMP.mmd instead."
      acc = addopT[1](acc, (dx[1]*sum(yy))*dx[2])
    end
  else
    error("Can't do higher dimensions yet")
  end
  return acc
end
