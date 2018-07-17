
function string(d::KernelDensityEstimate.BallTreeDensity)
  # TODO only supports single bandwidth per dimension at this point
  pts = getPoints(d)
  return "KDE:$(size(pts,2)):$(getBW(d)[:,1]):$(pts)"
end

function parsestringvector(str::AS; dlim=',') where {AS <: AbstractString}
  sstr = split(split(strip(str),'[')[end],']')[1]
  ssstr = strip.(split(sstr,dlim))
  parse.(Float64, ssstr)
end

function convert(::Type{BallTreeDensity}, str::AS) where {AS <: AbstractString}
  @assert ismatch(r"KDE:", str)
  sstr = strip.(split(str, ':'))
  N = parse(Int, sstr[2])
  bw = parsestringvector(sstr[3])
  dims = length(bw)
  ptrowstrs = split(sstr[4],';')
  @assert dims == length(ptrowstrs)
  pts = zeros(dims, N)
  for i in 1:dims
    pts[i,:] = parsestringvector(ptrowstrs[i], dlim=' ')
  end
  kde!(pts, bw)
end
# psubs = split(psubs, '[')[end]
# psubsub = split(psubs, ']')[1]
# pw = split(psubsub, ',')




function minDistGauss!(restmp::Array{Float64, 1}, bd::BallTreeDensity, dRoot::Int, atTree::BallTreeDensity, aRoot::Int)
  #@show "minDistGauss",  dRoot, aRoot
  #atCenter = center(atTree, aRoot)
  #densCenter = center(bd, dRoot)
  #bw = bwMax(bd, dRoot)
  @fastmath @inbounds begin
    restmp[1] = 0.0
    #tmp = 0.0
    for k=1:Ndim(atTree.bt)
      restmp[2] = abs( center(atTree.bt, aRoot, k) - center(bd.bt, dRoot, k) )
      restmp[2] -= rangeB(atTree.bt, aRoot, k) + rangeB(bd.bt, dRoot, k)
      restmp[2] = (restmp[2] > 0) ? restmp[2] : 0.0
      if ( bwUniform(bd) )
        restmp[1] -= (restmp[2]*restmp[2])/bwMax(bd, dRoot, k)
      else
        restmp[1] -= (restmp[2]*restmp[2])/bwMax(bd, dRoot, k) + log(bwMin(bd, dRoot, k))
      end
    end
    restmp[1] = exp(restmp[1]/2.0)
  end
  nothing
end

function maxDistGauss!(rettmp::Array{Float64, 1}, bd::BallTreeDensity, dRoot::Int, atTree::BallTreeDensity, aRoot::Int)
  #atCenter = center(atTree, aRoot)
  #densCenter = center(bd, dRoot)
  #bw = bwMin(bd, dRoot)

    #tmp = 0.0
    rettmp[1] = 0.0
    for k in 1:Ndim(atTree.bt)
      @fastmath begin
        #@inbounds tmp = abs( center(atTree.bt, aRoot)[k] - center(bd.bt, dRoot)[k] )
        #@inbounds tmp += rangeB(atTree.bt, aRoot)[k] + rangeB(bd.bt, dRoot)[k]
        #@inbounds tmp = abs( atTree.bt.centers[((aRoot-1)*atTree.bt.dims+k)] - bd.bt.centers[((dRoot-1)*bd.bt.dims+k)] )
        #@inbounds tmp += atTree.bt.ranges[((aRoot-1)*atTree.bt.dims+k)] + bd.bt.ranges[((dRoot-1)*bd.bt.dims+k)]
        @inbounds rettmp[2] = abs( center(atTree.bt, aRoot, k) - center(bd.bt, dRoot, k) )
        @inbounds rettmp[2] += rangeB(atTree.bt, aRoot,k) + rangeB(bd.bt, dRoot, k)
        if ( bwUniform(bd) )
            #@inbounds result = result - (tmp*tmp)/(bwMin(bd, dRoot)[k])
            @inbounds rettmp[1] -= (rettmp[2]*rettmp[2])/bwMin(bd, dRoot, k)
        else
            #@inbounds result = result - (tmp*tmp)/(bwMin(bd, dRoot)[k]) + (log(bwMax(bd.bt, Root))[k])
            @inbounds rettmp[1] -= (rettmp[2]*rettmp[2])/(bwMin(bd, dRoot, k)) + (log(bwMax(bd.bt, dRoot, k))) # TODO - Root not defined here
        end
      end
    end
    rettmp[1] = exp(rettmp[1]/2.0)
    nothing
end

#   Bounds on kernel values between points in this subtree & another
function maxDistKer!(rettmp, bd::BallTreeDensity, dRoot::Int, atTree::BallTreeDensity, aRoot::Int)
#switch(getType(atTree))
#  { case Gaussian:
   maxDistGauss!(rettmp, bd, dRoot, atTree, aRoot)
#    case Laplacian:      return maxDistLaplace(dRoot,atTree,aRoot);
#    case Epanetchnikov:  return maxDistEpanetch(dRoot,atTree,aRoot);
#  }
  nothing
end

function minDistKer!(rettmp, bd::BallTreeDensity, dRoot::Int, atTree::BallTreeDensity, aRoot::Int)
#switch(getType())
#  { case Gaussian:
   minDistGauss!(rettmp, bd, dRoot, atTree, aRoot)
#    case Laplacian:      return minDistLaplace(dRoot,atTree,aRoot);
#    case Epanetchnikov:  return minDistEpanetch(dRoot,atTree,aRoot);
#  }
  nothing
end

# need to declare these here, for kernel
type pArrHdls
    pMin::Array{Float64,1}
    pMax::Array{Float64,1}
    pAdd::Array{Float64,2}
    pErr::Array{Float64,1}
    mini::Array{Float64,1}
    maxi::Array{Float64,1}
end

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

function evalDirect(bd::BallTreeDensity, dRoot::Int, atTree::BallTreeDensity, aRoot::Int, hdl::pArrHdls)
  firstFlag = true;
  minVal=2e22;
  maxVal=0.0;
  restmp = Array{Float64,1}(2)
  d = 0.0
  for j in leafFirst(atTree.bt, aRoot):leafLast(atTree.bt, aRoot)
    for i in leafFirst(bd.bt,dRoot):leafLast(bd.bt, dRoot)
      if (bd != atTree || i != j)                              # Check leave-one-out condition;
        #d = weight(bd.bt, i) * maxDistKer(bd, i, atTree, j)   #  Do direct N^2 kernel evaluation
        maxDistKer!(restmp, bd, i, atTree, j)
        d = bd.bt.weights[i] * restmp[1]                       #  Do direct N^2 kernel evaluation
        hdl.pMin[j] = hdl.pMin[j] + d
        hdl.pMax[j] = hdl.pMax[j] + d
      end
    end
    #D = 0.0
    #D = @parallel (+) for i=leafFirst(bd,dRoot):leafLast(bd, dRoot)
    #  (bd != atTree || i != j) ? weight(bd, i) * maxDistKer(bd, i, atTree, j) : 0.0               # Check leave-one-out condition;
    #end
    #hdl.pMin[j] += D
    #hdl.pMax[j] += D

    @inbounds if (hdl.pMin[j] < minVal) minVal = hdl.pMin[j]; end  # determine min & max value in this block
    @inbounds if (hdl.pMax[j] > maxVal) maxVal = hdl.pMax[j]; end
  end
  @inbounds hdl.pMin[aRoot] = minVal; hdl.pMax[aRoot] = maxVal;
  nothing
end


# Recursively evaluate the density implied by the samples of the
# subtree (rooted at dRoot) of densTree at the locations given by
# the subtree (rooted at aRoot) of *this, to within the error
# percentage "maxErr"
function evaluate(bd::BallTreeDensity, dRoot::Int, atTree::BallTreeDensity, aRoot::Int, maxErr::Float64, hdl::pArrHdls)

  #result = 0.0
  #tmp = 0.0
  restmp = Array{Float64,1}(2)

  # find the minimum and maximum effect of these two balls on each other
  minDistKer!(restmp, bd, dRoot, atTree, aRoot)
  Kmax = restmp[1]
  maxDistKer!(restmp, bd, dRoot, atTree, aRoot)
  Kmin = restmp[1]

  total = hdl.pMin[ aRoot ];		   	        # take pmin of data below this level
  #total += hdl.pAdd[1,aRoot] - hdl.pErr[aRoot]; # add lower bound from local expansion
  total += weight(bd, dRoot)*Kmin;              # also add minimum for this block

  # if the weighted contribution of this multiply is below the
  #    threshold, no need to recurse; just treat as constant
  if ( Kmax - Kmin <= maxErr * total)                     # APPROXIMATE: PERCENT
    Kmin *= weight(bd, dRoot); Kmax *= weight(bd, dRoot);

    if (bd == atTree && aRoot==dRoot)                  # LEAVE-ONE-OUT (and same subtree)
      for k in leafFirst(atTree, aRoot):(leafLast(atTree, aRoot))
        hdl.pMin[k] += Kmin * (1.0 - weight(bd, k)/weight(bd, dRoot))   # leave our weight out of it
        hdl.pMax[k] += Kmax * (1.0 - weight(bd, k)/weight(bd, dRoot))
      end
      recurseMinMax(atTree, aRoot, hdl)
    else                                                 #     NO L-O-O => just add away
      # hdl.pAdd[1,aRoot] += (Kmin + Kmax)/2.0; hdl.pErr[aRoot] = (Kmax-Kmin)/2.0;
      # !!! Should *not* do this -- instead add to local expansion (constant term)
      for k in leafFirst(atTree, aRoot):(leafLast(atTree, aRoot))
        hdl.pMin[k] += Kmin;
        hdl.pMax[k] += Kmax;
      end

      if !(isLeaf(atTree, aRoot))
          hdl.pMin[aRoot] += Kmin; hdl.pMax[aRoot] += Kmax
      end
    end

  else
    if (Npts(bd, dRoot)*Npts(atTree, aRoot)<=DirectSize)  # DIRECT EVALUATION
        evalDirect(bd, dRoot, atTree, aRoot, hdl)
    else
        # RECURSE ON SUBTREES

        # Find the subtree in closest to the other tree's left child and do
        # that first so that the values are higher and there is a better
        # chance of being able to skip a recursion.
        close_ = closer!(atTree, left(atTree, aRoot), right(atTree, aRoot), bd, left(bd, dRoot))
        if (left(bd, dRoot) != NO_CHILD && close_ != NO_CHILD)
          evaluate(bd, left(bd, dRoot), atTree, close_, maxErr, hdl);
        end
        far   = (close_ == left(atTree, aRoot)) ? right(atTree, aRoot) : left(atTree, aRoot);
        if (left(bd, dRoot) != NO_CHILD && far != NO_CHILD)
          evaluate(bd, left(bd, dRoot), atTree, far, maxErr, hdl);
        end

        # Now the same thing for the density's right child
        close_ = closer!(atTree, left(atTree, aRoot), right(atTree, aRoot), bd, right(bd, dRoot))
        if (right(bd, dRoot) != NO_CHILD && close_ != NO_CHILD)
          evaluate(bd, right(bd, dRoot) , atTree, close_, maxErr, hdl)
        end
        far   = (close_ == left(atTree, aRoot)) ? right(atTree, aRoot) : left(atTree, aRoot);
        if (right(bd, dRoot) != NO_CHILD && far != NO_CHILD)
          evaluate(bd, right(bd, dRoot), atTree, far, maxErr, hdl);
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
    end
  end#?? this one is extra

  return Union{}
end

# Dual Tree evaluation: estimate the values at this ball tree's
# points given the other tree as the samples from a distribution.
function evaluate(bd::BallTreeDensity, locations::BallTreeDensity, p::Array{Float64,1}, maxErr::Float64)
  #println("BallTreeDensity::evaluate(bd, locations, double*, maxerr) -- starting")

  if bd.bt.dims != locations.bt.dims
    error("evaluate -- dimensions of two BallTreeDensities must match")
  end
  if (p == Union{})
    error("evaluate -- density p must exist")
  end

  # TODO -- this internal declaration is making it slow!!
  hdl = pArrHdls(zeros(2*locations.bt.num_points),
                 zeros(2*locations.bt.num_points),
                 zeros(1,0),#zeros(1,2*locations.bt.dims),
                 zeros(0),#zeros(2*locations.bt.num_points),
                 [0.], [0.])

  evaluate(bd, root(), locations, root(), 2.0*maxErr, hdl)

  # Compute & account for the kernel f'ns normalization constant
  #norm = 1.0
  #switch(getType()) {
    #case Gaussian:
    norm = (2.0*pi)^((bd.bt.dims)/2.0)
    if (bwUniform(bd))
      for i in 1:bd.bt.dims
          norm *= sqrt(bd.bandwidthMax[i])
      end
    end

    #case Laplacian: norm = pow(2, ((double)Ndim()) );
    #                if (bwUniform())
    #                  for (unsigned int i=0;i<Ndim();i++) norm *= bandwidthMax[i];
    #                break;
    #case Epanetchnikov: norm = pow(4.0/3, ((double)Ndim()) );
    #                if (bwUniform())
    #                  for (unsigned int i=0;i<Ndim();i++) norm *= bandwidthMax[i];
    #                break;
  #}
  lRoot = root()
  #@show size(hdl.pMin), size(hdl.pMax)
  if (bd == locations)    # if we need to do leave-one-out
    for j in leafFirst(locations, lRoot):(leafLast(locations, lRoot))
      #println("if bd==locations j=$(j)")
      p[getIndexOf(locations, j)] = .5*(hdl.pMin[j]+hdl.pMax[j])/norm/(1-weight(bd, j))
    end
  else
    for j in leafFirst(locations, lRoot):(leafLast(locations, lRoot))
      #println("else bd==locations j=$(j)")
      p[getIndexOf(locations, j)] = .5*(hdl.pMin[j]+hdl.pMax[j])/norm;
    end
  end
  hdl.pMin = zeros(0)
  hdl.pMax = zeros(0)
  Union{}
end


function makeDualTree(bd1::BallTreeDensity, bd2::BallTreeDensity, errTol::Float64)
    #println("makeDualTree(::BTD,::BTD) -- UNTESTED!")
    #densTree = bd1
    #atTree = bd2
    pRes = zeros(Npts(bd2))
    evaluate(bd1, bd2, pRes, errTol)
    return pRes
end

function evaluateDualTree(bd::BallTreeDensity, pos::Array{Float64,2}, lvFlag::Bool=false, errTol::Float64=1e-3)
    #dim = size(pos,1)
    if (bd.bt.dims != size(pos,1)) error("bd and pos must have the same dimension") end
    if (lvFlag)
        p = makeDualTree(bd, errTol)
    else
      posKDE = makeBallTreeDensity(pos, ones(size(pos,2))./size(pos,2));
      p = makeDualTree(bd,posKDE,errTol)
    end
    return p
end
# should make this be the Union again TODO ??
function evaluateDualTree(bd::BallTreeDensity, pos::AbstractArray{Float64,1}, lvFlag::Bool=false, errTol::Float64=1e-3)
    warn("evaluateDualTree vector evaluation API is changing for single point evaluation across multiple dimensions rather than assuming multiple points on a univariate kde.")
    pos2 = zeros(1,length(pos))
    pos2[1,:] = pos[:]
    return evaluateDualTree(bd, pos2, lvFlag, errTol)
end

function makeDualTree(bd::BallTreeDensity, errTol::Float64)
    #densTree = bd ## just a new pointer to the same data...
    pRes = zeros(Npts(bd))
    #println("makeDualTree -- leave one out, going to call evaluate(densTree")
    evaluate(bd, bd, pRes, errTol) # this seems excessive

    return pRes
end

function (bd::BallTreeDensity)(pos::Array{Float64,2}, lvFlag::Bool=false, errTol::Float64=1e-3)
  evaluateDualTree(bd, pos, lvFlag, errTol)
end

function (bd::BallTreeDensity)(pos::Array{Float64,1}, lvFlag::Bool=false, errTol::Float64=1e-3)
  evaluateDualTree(bd, reshape(pos,:,1), lvFlag, errTol)
end

function evaluateDualTree(bd::BallTreeDensity, pos::BallTreeDensity, lvFlag::Bool=false, errTol::Float64=1e-3)
    if (bd.bt.dims != pos.bt.dims) error("bd and pos must have the same dimension") end
    if (lvFlag)
        p = makeDualTree(bd, errTol)
    else
        posKDE = pos
        p = makeDualTree(bd,posKDE,errTol)
    end
    return p
end

function evalAvgLogL(bd1::BallTreeDensity, bd2::BallTreeDensity)
  L = evaluateDualTree(bd1, bd2, false) # true
  #printBallTree(bd1)
  W = getWeights(bd2)
  ind = find(L.==0.0)
  ll = nothing
  if sum(find(W[ind])) > 0
    # println("evalAvgLogL -- in if")
    ll=-Inf
  else
    # println("evalAvgLogL -- in else")
    L[ind] = 1
    ll = (log.(L)')*W
  end
  return ll
end

function evalAvgLogL(bd1::BallTreeDensity, at::Array{Float64,1})
    error("evalAvgLogL(bd1::BallTreeDensity, at::Array{Float64,1}) -- not implemented yet")
end

# estimate KL-divergence D_{KL}(p1 || p2)
function kld(p1::BallTreeDensity, p2::BallTreeDensity; method::Symbol=:direct)
  if method == :direct
    return evalAvgLogL(p1,p1) - evalAvgLogL(p2,p1)
  elseif method == :unscented
    D = Ndim(p1)
    N = Npts(p1)
    ptsE = getPoints(p1)
    ptsE = repmat(ptsE,1,2*D+1)
    bw = getBW(p1);
    for i in 1:D
      ptsE[i,(i-1)*N+(1:N)] = ptsE[i,(i-1)*N+(1:N)] + bw[i,:];
      ptsE[i,(2*i-1)*N+(1:N)] = ptsE[i,(2*i-1)*N+(1:N)] - bw[i,:];
    end;
    pE = kde!(ptsE);
    return evalAvgLogL(p1,pE) - evalAvgLogL(p2,pE);
  end
end

function entropy(bd::BallTreeDensity)
    H = -evalAvgLogL(bd,bd)
    return H[1]
end

function updateBandwidth!(bd::BallTreeDensity, bw::Array{Float64, 1})
    if (bd.multibandwidth==0)
        bd.bandwidth = bw
        bd.bandwidthMax = bd.bandwidthMin = bd.bandwidth[(bd.bt.num_points*bd.bt.dims+1):end]
    else
        error("updateBandwidth! -- multibandwidth==0 ELSE not implemented yet")
    end
end

function nLOO_LL(alpha::Float64, bd::BallTreeDensity)
  alpha = alpha.^2 # assume always a Gaussian
  updateBandwidth!(bd,bd.bandwidth*alpha)
  H = entropy(bd)
  updateBandwidth!(bd,bd.bandwidth/alpha)
  return H
end

function golden(bd::BallTreeDensity, ax::Float64, bx::Float64, cx::Float64, tol::Float64=1e-2)
#GOLDEN   Minimize function of one variable using golden section search
#
#   xmin, fmin = golden(npd, f, ax, bx, cx, tol) computes a local minimum
#   of f. xmin is the computed local minimizer of f and fmin is
#   f(xmin). xmin is computed to an relative accuracy of TOL.
#
#   The parameters ax, bx and cx must satisfy the following conditions:
#   ax < bx < cx, f(bx) < f(ax) and f(bx) < f(cx).
#
#   xmin satisfies ax < xmin < cx. golden is guaranteed to succeed if f
#   is continuous between ax and cx
#
#   Roman Geus, ETH Zuerich, 9.12.97

    C = (3.0-sqrt(5.0))/2.0;  R = 1.0-C;

    x0 = ax;  x3 = cx;
    if (abs(cx-bx) > abs(bx-ax))
      x1 = bx;
      x2 = bx + C*(cx-bx)
    else
      x2 = bx;
      x1 = bx - C*(bx-ax)
    end

    #@show x1, x2
    f1 = nLOO_LL(x1,bd)
    f2 = nLOO_LL(x2,bd)

    k = 1;
    #tic()
    while abs(x3-x0) > tol*(abs(x1)+abs(x2))
    #  fprintf(1,'k=%4d, |a-b|=%e\n', k, abs(x3-x0));
      if f2 < f1
        x0 = x1;  x1 = x2;
        x2 = R*x1 + C*x3   # x2 = x1+c*(x3-x1)
        f1 = f2;
        f2 = nLOO_LL(x2,bd)
      else
        x3 = x2; x2 = x1;
        x1 = R*x2 + C*x0   # x1 = x2+c*(x0-x2)
        f2 = f1
        f1 = nLOO_LL(x1,bd)
      end
      k += 1
    #  [x0,x1,x2,x3,f1,f2]
    end
    #println("Time for while loop $(toc())")

    if f1 < f2
      xmin = x1
      fmin = f1
    else
      xmin = x2
      fmin = f2
    end
    return xmin, fmin
end

function neighborMinMax(bd::BallTreeDensity)
    tmp = (2*bd.bt.ranges).^2
    rang = reshape(bd.bt.ranges[1:(floor(Int,end/2.0))],bd.bt.dims,bd.bt.num_points)
    maxm = sqrt(sum( (2.0*rang[:,1]).^2 ))
    ssumt = sqrt.(sum( (2.0*rang[:,1:(bd.bt.num_points-1)]).^2 ,1))
    minm = minimum(ssumt)
    minm = max(minm, 1e-6)
    return minm, maxm
end

function ksize(bd::BallTreeDensity, t::String="lcv")
    # sticking with lcv to start
    Nd = bd.bt.dims; Np = bd.bt.num_points;
    ##if (t=="lcv" || 'unif','lcvp','unifp')
    minm,maxm = neighborMinMax(bd)
    p = kde!(getPoints(bd),[(minm+maxm)/2.0],getWeights(bd))#,getType(bd));
    ks, dummy =  golden(p,2.0*minm/(minm+maxm),1.0,2.0*maxm/(minm+maxm)) #(p,nLOO_LL
    ks = ks * ( minm + maxm )/2.0
    ##end
    npd = kde!(getPoints(p),[ks],getWeights(p))
    return npd
end

function kde!(points::A, autoselect::String="lcv") where {A <: AbstractArray{Float64,2}}
  p = kde!(points, [1.0])
  #BEFORE
  # p = ksize(p, autoselect)
  #AFTER with independent weights on each dimension
  dims = size(points,1)
  bwds = zeros(dims)
  for i in 1:dims
    pp = ksize(marginal(p,[i]), autoselect)
    # pp2 = kde!(samplePts[i,:], "lcv")
    bwds[i] = getBW(pp)[1]
  end
  p = kde!(points, bwds)

  return p
end

#function kde!(points::Array{Float64,2}, autoselect::String="lcv")
#  p = kde!(points, [retrieveTCP( points)])
#  return p
#end

function kde!(points::Array{Float64,1}, autoselect::String="lcv")
  return kde!(reshape(points, 1, length(points)), autoselect)
end

function getKDERange(bd::BallTreeDensity; extend::Float64=0.1)
  rangeV = nothing
  pts = getPoints(bd)
  if (bd.bt.dims == 1) && false
    rangeV = [minimum(pts),maximum(pts)]
    dr = extend*(rangeV[2]-rangeV[1])
    rangeV[1] = rangeV[1] - dr;
    rangeV[2] = rangeV[2] + dr;
  else
    rangeV = zeros(bd.bt.dims,2)
    for i in 1:bd.bt.dims
      rangeV[i,1], rangeV[i,2] = minimum(pts[i,:]), maximum(pts[i,:])
      dr = extend*(rangeV[i,2]-rangeV[i,1])
      rangeV[i,1] = rangeV[i,1] - dr;
      rangeV[i,2] = rangeV[i,2] + dr;
    end
  end

  return rangeV
end

function getKDERange(BD::Vector{BallTreeDensity}; extend::Float64=0.1)
  rangeV = getKDERange(BD[1], extend=extend)
  for i in 2:length(BD)
    tmprv = getKDERange(BD[i], extend=extend)
    for j in 1:Ndim(BD[i])
      rangeV[j,1] = rangeV[j,1] < tmprv[j,1] ? rangeV[j,1] : tmprv[j,1]
      rangeV[j,2] = rangeV[j,2] > tmprv[j,2] ? rangeV[j,2] : tmprv[j,2]
    end
  end
  return rangeV
end

function getKDERangeLinspace(bd::BallTreeDensity; extend::Float64=0.1, N::Int=201)
  v = getKDERange(bd,extend=extend)
  return linspace(v[1],v[2],N)
end

function getKDEMax(p::BallTreeDensity;N=200)
  m = zeros(p.bt.dims)
  for i in 1:p.bt.dims
    mm = marginal(p,[i])
    rangeV = getKDERange(mm)
    X = linspace(rangeV[1],rangeV[2],N)
    yV = evaluateDualTree(mm,X)
    m[i] = X[findfirst(yV,maximum(yV))]
  end
  return m
end

function getKDEMean(p::BallTreeDensity)
  return vec(Base.mean(getPoints(p),2))
end
function getKDEfit(p::BallTreeDensity; distribution=MvNormal)
  fit(distribution, getPoints(p))
end


function intersIntgAppxIS(p::BallTreeDensity, q::BallTreeDensity;N=201)
  ndims = Ndim(p)
  xx = zeros(ndims, N)
  dx = zeros(ndims)
  LD = Array{LinSpace,1}(ndims)
  for d in 1:ndims
    LD[d] = getKDERangeLinspace(marginal(p,[d]), N=N, extend=0.3)
    dx[d] = LD[d][2]-LD[d][1]
  end
  xx[1,:] = collect(getKDERangeLinspace(marginal(p,[1]), N=N, extend=0.3))
  acc = 0.0
  if ndims == 1
    yy=evaluateDualTree(p,xx[:])
    yy= yy .* evaluateDualTree(q,xx[:])
    acc = sum(yy)*dx[1]
  elseif ndims == 2
    for i in 1:N
      for j in 1:N
        @inbounds xx[2,j] = LD[2][i]
      end
      yy=evaluateDualTree(p,xx)
      yy = yy .* evaluateDualTree(q,xx)
      acc += (dx[1]*sum(yy))*dx[2]
    end
  else
    error("Can't do higher dimensions yet")
  end
  return acc
end

function test03()
    spls = zeros(1,3)
    spls[1,:] = [0.5172, 0.7169, 0.4049]
    p1000 = kde!(spls,"lcv");
    printBallTree(p1000);

    # should get
    #D= 1
    #N= 3
    #centers= 0.5609     0.46105           0      0.4049      0.5172      0.7169
    #weights= 1     0.66667           0     0.33333     0.33333     0.33333
    #ranges= 0.156     0.05615           0           0           0           0
    #leftch= 1  3  0  3  4  5
    #rightch= 5           4           0  4294967295  4294967295  4294967295
    #highest= 5  4  0  3  4  5
    #lowest= 3  3  0  3  4  5
    #perm= 0  0  0  2  0  1
    #means= 0.54633     0.46105           0      0.4049      0.5172      0.7169
    #bw= 0.05517    0.041674           0    0.038521    0.038521    0.038521
end
function test04()
    spls = zeros(2,3)
    spls[1,:] = [0.5172, 0.7169, 0.4049]
    spls[2,:] = [0.0312, 1.0094, 2.0204]
    p1000 = kde!(spls,"lcv");
    printBallTree(p1000);
    # should get
    #BallTreeDensity
    #D= 2
    #N= 3
    #centers= 0.5609      1.0258     0.61705      0.5203           0           0      0.5172      0.0312      0.7169      1.0094      0.4049      2.0204
    #weights= 1     0.66667           0     0.33333     0.33333     0.33333
    #ranges= 0.156      0.9946     0.09985      0.4891           0           0           0           0           0           0           0           0
    #leftch= 1  3  0  3  4  5
    #rightch= 5           4           0  4294967295  4294967295  4294967295
    #highest= 5  4  0  3  4  5
    #lowest= 3  3  0  3  4  5
    #perm= 0  0  0  0  1  2
    #means= 0.54633      1.0203     0.61705      0.5203           0           0      0.5172      0.0312      0.7169      1.0094      0.4049      2.0204
    #bw= 1.0268      1.6697      1.0201      1.2494           0           0      1.0101      1.0101      1.0101      1.0101      1.0101      1.0101
end
