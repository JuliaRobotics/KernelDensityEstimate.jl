

function kde!(points::A,
              addop=(+,),
              diffop=(-,) ) where {A <: AbstractArray{Float64,2}}
  #
  dims = size(points,1)

  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=dims ? ([ (addop[1]) for i in 1:dims]...,) : addop
  diffopT = length(diffop)!=dims ? ([ (diffop[1]) for i in 1:dims]...,) : diffop

  p = kde!(points, [1.0], addopT, diffopT)
  bwds = zeros(dims)

  # TODO convert to @threads after memory allocations are avoided
  for i in 1:dims
    # TODO implement ksize! method to avoid memory allocation with pp
    pp = ksize(marginal(p,[i]), (addopT[i],), (diffopT[i],) )
    bwds[i] = getBW(pp)[1]

    # TODO add if for circular dimensions, until KDE.golden is manifold ready
  end
  p = kde!(points, bwds, addopT,  diffopT )

  return p
end


function kde!(points::Array{Float64,1}, addop=(+,), diffop=(-,) )
  return kde!(reshape(points, 1, length(points)), addop, diffop )
end

function kde!(points::A,
              ks::Array{Float64,1},
              weights::Array{Float64,1},
              addop=(+,),
              diffop=(-,)  ) where {A <: AbstractArray{Float64,2}}
  #
  Nd, Np = size(points)
  if (length(ks) == 1)
    ks = repeat(ks,Nd)
  end

  ks = ks.^2 # Guassian only at this point, taking covariance
  weights = weights./sum(weights);
  #bwsize = length(ks);
  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=Nd ? ([ (addop[1]) for i in 1:Nd]...,) : addop
  diffopT = length(diffop)!=Nd ? ([ (diffop[1]) for i in 1:Nd]...,) : diffop
  # getMuT = length(getMu)!=Ndim ? ([ getMu[1] for i in 1:Ndim]...,) : getMu
  # getLambdaT = length(getLambda)!=Ndim ? ([ getLambda[1] for i in 1:Ndim]...,) : getLambda

  makeBallTreeDensity(points, weights, ks, GaussianKer, addopT, diffopT)

  #if (length())
end

"""
    $(SIGNATURES)

Construct a BallTreeDensity object using `points` for centers and bandwidth `ks`.
"""
function kde!(points::A, ks::Array{Float64,1}, addop=(+,), diffop=(-,)) where {A <: AbstractArray{Float64,2}}

  Nd, Np = size(points)
  weights = ones(Np)

  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=Nd ? ([ (addop[1]) for i in 1:Nd]...,) : addop
  diffopT = length(diffop)!=Nd ? ([ (diffop[1]) for i in 1:Nd]...,) : diffop
  # getMuT = length(getMu)!=Ndim ? ([ getMu[1] for i in 1:Ndim]...,) : getMu
  # getLambdaT = length(getLambda)!=Ndim ? ([ getLambda[1] for i in 1:Ndim]...,) : getLambda

  kde!(points, ks, weights, addop, diffop)
end

function kde!(points::Array{Float64,1}, ks::Array{Float64,1}, addop=(+,), diffop=(-,))
  Np = length(points)
  pts = zeros(1,Np)
  pts[1,:] = points
  weights = ones(Np)
  kde!(pts, ks, weights, addop, diffop)
end

"""
    $(SIGNATURES)

Return the points (centers) used to construct the KDE.
"""
function getPoints(bd::BallTreeDensity)
    pts=zeros(bd.bt.dims, bd.bt.num_points)
    perm = bd.bt.permutation[(bd.bt.num_points + 1):end]
    bd.bt.num_points, bd.bt.dims
    res = reshape(bd.bt.centers[(bd.bt.dims*bd.bt.num_points+1):end], bd.bt.dims, bd.bt.num_points)
    ##for i in 1:bd.bt.dims
      pts[:, perm ] = res[:,:]
    ##end
    pts = pts[:,1:bd.bt.num_points]
    return pts
end


"""
    $(SIGNATURES)

Return the bandwidths used for each kernel in the density estimate.
"""
function getBW(bd::BallTreeDensity, ind::Array{Int,1}=zeros(Int,0))
  if length(ind)==0
    ind=1:bd.bt.num_points
  end
  s = zeros(bd.bt.dims,bd.bt.num_points)
  perm = bd.bt.permutation[(bd.bt.num_points + 1):end]
  res = reshape(bd.bandwidth[(bd.bt.dims*bd.bt.num_points + 1):end],bd.bt.dims,bd.bt.num_points)
  s[:, perm ] = res[:,:]
  s = s[:,ind]
  s = sqrt.(s)  # stddev gaussian covariance
  return s
end

"""
    $(SIGNATURES)

Return the weights used for each kernel in the density estimate.
"""
function getWeights(bd::BallTreeDensity, ind::Array{Int,1}=zeros(Int,0))
  if length(ind)==0
    ind=1:bd.bt.num_points
  end
  perm = bd.bt.permutation[(bd.bt.num_points+1):end]
  wts = zeros(bd.bt.num_points)
  wts[perm] = bd.bt.weights[(bd.bt.num_points+1):end]
  wts = wts[ind]
  return wts
end

"""
    $(SIGNATURES)

Extract the marginal distribution from the given higher dimensional kernel density estimate object.
"""
function marginal(bd::BallTreeDensity, ind::Array{Int,1})
  pts = getPoints(bd)
  if size(bd.bandwidth,2) > 2*bd.bt.num_points
    sig = getBW(bd)
  else
    # TODO avoid memory allocation here
    sig = getBW(bd,[1])
  end
  wts = getWeights(bd)
  p = kde!(pts[ind,:],sig[ind], wts)
end

function randKernel(N::Int, M::Int, ::Type{KernelDensityEstimate.GaussianKer}) #t::Int)
  return randn(N,M)
end

"""
   $(SIGNATURES)

Randomly sample points from the KernelDensityEstimate object.
"""
function sample(npd::BallTreeDensity, Npts::Int)
  pts = getPoints(npd)
  points = zeros(npd.bt.dims, Npts)
  ind = zeros(Int,Npts)
  bw  = getBW(npd)
  w = getWeights(npd);
  w = cumsum(w)
  w = w./w[end]
  randnums = randKernel(npd.bt.dims, Npts, getType(npd));
  t = [sort(rand(Npts));10];
  ii = 1
  for i in 1:size(pts,2)
    while (w[i] > t[ii])
      points[:,ii] = pts[:,i] + bw[:,i].*randnums[:,ii];
      ind[ii] = i
      ii += 1
    end
  end
  return points, ind
end

function sample(npd::BallTreeDensity, Npts::Int, ind::Array{Int,1})
  pts = getPoints(npd)
  points = pts[:,ind] + getBW(npd,ind).*randKernel(npd.bt.dims, length(ind), getType(npd))
  return points, ind
end

"""
    $(SIGNATURES)

Randomly sample points from the KernelDensityEstimate object.
"""
function rand(p::BallTreeDensity, N::Int=1)
    return KernelDensityEstimate.sample(p,N)[1]
end
