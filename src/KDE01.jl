function kde!(points::A, ks::Array{Float64,1}, weights::Array{Float64,1}) where {A <: AbstractArray{Float64,2}}
  Nd, Np = size(points)
  if (length(ks) == 1)
    ks = repmat(ks,Nd)
  end

  ks = ks.^2 # Guassian only at this point, taking covariance
  weights = weights./sum(weights);
  #bwsize = length(ks);

  makeBallTreeDensity(points,weights,ks)

  #if (length())
end

function kde!(points::A, ks::Array{Float64,1}) where {A <: AbstractArray{Float64,2}}
  Nd, Np = size(points)
  weights = ones(Np)
  kde!(points, ks, weights)
end

function kde!(points::Array{Float64,1}, ks::Array{Float64,1})
  Np = length(points)
  pts = zeros(1,Np)
  pts[1,:] = points
  weights = ones(Np)
  kde!(pts, ks, weights)
end

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

function getType(bd::BallTreeDensity)
    return 0
end

#s = zeros(dens.D,dens.N);
#s(:,double(dens.perm(dens.N + (1:dens.N)))+1) = dens.bandwidth(:,dens.N + (1:dens.N));
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

function marginal(bd::BallTreeDensity, ind::Array{Int,1})
  pts = getPoints(bd)
  if size(bd.bandwidth,2) > 2*bd.bt.num_points
    sig = getBW(bd)
  else
    sig = getBW(bd,[1])
  end
  wts = getWeights(bd)
  p = kde!(pts[ind,:],sig[ind], wts)
end

function randKernel(N::Int, M::Int, t::Int)
  return randn(N,M)
end

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

function rand(p::BallTreeDensity, N::Int=1)
    return KernelDensityEstimate.sample(p,N)[1]
end
