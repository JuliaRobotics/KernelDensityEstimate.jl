#not available in Julia 0.3.8, but is in 0.4dev
#@enum(KernelType,"Gaussian", "Epanetchnikov", "Laplacian")
#@enum(Gradient,WRTMean, WRTVariance, WRTWeight)

struct GaussianKer
    val::Float64
end

global DirectSize = 100;        # if N*M is less than this, just compute.

abstract type MixtureDensity end

mutable struct BallTreeDensity <: MixtureDensity
  bt::BallTree

  KernelType
  multibandwidth::Int               # flag: is bandwidth uniform?

  means::Array{Float64,1}                  # Weighted mean of points from this level down
  bandwidth::Array{Float64,1}              # Variance or other multiscale bandwidth
  bandwidthMin::Array{Float64,1}
  bandwidthMax::Array{Float64,1}           # Bounds on BW in non-uniform case

  calcStatsHandle::Function
  swapHandle::Function
end

getType(bd::BallTreeDensity) = bd.KernelType

#function Ndim(bd::BallTreeDensity)
#    return Ndim(bd.bt)
#end
Ndim(bd::BallTreeDensity) = Ndim(bd.bt)
Npts(bd::BallTreeDensity) = Npts(bd.bt)

Npts(bd::BallTreeDensity, i::Int) = Npts(bd.bt, i)

## todo make sure these two are working properly
#function center(bd::BallTreeDensity, i::Int)
#    return center(bd.bt, i)
#end
center(bd::BallTreeDensity, i::Int) = center(bd.bt, i)

#function rangeB(bd::BallTreeDensity, i::Int)
#    return rangeB(bd.bt, i)
#end
rangeB(bd::BallTreeDensity, i::Int) = rangeB(bd.bt, i)

weight(bd::BallTreeDensity, i::Int) = weight(bd.bt,i)

isLeaf(bd::BallTreeDensity, ind::Int) = isLeaf(bd.bt, ind)

validIndex(bd::BallTreeDensity, ind::Int) = validIndex(bd.bt, ind)

left(bd::BallTreeDensity, i::Int) = left(bd.bt, i)

right(bd::BallTreeDensity, i::Int) = right(bd.bt, i)

leafFirst(bd::BallTreeDensity, i::Int) = leafFirst(bd.bt, i)

leafLast(bd::BallTreeDensity, i::Int) = leafLast(bd.bt, i)

getIndexOf(bd::BallTreeDensity, i::Int) = getIndexOf(bd.bt, i)

mean(bd::BallTreeDensity, i::Int) = bd.means[((i-1)*bd.bt.dims+1):end]
mean(bd::BallTreeDensity, i::Int, k::Int) = bd.means[((i-1)*bd.bt.dims+k)]

variance(bd::BallTreeDensity, i::Int) = bd.bandwidth[((i-1)*bd.bt.dims+1):end] # !!! only works for Gaussian
variance(bd::BallTreeDensity, i::Int, k::Int) = bd.bandwidth[((i-1)*bd.bt.dims+k)]

bw(bd::BallTreeDensity, i::Int) = bd.bandwidth[((i-1)*bd.bt.dims+1):end]
bw(bd::BallTreeDensity, i::Int, k::Int) = bd.bandwidth[((i-1)*bd.bt.dims+k)]

bwMax(bd::BallTreeDensity, i::Int) = bd.bandwidthMax[((i-1)*bd.bt.dims*bd.multibandwidth+1):end]
bwMax(bd::BallTreeDensity, i::Int, k::Int) = bd.bandwidthMax[((i-1)*bd.bt.dims*bd.multibandwidth+k)]

bwMin(bd::BallTreeDensity, i::Int) = bd.bandwidthMin[((i-1)*bd.bt.dims*bd.multibandwidth+1):end]
bwMin(bd::BallTreeDensity, i::Int, k::Int) = bd.bandwidthMin[((i-1)*bd.bt.dims*bd.multibandwidth+k)]

bwUniform(bd::BallTreeDensity) = bd.multibandwidth==0

root(bd::BallTreeDensity) = root()

function buildTree!(bd::BallTreeDensity)
  #println("buildTree!(::BallTreeDensity)")
  buildTree!(bd.bt)
  nothing
end

# Swap the ith leaf with the jth leaf.
function swapDensity!(bd::BallTreeDensity, i::Int, j::Int)
  if (i==j)
      return Union{}
  end

  swapBall!(bd.bt,i,j);

  i = (i-1)*bd.bt.dims+1;
  j = (j-1)*bd.bt.dims+1;
  for k in 1:bd.bt.dims
    tmp = bd.means[i];
    bd.means[i] = bd.means[j];
    bd.means[j] = tmp;
    tmp = bd.bandwidth[i];
    bd.bandwidth[i]  = bd.bandwidth[j];
    bd.bandwidth[j]  = tmp;
    if (!bwUniform(bd))
      tmp = bd.bandwidthMax[i];bd.bandwidthMax[i]=bd.bandwidthMax[j];bd.bandwidthMax[j]=tmp;
      tmp = bd.bandwidthMin[i];bd.bandwidthMin[i]=bd.bandwidthMin[j];bd.bandwidthMin[j]=tmp;
    end
    i+=1
    j+=1
  end
end

function calcStatsDensity!(bd, root::Int)
  #println("calcStatsDensity! -- is running")
  calcStatsBall!(bd.bt, root)

  k = 0
  wtL=0.
  wtR = 0.
  wtT = 0.

  leftI = left(bd.bt, root);
  rightI=right(bd.bt, root);        # get children indices
  if (!validIndex(bd.bt,leftI) || !validIndex(bd.bt, rightI))
    return Union{}   # nothing to do if this
  end    #   isn't a parent node

  Ni  = bd.bt.dims*(root-1)
  NiL = bd.bt.dims*(leftI-1)
  NiR = bd.bt.dims*(rightI-1)
  wtL = bd.bt.weights[leftI]
  wtR = bd.bt.weights[rightI]
  wtT = wtL + wtR + eps(Float64);
  wtL /= wtT
  wtR /= wtT
  if (!bwUniform(bd))
    for k in 1:bd.bt.dims
      if (bd.bandwidthMax[NiL+k] > bd.bandwidthMax[NiR+k])
        bd.bandwidthMax[Ni+k] =  bd.bandwidthMax[NiL+k]
      else
        bd.bandwidthMax[Ni+k] =   bd.bandwidthMax[NiR+k]
      end
      if (bd.bandwidthMin[NiL+k] < bd.bandwidthMin[NiR+k])
        bd.bandwidthMin[Ni+k] = bd.bandwidthMin[NiL+k]
      else
        bd.bandwidthMin[Ni+k] = bd.bandwidthMin[NiR+k]
      end
    end
  end
  #switch(type) {
  #  case Gaussian:
    for k in 1:bd.bt.dims
      #@show k, Ni, NiL, NiR, round([wtL, wtR, wtT],2), round(bd.means[NiL+k],2), round(bd.means[NiR+k],2)
      bd.means[Ni+k]     = wtL * bd.means[NiL+k] + wtR * bd.means[NiR+k];
      bd.bandwidth[Ni+k] = wtL* (bd.bandwidth[NiL+k] + bd.means[NiL+k]*bd.means[NiL+k]) +
                        wtR* (bd.bandwidth[NiR+k] + bd.means[NiR+k]*bd.means[NiR+k]) -
                        bd.means[Ni+k]*bd.means[Ni+k];
      #@show round([bd.means[Ni+k], bd.bandwidth[Ni+k]],2)
    end
  #  break;
  #case Laplacian:
  #  for(unsigned int k=0; k < dims; k++) {
  #    means[Ni+k]     = wtL * means[NiL+k] + wtR * means[NiR+k];
  #    bandwidth[Ni+k] = wtL* (2*bandwidth[NiL+k]*bandwidth[NiL+k] + means[NiL+k]*means[NiL+k]) +
  #                      wtR* (2*bandwidth[NiR+k]*bandwidth[NiR+k] + means[NiR+k]*means[NiR+k]) -
  #                      means[Ni+k]*means[Ni+k];     // compute in terms of variance
  #    bandwidth[Ni+k] = sqrt(.5*bandwidth[Ni+k]);    //  then convert back to normal BW rep.
  #  }; break;
  #case Epanetchnikov:
  #  for(unsigned int k=0; k < dims; k++) {
  #    means[Ni+k]     = wtL * means[NiL+k] + wtR * means[NiR+k];
  #    bandwidth[Ni+k] = wtL* (.2*bandwidth[NiL+k]*bandwidth[NiL+k] + means[NiL+k]*means[NiL+k]) +
  #                      wtR* (.2*bandwidth[NiR+k]*bandwidth[NiR+k] + means[NiR+k]*means[NiR+k]) -
  #                     means[Ni+k]*means[Ni+k];     // compute in terms of variance
  #   bandwidth[Ni+k] = sqrt(5*bandwidth[Ni+k]);     //  then convert back to normal BW rep.
  #  }; break;
 #}
  return Union{}
end



# Create new matlab arrays and put them in the given structure
function makeBallTreeDensity(_pointsMatrix::Array{Float64,2}, _weightsMatrix::Array{Float64,1}, _bwMatrix::Array{Float64,1}, _type::DataType=GaussianKer)

  bt = makeBallTree(_pointsMatrix, _weightsMatrix, true)
  bt.calcStatsHandle = calcStatsDensity!
  bt.swapHandle = swapDensity!

  Np::Int = bt.num_points
  Nd::Int = bt.dims

  means = zeros(Nd*2*Np)
  means[(Np*Nd+1):end] = reshape(_pointsMatrix, Nd*Np, 1)

  if (size(_bwMatrix,2) == 1)
    multibw = 0
    bandwidth = zeros(Nd*2*Np)
    # @show Nd, Np, size(_bwMatrix), size(repeat(_bwMatrix,Np))
    bandwidth[(Np*Nd+1):end] = repeat(_bwMatrix,Np)
    bandwidthMin = bandwidth[(Np*Nd+1):end]
    bandwidthMax = bandwidth[(Np*Nd+1):end]
  else
    multibw = 1
    bandwidth = zeros(Nd*2*Np)
    bandwidthMin = zeros(Nd*2*Np)
    bandwidthMax = zeros(Nd*2*Np)
    bandwidth[(Np*Nd+1):end] = _bwMatrix
    bandwidthMin[(Np*Nd+1):end] = _bwMatrix
    bandwidthMax[(Np*Nd+1):end] = _bwMatrix
  end
  bd = BallTreeDensity(bt, _type, multibw, means, bandwidth, bandwidthMin, bandwidthMax, bt.calcStatsHandle, bt.swapHandle)
  bd.bt.data = bd # we need the circular reference for emulating polymorphism

  buildTree!(bd)

  return bd
end

function makeBallTreeDensity(_pointsMatrix::Array{Float64,2}, _weightsMatrix::Array{Float64,1}, _type::DataType=GaussianKer)
  makeBallTreeDensity(_pointsMatrix, _weightsMatrix, ones(size(_pointsMatrix,1)), _type)
end


# Figure out which of two children in this tree is closest to a given
# ball in another tree.  Returns the index in this tree of the closer
# child.
function closer!(bd::BallTreeDensity, myLeft::Int, myRight::Int, otherTree::BallTreeDensity, otherRoot::Int)
  if (myRight==NO_CHILD || otherRoot==NO_CHILD)
    return myLeft
  end
  dist_sq_l = 0.0
  dist_sq_r = 0.0
  #@fastmath @inbounds begin
    for i in 1:bd.bt.dims
      dist_sq_l += (center(otherTree.bt, otherRoot, i) - center(bd.bt, myLeft, i)) *
        (center(otherTree.bt, otherRoot, i) - center(bd.bt, myLeft, i));
      dist_sq_r += (center(otherTree.bt, otherRoot, i) - center(bd.bt, myRight, i)) *
        (center(otherTree.bt, otherRoot, i) - center(bd.bt, myRight, i));
    end
  #end

  if (dist_sq_l < dist_sq_r)
    return myLeft;
  else
    return myRight;
  end
end

function closer!(bd::BallTreeDensity, i::Int, j::Int, other_tree::BallTreeDensity)
  return closer!(bd,i,j,other_tree,root()) #othertree.root() = 1 anyway
end

# Perform a *slight* adjustment of the tree: move the points by delta, but
#   don't reform the whole tree; just fix up the statistics.
#
function movePoints!(bd::BallTreeDensity, delta::Array{Float64,1})
  for i in leafFirst(bd,root()):leafLast(bd,root())
    for k in 1:bt.dims                      # first adjust locations by delta
      bd.bt.centers[bd.bt.dims*(i-1)+k] += delta[ (getIndexOf(bd,i)-1)*bd.bt.dims + k ];
    end
  end
  for i in bd.bt.num_points:-1:1               # then recompute stats of
    calcStats!(bd.bt.data, i)                       #   parent nodes
  end
  calcStats!(bd.bt.data, root())                    #   and finally root node
  return Union{}
end

# Assumes newWeights is the right size (num_points)
function changeWeights!(bd::BallTreeDensity, newWeights::Array{Float64,1})
  for i in (bd.bt.num_points+1):bd.bt.num_points*2
    bd.bt.weights[i] = newWeights[ getIndexOf(bt, i) ];
  end

  for i in bd.bt.num_points:-1:1
    calcStats!(bd.bt.data,i)
  end
  calcStats!(bd.bt.data, root());
  return Union{}
end


function resample(p::BallTreeDensity, Np::Int=-1, ksType::String="lcv")
# resample(p,Np,KSType) -- construct a new estimate of the KDE p by sampling
#                      Np new points; determines a bandwidth by ksize(pNew,KSType)
#                      NOTE: KStype = 'discrete' resamples points by weight &
#                            preserves original kernel size
  if (Np==-1)
    Np = getNpts(p)
  end
  if (ksType == "discrete")
    q = kde!(getPoints(p),zeros(getDim(p),1),getWeights(p))
    samplePts,ind = sample(q,Np)
    if (size(p.bandwidth,2)>2*p.bt.num_points)
      ks = getBW(p,ind)
    else
      ks = getBW(p,1)
    end
    p2 = kde!(samplePts,ks);
  else
    samplePts, = sample(p,Np)
    p2 = kde!(samplePts, ksType)
  end
  return p2
end


function printBallTree(bd::BallTreeDensity, rnd=15)
  printBallTree(bd.bt);
  @show round(bd.means,rnd);
  @show round(bd.bandwidth,rnd);
  @show round(bd.bandwidthMin,rnd);
  @show round(bd.bandwidthMax,rnd);

  nothing
end

function test02()
  ## matlab results
  #centers=[30]=4.4, 3.6, 5, 3.2, 0.18, 2, 5.5, 6.5, 7.5, 2.6, 0.18, 2.5, 0, 0, 0, 1.9, 1.2, 3, 3.3, -0.81, 2, 4.6, -0.56, 1, 4, 5, 6, 7, 8, 9,
  #weights=[10]=1, 0.6, 0.4, 0.4, 0, 0.2, 0.2, 0.2, 0.2, 0.2,
  #ranges=[30]=2.6, 4.4, 4, 1.4, 0.98, 1, 1.5, 1.5, 1.5, 0.7, 0.98, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  #highest_leaf=[10]=9, 7, 9, 6, 0, 5, 6, 7, 8, 9,
  #lowest_leaf=[10]=5, 5, 8, 5, 0, 5, 6, 7, 8, 9,
  #permutation=[10]=0, 0, 0, 0, 0, 2, 1, 0, 3, 4,

  mus = [[4.6173,    3.2641,    1.8729, 4, 7 ]',
     [-0.5592,   -0.8088,    1.1610, 5, 8]',
     [1.,2.,3, 6, 9]']


  bd = makeBallTreeDensity(mus, 0.2*ones(5), 0.04*ones(3))
  printBallTree(bd)
end
