

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

function buildTree!(bd::BallTreeDensity, addop=(+,), diffop=(-,))
  #println("buildTree!(::BallTreeDensity)")
  buildTree!(bd.bt, addop, diffop)
  nothing
end

# Swap the ith leaf with the jth leaf.
function swapDensity!(bd::BallTreeDensity, i::Int, j::Int)
  if (i==j)
      return nothing
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
      tmp = bd.bandwidthMax[i];
      bd.bandwidthMax[i]=bd.bandwidthMax[j];
      bd.bandwidthMax[j]=tmp;
      tmp = bd.bandwidthMin[i];
      bd.bandwidthMin[i]=bd.bandwidthMin[j];
      bd.bandwidthMin[j]=tmp;
    end
    i+=1
    j+=1
  end
end

function calcStatsDensity!(bd, root::Int, addop=(+,), diffop=(-,))
  #println("calcStatsDensity! -- is running")
  calcStatsBall!(bd.bt, root, addop, diffop)

  k = 0
  wtL=0.
  wtR = 0.
  wtT = 0.

  leftI = left(bd.bt, root);
  rightI=right(bd.bt, root);        # get children indices
  if (!validIndex(bd.bt,leftI) || !validIndex(bd.bt, rightI))
    return nothing   # nothing to do if this
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
      bd.means[Ni+k]     = wtL * bd.means[NiL+k] + wtR * bd.means[NiR+k];
      bd.bandwidth[Ni+k] = wtL* (bd.bandwidth[NiL+k] + bd.means[NiL+k]*bd.means[NiL+k]) +
                           wtR* (bd.bandwidth[NiR+k] + bd.means[NiR+k]*bd.means[NiR+k]) -
                           bd.means[Ni+k]*bd.means[Ni+k];
    end
  return nothing
end



# Create new matlab arrays and put them in the given structure
function makeBallTreeDensity(_pointsMatrix::Array{Float64,2},
                             _weightsMatrix::Array{Float64,1},
                             _bwMatrix::Array{Float64,1},
                             _type::DataType=GaussianKer,
                             addop=(+,),
                             diffop=(-,)  )
  #
  bt = makeBallTree(_pointsMatrix, _weightsMatrix, true, addop, diffop)
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

  buildTree!(bd, addop, diffop)

  return bd
end

function makeBallTreeDensity(_pointsMatrix::Array{Float64,2},
                             _weightsMatrix::Array{Float64,1},
                             _type::DataType=GaussianKer,
                             addop=(+,),
                             diffop=(-,)  )
  makeBallTreeDensity(_pointsMatrix, _weightsMatrix, ones(size(_pointsMatrix,1)), _type, addop, diffop)
end


# Figure out which of two children in this tree is closest to a given
# ball in another tree.  Returns the index in this tree of the closer
# child.
function closer!(bd::BallTreeDensity,
                 myLeft::Int,
                 myRight::Int,
                 otherTree::BallTreeDensity,
                 otherRoot::Int,
                 addop=(+,),
                 diffop=(-,)  )
  #
  if (myRight==NO_CHILD || otherRoot==NO_CHILD)
    return myLeft
  end
  dist_sq_l = 0.0
  dist_sq_r = 0.0
  #@fastmath @inbounds begin
    for i in 1:bd.bt.dims
      dist_sq_l = addop[i](dist_sq_l, diffop[i](center(otherTree.bt, otherRoot, i), center(bd.bt, myLeft, i))^2  );
                                  # * diffop[i](center(otherTree.bt, otherRoot, i), center(bd.bt, myLeft, i))   );
      dist_sq_r = addop[i](dist_sq_r, diffop[i](center(otherTree.bt, otherRoot, i), center(bd.bt, myRight, i))^2  );
                                  # * diffop[i](center(otherTree.bt, otherRoot, i), center(bd.bt, myRight, i)));
    end
  #end

  if (dist_sq_l < dist_sq_r)
    return myLeft;
  else
    return myRight;
  end
end

function closer!(bd::BallTreeDensity, i::Int, j::Int, other_tree::BallTreeDensity, addop=(+,), diffop=(-,))
  return closer!(bd, i, j, other_tree, root(), addop, diffop) #othertree.root() = 1 anyway
end

# Perform a *slight* adjustment of the tree: move the points by delta, but
#   don't reform the whole tree; just fix up the statistics.
#
function movePoints!(bd::BallTreeDensity,
                     delta::Array{Float64,1},
                     addop=(+,),
                     diffop=(-,) )
  #
  for i in leafFirst(bd,root()):leafLast(bd,root())
    for k in 1:bt.dims                      # first adjust locations by delta
      bd.bt.centers[bd.bt.dims*(i-1)+k] = addop[k]( bd.bt.centers[bd.bt.dims*(i-1)+k], delta[ (getIndexOf(bd,i)-1)*bd.bt.dims + k ] );
    end
  end
  for i in bd.bt.num_points:-1:1               # then recompute stats of
    calcStats!(bd.bt.data, i, addop, diffop)                       #   parent nodes
  end
  calcStats!(bd.bt.data, root(), addop, diffop)                    #   and finally root node
  return nothing
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
  return nothing
end


function resample(p::BallTreeDensity, Np::Int=-1, ksType::Symbol=:lcv)
# resample(p,Np,KSType) -- construct a new estimate of the KDE p by sampling
#                      Np new points; determines a bandwidth by ksize(pNew,KSType)
#                      NOTE: KStype = 'discrete' resamples points by weight &
#                            preserves original kernel size
  if (Np==-1)
    Np = getNpts(p)
  end
  if (ksType == :discrete)
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
    p2 = kde!(samplePts)
  end
  return p2
end


function printBallTree(bd::BallTreeDensity, rnd=15)
  printBallTree(bd.bt);
  @show round.(bd.means,digits=rnd);
  @show round.(bd.bandwidth,digits=rnd);
  @show round.(bd.bandwidthMin,digits=rnd);
  @show round.(bd.bandwidthMax,digits=rnd);

  nothing
end
