

# Ball Tree

global NO_CHILD = -1
# FIELD_NAMES = ["D", "N", "centers", "ranges", "weights",
#             "lower", "upper", "leftch", "rightch", "perm"]
# nfields = 10;

mutable struct BallTree
  dims::Int                     # dimension of data
  num_points::Int               # of points
  centers::Array{Float64,1}     # ball centers, dims numbers per ball
  ranges::Array{Float64,1}      # bounding box ranges, dims per ball, dist from center to one side
  weights::Array{Float64,1}     # total weight in each ball

  left_child::Array{Int,1}
  right_child::Array{Int,1}     # left, right children; no parent indices
  lowest_leaf::Array{Int,1}
  highest_leaf::Array{Int,1}    # lower & upper leaf indices for each ball
  permutation::Array{Int,1}     # point's position in the original data

  next::Int                     # internal var for placing the non-leaf nodes

  swapHandle::Function
  calcStatsHandle::Function
  data
end


root() = 1
Ndim(bt::BallTree) = bt.dims
Npts(bt::BallTree) = bt.num_points
Npts(bt::BallTree, i::Int) = bt.highest_leaf[i]-bt.lowest_leaf[i]+1


center(bt::BallTree, i::Int) = bt.centers[((i-1)*bt.dims+1):end]
center(bt::BallTree, i::Int, k::Int) = bt.centers[((i-1)*bt.dims+k)]


rangeB(bt::BallTree, i::Int) = bt.ranges[((i-1)*bt.dims+1):end]
rangeB(bt::BallTree, i::Int, k::Int) = bt.ranges[((i-1)*bt.dims+k)]


weight(bt::BallTree, i::Int) = bt.weights[i]


isLeaf(bt::BallTree, ind::Int) = ind >= bt.num_points

validIndex(bt::BallTree, ind::Int) = ((0<ind) && (ind <= 2*bt.num_points))

left(bt::BallTree, i::Int) = bt.left_child[i]

right(bt::BallTree, i::Int) = bt.right_child[i]

leafFirst(bt::BallTree, i::Int) = bt.lowest_leaf[i]

leafLast(bt::BallTree, i::Int) = bt.highest_leaf[i]

# Convert a BallTree::index to the numeric index in the original data
getIndexOf(bt::BallTree, i::Int) = bt.permutation[i]

function swap!(data, _i::Int, _j::Int)
  return data.swapHandle(data, _i, _j)
end
function calcStats!(data, root::Int, addop=(+,), diffop=(-,))
  #@show "Fancy calcStats"
  return data.calcStatsHandle(data, root , addop, diffop )
end


# Swap the ith leaf with the jth leaf.  Actually, only swap the
# weights, permutation, and centers, so only for swapping
# leaves. Will not swap ranges correctly and will not swap children
# correctly.
function swapBall!(bt::BallTree, _i::Int, _j::Int)
  #println("swapBall! -- (i,j)=$((_i,_j))")
  i = _i
  j = _j
  tmp = 0.
  if (i==j)
    return nothing
  end

  # swap weights
  tmp = bt.weights[i]
  bt.weights[i] = bt.weights[j]
  bt.weights[j] = tmp

  # swap perm
  k = bt.permutation[i]
  bt.permutation[i] = bt.permutation[j]
  bt.permutation[j] = k;

  # swap centers
  i = (i-1)*bt.dims
  j = (j-1)*bt.dims;
  for k in 1:bt.dims
    i+=1
    j+=1
    tmp = bt.centers[i];
    bt.centers[i] = bt.centers[j];
    bt.centers[j] = tmp;
  end
end

# Find the dimension along which the leaves between low and high
# inclusive have the greatest variance
function  most_spread_coord(bt::BallTree, low::Int, high::Int, addop=(+,), diffop=(-,))
  #BallTree::index dimension, point, max_dim;
  #double mean, variance, max_variance;
  #println("most_spread_coord -- low, high = $((low, high))")
  max_variance = 0
  max_dim = 1
  w = 1.0/(high-low)

  for dimension in 1:bt.dims
    # compute mean
    mean = 0
    for point = (bt.dims*(low-1) + dimension):bt.dims:(bt.dims*(high-1))
      # scale each value by 1/N
      mean = addop[dimension](mean, w*bt.centers[point])
    end

    # now that we have the mean, compute variance
    variance = 0
    for point in (bt.dims*(low-1) + dimension):bt.dims:(bt.dims*(high-1))
      variance += (diffop[dimension](bt.centers[point], mean))^2 # * (bt.centers[point] - mean);
    end

    # update variance if needed
    if (variance > max_variance)
      max_variance = variance;
      max_dim = dimension;
    end
  end

  #println("most_spread_coord -- max_var=$(max_variance), return=$(max_dim)")
  return max_dim;
end

# """
#     $SIGNATURES
#
# unrandomized partition algorithm for quicksort.
# Partitions the leaves from low to high inclusive around
# a random pivot in the given dimension.  Does not affect non-leaf
# nodes, but does relabel the leaves from low to high.
#
# straight from CLR (? Intro to algorithms - Leierson, Rivest),
# """
# function partition!(bt::BallTree, dimension::Int, low::Int, high::Int, diffop)
#   pivot = low;  # not randomized, could set pivot to a random element
#
#   while (low < high)
#     while(bt.centers[bt.dims*(high-1) + dimension] >= bt.centers[bt.dims*(pivot-1) + dimension])
#       high-=1
#     end
#     while(bt.centers[bt.dims*(low-1) + dimension] < bt.centers[bt.dims*(pivot-1) + dimension])
#       low+=1
#     end
#
#     bt.swapHandle(bt, low, high)
#     pivot = high
#   end
#
#   return high;
# end
# function partition!(bt::BallTree, dimension::Int, low::Int, high::Int, diffop)
#   pivot = low;  # not randomized, could set pivot to a random element
#
#   while (low < high)
#     while diffop(bt.centers[bt.dims*(high-1) + dimension], bt.centers[bt.dims*(pivot-1) + dimension]) >= 0.0
#       high-=1
#     end
#     while diffop(bt.centers[bt.dims*(low-1) + dimension], bt.centers[bt.dims*(pivot-1) + dimension]) < 0.0
#       low+=1
#     end
#
#     bt.swapHandle(bt, low, high)
#     pivot = high
#   end
#
#   return high;
# end

# Function to partition the data into two (equal-sized or near as possible)
#   sets, one of which is uniformly "greater" than the other in the given
#   dimension.
function select!(bt::BallTree, dimension::Int, position::Int, low::Int, high::Int, diffop::T1) where {T1 <: Tuple}
    m = 0
    r = 0
    i = 0
  while low < high
    r = (floor(Int,(low + high)/2))
    swap!(bt.data, r, low)
    m = low;
    for i in (low):high
      if diffop[dimension](bt.centers[dimension+bt.dims*(i-1)], bt.centers[dimension+bt.dims*(low-1)]) < 0.0
        m+=1
        swap!(bt.data, m, i);
      end
    end
    swap!(bt.data, low, m);
    if (m <= position) low=m+1 end
    if (m >= position) high=m-1 end
  end
  return nothing
end

"""
    $SIGNATURES

Return "smaller" and "larger" of two child nodes.
"""
function getMiniMaxi(bt::BallTree,
                     leftI::Int,
                     rightI::Int,
                     d::Int,
                     addop,
                     diffop  )::Tuple{Float64, Float64}
  #

  mini = 0.0
  maxi = 0.0

  # whos the most positive
  a = addop[d]( center(bt, leftI, d), rangeB(bt, leftI, d) )
  b = addop[d]( center(bt, rightI, d), rangeB(bt, rightI, d) )
  if (a > b)
    maxi = addop[d]( center(bt, leftI, d), rangeB(bt, leftI, d) )
  else
    maxi = addop[d]( center(bt, rightI, d), rangeB(bt, rightI, d) )
  end

  c = diffop[d](center(bt, leftI, d), rangeB(bt, leftI, d) )
  c2 = diffop[d](center(bt, rightI, d), rangeB(bt, rightI, d))
  if c < c2
    mini = diffop[d]( center(bt, leftI, d), rangeB(bt, leftI, d) )
  else
    mini = diffop[d]( center(bt, rightI, d), rangeB(bt, rightI, d) )
  end

  return mini, maxi
end

# Calculate the statistics of level "root" based on the statistics of
#   its left and right children.
function calcStatsBall!(bt::BallTree, root::Int, addop=(+,), diffop=(-,))
  #println("calcStatsBall! -- root=$(root)")
  Ni = 0
  NiL = 0
  NiR = 0
  d = 0

  # get children indices
  leftI = left(bt, root)
  rightI=right(bt, root)

  # @show round.(bt.centers, digits=3)

  # nothing to do if this isn't a parent node
  if (!(validIndex(bt, leftI)) || !(validIndex(bt, rightI)))
    return nothing
  end

  # figure out the center and ranges of this ball based on it's children
  maxi = 0.
  mini = 0.
  for d=1:bt.dims

    # get which child is mini or maxi
    mini, maxi = getMiniMaxi(bt, leftI, rightI, d, addop, diffop)

    #@show (d, mini,maxi);
    # @show (root-1)*bt.dims+d
    # @show maxi, mini

    # TODO implicit Euclidean comparison (not right!)
    # computing the parent node halfspan
    # bt.ranges[(root-1)*bt.dims+d] = diffop[d](maxi, mini) / 2.0;  # Basic Euclidean
    halfspan = diffop[d](maxi, mini) / 2.0; # Better on-manifold
    bt.ranges[(root-1)*bt.dims+d] = halfspan

    # Computing the parent node center
    # @show "naive Euclidean mean", addop[d](maxi, mini) / 2.0
    thecenter = addop[d](mini, halfspan)
    bt.centers[(root-1)*bt.dims+d] = thecenter # addop[d](maxi, mini) / 2.0;
  end

  # if the left ball is the same as the right ball (should only
  # happen when calling the function directly with the same argument
  # twice), don't count the weight twice
  if (leftI != rightI)
    bt.weights[root] = bt.weights[leftI] + bt.weights[rightI];
  else
    bt.weights[root] = bt.weights[leftI]
  end
  # error("finishing with root=$root")


  return nothing
end


# Given the leaves, build the rest of the tree from the top down.
# Split the leaves along the most spread coordinate, build two balls
# out of those, and then build a ball around those two children.
function buildBall!(bt::BallTree,
                    low::Int,
                    high::Int,
                    root::Int,
                    addop=(+,),
                    diffop=(-,)  )::Nothing
  global NO_CHILD
  #println("buildBall! -- (low, high, root)=$((low, high, root))")
  # special case for N=1 trees
  if (low == high)
    bt.lowest_leaf[root] = low;
    bt.highest_leaf[root] = high;
    bt.left_child[root] = low;

    # point right child to the same as left for calc stats, and then
    # point it to the correct NO_CHILD afterwards.  kinda kludgey
    bt.right_child[root] = high;
    calcStats!(bt.data, root, addop, diffop)
    bt.right_child[root] = NO_CHILD;
    return nothing
  end

  #BallTree::index coord, split, left, right;
  coord = most_spread_coord(bt, low, high, addop, diffop); # find dimension of widest spread

  # split the current leaves into two groups, to build balls on them.
  # Choose the most spread coordinate to split them on, and make sure
  # there are the same number of points in each (+-1 for round off
  # error).
  split = (floor(Int,(low + high) / 2))
  #@show coord, split, low, high
  select!(bt, coord, split, low, high, diffop)

  # an alternative is to use partition, but that doesn't deal well
  # with repeated numbers and it doesn't split into balanced sets.
  # split = partition(coord, low, high);

  # if the left sub-tree is just one leaf, don't make a new non-leaf
  # node for it, just point left_idx directly to the leaf itself.
  if (split <= low)
    left = low
  else
    left = bt.next;
    bt.next+=1
  end

  # same for the right
  if (split+1 >= high)
    right = high
  else
    right = bt.next
    bt.next+=1
  end

  bt.lowest_leaf[root]  = low;
  bt.highest_leaf[root] = high;
  bt.left_child[root]   = left;
  bt.right_child[root]  = right;

  # build sub-trees if necessary
  if(left != low)
    buildBall!(bt, low, split, left, addop, diffop)
  end
  if(right != high)
    buildBall!(bt, split+1, high, right, addop, diffop)
  end

  calcStats!(bt.data, root, addop, diffop);
  return nothing
end

# Public method to build the tree, just calls the private method with
# the proper starting arguments.
function buildTree!(bt::BallTree, addop=(+,), diffop=(-,))
  global NO_CHILD
  #println("buildTree!(::BallTree) -- is running")
  i=bt.num_points
  @inbounds for j in 1:bt.num_points
    for k in 1:bt.dims
      bt.ranges[i*bt.dims+k] = 0
    end
    i+=1
    bt.lowest_leaf[i] = i
    bt.highest_leaf[i] = i
    bt.left_child[i] = i
    bt.right_child[i] = NO_CHILD #int64(uint32(NO_CHILD))
    bt.permutation[i] = j
  end
  bt.next = 2

  buildBall!(bt, bt.num_points+1, 2*bt.num_points, 1, addop, diffop)
  return nothing
end

# New BallTree
function makeBallTree(_pointsMatrix::Array{Float64,2},
                      _weights::Array{Float64,1},
                      suppressBuildTree=false,
                      addop=(+,),
                      diffop=(-,)  )
  # get fields from input arguments
  Nd = size(_pointsMatrix,1);
  Np = size(_pointsMatrix,2);

  _points = reshape(_pointsMatrix, Nd*Np, 1)
  centers = zeros(Nd*2*Np)
  ranges = zeros(Nd*2*Np)
  weights = zeros(2*Np)
  centers[(Np*Nd+1):end] = _points[1:end]
  weights[(Np+1):end] = _weights[1:end]

  bt = BallTree(Nd, Np, centers, ranges, weights,
                  ones(Int,2*Np), ones(Int,2*Np), ones(Int,2*Np),
                  ones(Int,2*Np), zeros(Int,2*Np),
                  0, swapBall!, calcStatsBall!, [])
  bt.data = bt

  if (Np > 0 && suppressBuildTree!=true)
    buildTree!(bt, addop, diffop)
  end
  return bt
end

function printBallTree(bt::BallTree)
    @show round.(bt.centers,dims=15);
    @show round.(bt.weights,dims=15);
    @show round.(bt.ranges,dims=15);
    @show bt.left_child;
    @show bt.right_child;
    @show bt.highest_leaf;
    @show bt.lowest_leaf;
    @show bt.permutation;
    return nothing
end
