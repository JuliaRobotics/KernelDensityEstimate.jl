

# Ball Tree

NO_CHILD = -1 # Maybe do int64(uint32(-1)) ??
FIELD_NAMES = ["D", "N", "centers", "ranges", "weights",
            "lower", "upper", "leftch", "rightch", "perm"]
nfields = 10;

type BallTree
  dims::Int64                     # dimension of data
  num_points::Int64               # of points
  centers::Array{Float64,1}       # ball centers, dims numbers per ball
  ranges::Array{Float64,1}        # bounding box ranges, dims per ball, dist from center to one side
  weights::Array{Float64,1}       # total weight in each ball

  left_child::Array{Int64,1}
  right_child::Array{Int64,1}     # left, right children; no parent indices
  lowest_leaf::Array{Int64,1}
  highest_leaf::Array{Int64,1}    # lower & upper leaf indices for each ball
  permutation::Array{Int64,1}     # point's position in the original data

  next::Int64                     # internal var for placing the non-leaf nodes

  swapHandle::Function
  calcStatsHandle::Function
  data
end


root() = 1
Ndim(bt::BallTree) = bt.dims
Npts(bt::BallTree) = bt.num_points
Npts(bt::BallTree, i::Int64) = bt.highest_leaf[i]-bt.lowest_leaf[i]+1


## todo make sure these two are working properly -- this is the problem, we are always assigning new memory for each return
#function center(bt::BallTree, i::Int64)
#    return bt.centers[((i-1)*bt.dims+1):end]
#end
center(bt::BallTree, i::Int64) = bt.centers[((i-1)*bt.dims+1):end]
center(bt::BallTree, i::Int64, k::Int) = bt.centers[((i-1)*bt.dims+k)]

#function rangeB(bt::BallTree, i::Int64)
#    return bt.ranges[((i-1)*bt.dims+1):end]
#end
rangeB(bt::BallTree, i::Int64) = bt.ranges[((i-1)*bt.dims+1):end]
rangeB(bt::BallTree, i::Int64, k::Int) = bt.ranges[((i-1)*bt.dims+k)]

#function weight(bt::BallTree, i::Int64)
#    return bt.weights[i]
#end
weight(bt::BallTree, i::Int64) = bt.weights[i]

#function isLeaf(bt::BallTree, ind::Int64)
#    return ind >= bt.num_points
#end
isLeaf(bt::BallTree, ind::Int64) = ind >= bt.num_points

validIndex(bt::BallTree, ind::Int64) = ((0<ind) && (ind <= 2*bt.num_points))

left(bt::BallTree, i::Int64) = bt.left_child[i]

right(bt::BallTree, i::Int64) = bt.right_child[i]

leafFirst(bt::BallTree, i::Int64) = bt.lowest_leaf[i]

leafLast(bt::BallTree, i::Int64) = bt.highest_leaf[i]

# Convert a BallTree::index to the numeric index in the original data
getIndexOf(bt::BallTree, i::Int64) = bt.permutation[i]

function swap!(data, _i::Int64, _j::Int64)
  return data.swapHandle(data, _i, _j)
end
function calcStats!(data, root::Int64)
  #@show "Fancy calcStats"
  return data.calcStatsHandle(data, root)
end


# Swap the ith leaf with the jth leaf.  Actually, only swap the
# weights, permutation, and centers, so only for swapping
# leaves. Will not swap ranges correctly and will not swap children
# correctly.
function swapBall!(bt::BallTree, _i::Int64, _j::Int64)
  #println("swapBall! -- (i,j)=$((_i,_j))")
  i = _i
  j = _j
  tmp = 0.
  if (i==j)
    return Union{}
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
    bt.centers[i]  = bt.centers[j];
    bt.centers[j] = tmp;
  end
end

# Find the dimension along which the leaves between low and high
# inclusive have the greatest variance
function  most_spread_coord(bt::BallTree, low::Int64, high::Int64)
  #BallTree::index dimension, point, max_dim;
  #double mean, variance, max_variance;
  #println("most_spread_coord -- low, high = $((low, high))")
  max_variance = 0
  max_dim = 1

  for dimension = 1:bt.dims
    mean = 0
    for point = (bt.dims*(low-1) + dimension):bt.dims:(bt.dims*(high-1))
      mean += bt.centers[point]
    end
    mean /= (high - low)

    variance = 0
    for point in (bt.dims*(low-1) + dimension):bt.dims:(bt.dims*(high-1))
      variance += (bt.centers[point] - mean) * (bt.centers[point] - mean);
    end
    if (variance > max_variance)
      max_variance = variance;
      max_dim = dimension;
    end
  end

  #println("most_spread_coord -- max_var=$(max_variance), return=$(max_dim)")
  return max_dim;
end

# straight from CLR, the unrandomized partition algorithm for
# quicksort.  Partitions the leaves from low to high inclusive around
# a random pivot in the given dimension.  Does not affect non-leaf
# nodes, but does relabel the leaves from low to high.
function partition!(bt::BallTree, dimension::Int64, low::Int64, high::Int64)
  pivot = low;  # not randomized, could set pivot to a random element

  while (low < high)
    while(bt.centers[bt.dims*(high-1) + dimension] >= bt.centers[bt.dims*(pivot-1) + dimension])
      high-=1
    end
    while(bt.centers[bt.dims*(low-1) + dimension] < bt.centers[bt.dims*(pivot-1) + dimension])
      low+=1
    end

    bt.swapHandle(bt, low, high)
    pivot = high
  end

  return high;
end


# Function to partition the data into two (equal-sized or near as possible)
#   sets, one of which is uniformly greater than the other in the given
#   dimension.
function select!(bt::BallTree, dimension::Int64, position::Int64, low::Int64, high::Int64)
    m = 0
    r = 0
    i = 0
  while low < high
    r = (floor(Int64,(low + high)/2))
    swap!(bt.data, r, low)
    m = low;
    for i in (low):high
      if (bt.centers[dimension+bt.dims*(i-1)] < bt.centers[dimension+bt.dims*(low-1)])
        m+=1
        swap!(bt.data, m, i);
      end
    end
    swap!(bt.data, low, m);
    if (m <= position) low=m+1 end
    if (m >= position) high=m-1 end
  end
  return Union{}
end


# Calculate the statistics of level "root" based on the statistics of
#   its left and right children.
function calcStatsBall!(bt::BallTree, root::Int64)
  #println("calcStatsBall! -- root=$(root)")
  Ni = 0
  NiL = 0
  NiR = 0
  d = 0

  leftI = left(bt, root)
  rightI=right(bt, root);   # get children indices
  if (!(validIndex(bt, leftI)) || !(validIndex(bt, rightI)))
    return Union{} # nothing to do if this
  end           #   isn't a parent node

  # figure out the center and ranges of this ball based on it's children
  maxi = 0.
  mini = 0.
  for d=1:bt.dims
    a = center(bt, leftI, d) + rangeB(bt, leftI, d)
    b = center(bt, rightI, d) + rangeB(bt, rightI, d)
    #@show (d, leftI, rightI, a, b)
    if (a > b)
      maxi = center(bt, leftI, d) + rangeB(bt, leftI, d)
    else
      maxi = center(bt, rightI, d) + rangeB(bt, rightI, d)
    end

    if (center(bt, leftI, d) - rangeB(bt, leftI, d) < center(bt, rightI, d) - rangeB(bt, rightI, d))
      mini = center(bt, leftI, d) - rangeB(bt, leftI, d)
    else
      mini = center(bt, rightI, d) - rangeB(bt, rightI, d)
    end

    #@show (d, mini,maxi);
    bt.centers[(root-1)*bt.dims+d] = (maxi+mini) / 2.0;
    bt.ranges[(root-1)*bt.dims+d] = (maxi-mini) / 2.0;
  end

  # if the left ball is the same as the right ball (should only
  # happen when calling the function directly with the same argument
  # twice), don't count the weight twice
  if (leftI != rightI)
    bt.weights[root] = bt.weights[leftI] + bt.weights[rightI];
  else
    bt.weights[root] = bt.weights[leftI]
  end
  return Union{}
end


# Given the leaves, build the rest of the tree from the top down.
# Split the leaves along the most spread coordinate, build two balls
# out of those, and then build a ball around those two children.
function buildBall!(bt::BallTree, low::Int64, high::Int64, root::Int64)
  #println("buildBall! -- (low, high, root)=$((low, high, root))")
  # special case for N=1 trees
  if (low == high)
    bt.lowest_leaf[root] = low;
    bt.highest_leaf[root] = high;
    bt.left_child[root] = low;

    # point right child to the same as left for calc stats, and then
    # point it to the correct NO_CHILD afterwards.  kinda kludgey
    bt.right_child[root] = high;
    calcStats!(bt.data, root)
    bt.right_child[root] = NO_CHILD;
    return Union{}
  end

  #BallTree::index coord, split, left, right;
  coord = most_spread_coord(bt, low, high); # find dimension of widest spread

  # split the current leaves into two groups, to build balls on them.
  # Choose the most spread coordinate to split them on, and make sure
  # there are the same number of points in each (+-1 for round off
  # error).
  split = (floor(Int64,(low + high) / 2))
  #@show coord, split, low, high
  select!(bt, coord, split, low, high)

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
    buildBall!(bt, low, split, left)
  end
  if(right != high)
    buildBall!(bt, split+1, high, right)
  end

  calcStats!(bt.data, root);
  return Union{}
end

# Public method to build the tree, just calls the private method with
# the proper starting arguments.
function buildTree!(bt::BallTree)
  #println("buildTree!(::BallTree) -- is running")
  i=bt.num_points
  for j in 1:bt.num_points
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

  buildBall!(bt, bt.num_points+1, 2*bt.num_points, 1); # chgd for indexing 1
  return Union{}
end

# New BallTree
function makeBallTree(_pointsMatrix::Array{Float64,2}, _weights::Array{Float64,1}, suppressBuildTree=false)
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
                  ones(Int64,2*Np), ones(Int64,2*Np), ones(Int64,2*Np),
                  ones(Int64,2*Np), zeros(Int64,2*Np),
                  0, swapBall!, calcStatsBall!, [])
  bt.data = bt

  if (Np > 0 && suppressBuildTree!=true)
    buildTree!(bt)
  end
  return bt
end

function printBallTree(bt::BallTree)
    @show round(bt.centers,15);
    @show round(bt.weights,15);
    @show round(bt.ranges,15);
    @show bt.left_child;
    @show bt.right_child;
    @show bt.highest_leaf;
    @show bt.lowest_leaf;
    @show bt.permutation;
    return Union{}
end

function test01()
  ## test point 01
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

  bt = makeBallTree(mus, 0.2*ones(5))
  printBallTree(bt);
end
