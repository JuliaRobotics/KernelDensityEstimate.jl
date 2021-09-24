mutable struct GbGlb
    # REMEMBER -- all multi-dim arrays are column-indexing!  [i,j] => [j*N+i]
    particles::Array{Float64,2} # Array{Float64,1} # [Ndim x Ndens]      // means of selected particles
    variance::Array{Float64,2}  # Array{Float64,1} # [Ndim x Ndens]      //   variance of selected particles
    p::Array{Float64,1}         # [Np]                // probability of ith kernel
    ind::Array{Int,1}           # current indexes of MCMC step
    Malmost::Array{Float64,1}
    Calmost::Array{Float64,1}   #[Ndim x 1]   // Means & Cov. of all indices but j^th
    # Random number callback
    randU::Array{Float64,1}     # [Npoints * Ndens * (Niter+1)] uniformly distrib'd random variables
    randN::Array{Float64,1}     # [Ndim * Npoints] normally distrib'd random variables
    Ndim::Int
    Ndens::Int
    Nlevels::Int
    dNp::Int
    newPoints::Array{Float64,1}
    newWeights::Array{Float64,1}
    newIndices::Array{Int,2}    # final labels choosen per [sample, density] 
    trees::Array{BallTreeDensity,1}
    levelList::Array{Int,2}
    levelListNew::Array{Int,2}
    dNpts::Array{Int,1}         # ? how many points at current level in tree
    ruptr::Int
    rnptr::Int
    mn::Vector{Float64}         # ? only used in samplePoint as [1] -- something seems wrong with [1]
    vn::Vector{Float64}
    calclambdas::Vector{Float64}
    calcmu::Vector{Float64}
    # labelsChoosen[sample][densityId][level]
    labelsChoosen::Dict{Int,Dict{Int,Dict{Int,Int}}}
    recordChoosen::Bool
    partialDimMask::Vector{BitArray{1}}
end

function makeEmptyGbGlb(;recordChoosen::Bool=false)
  return GbGlb( zeros(0,0),
                zeros(0,0),
                zeros(0),
                zeros(Int,0),
                zeros(0),
                zeros(0),
                zeros(0),
                zeros(0),
                0,0,0,0,
                zeros(0),
                zeros(0),
                ones(Int,0,0),
                Vector{BallTreeDensity}(undef, 1),
                ones(Int,1,0),
                ones(Int,1,0),
                zeros(Int,0),
                0, 
                0,
                Float64[0.0;], Float64[0.0;],
                zeros(0), 
                zeros(0),
                Dict{Int,Dict{Int,Dict{Int,Int}}}(),
                recordChoosen,
                Vector{BitVector}()  )
end



function printGlbs(g::GbGlb, tag::String="")
  println(tag*"================================================================")
  println("Ndim=$(g.Ndim), Ndens=$(g.Ndens), Nlevels=$(g.Nlevels), dNp=$(g.dNp), dNpts=$(g.dNpts)")
  @show g.ind
  @show round.(g.particles, digits=2)
  @show round.(g.variance, digits=2)
  @show round.(g.p, digits=2)
  @show round.(g.Malmost, digits=2)
  @show round.(g.Calmost, digits=2)
  @show g.levelList
  @show g.levelListNew
  @show round.(g.newPoints, digits=4)
  @show g.newIndices
  @show g.partialDimMask
  nothing
end


mutable struct MSCompOpt
  pT::Float64
  tmpC::Float64
  tmpM::Float64
end


function updateGlbParticlesVariance!( glb::GbGlb, 
                                      j::Int, 
                                      idx::Int=0, 
                                      level::Int=0, 
                                      recordChoosen::Bool=false  )
  #
  for dim in 1:glb.Ndim
    # TODO should this skip partial dims???
    if !glb.partialDimMask[j][dim]
      # Skip inactive dim (partial) -- zero so that no information is added info = 1/var
      glb.particles[dim,j] = 0.0
      glb.variance[dim,j]  = 0.0
    else
      # regular case for active dimension
      glb.particles[dim,j] = mean(glb.trees[j], glb.ind[j], dim)
    # glb.particles[dim+glb.Ndim*(j-1)] = mean(glb.trees[j], glb.ind[j], dim)
      glb.variance[dim,j]  = bw(glb.trees[j], glb.ind[j], dim)
    end
  end

  if recordChoosen
    # Store the label selection
    glb.labelsChoosen[idx][j][level] = glb.trees[j].bt.permutation[glb.ind[j]]
  end

  nothing
end


"""
    $SIGNATURES

Recompute particles, variance.
"""
function calcIndices!(glb::GbGlb)::Nothing
  @fastmath @inbounds begin
    for j in 1:glb.Ndens
      updateGlbParticlesVariance!(glb, j)
    end
  end
  nothing
end

"""
    $SIGNATURES

Returns `Λ=Σ^{-1}` as sum of individual information matrices (inverse covariances -- i.e. bandwidths).

```math
Λ = Σ_i Λ_i
```
"""
getEuclidLambda(lambdas::Vector{Float64}) = sum(lambdas)

"""
    $SIGNATURES

Returns `Λμ`.

```math
Λμ = Σ_i (Λ_i*μ_i)
```
"""
function getEuclidMu( mus::Vector{Float64}, 
                      lambdas::Vector{Float64}, 
                      scale::Float64=1.0  )
  #
  lambdamu = 0.0
  for z in 1:length(mus)
    lambdamu += mus[z]*lambdas[z]
  end
  return scale*lambdamu
end


"""
    $SIGNATURES

Multiplication of Gaussians using leave in densities (if skip `skip` > 0).  
For on-manifold operations, set `getMu` and `getLambda` operations accordingly.

Notes
-----
- Used twice in samplePoint! (won't skip) and sampleIndex (will skip LOO).
- Assumes manifold `diffop` baked into `getMu`.
- use `j`th dimension
"""
function gaussianProductMeanCov!( glb::GbGlb,
                                  dim::Int,
                                  destMu::Vector{Float64},
                                  destCov::Vector{Float64},
                                  offset::Int,
                                  skip::Int,
                                  addop=+,  # currently not required -- baked into getMu
                                  getMu::Function=getEuclidMu,
                                  getLambda::Function=getEuclidLambda  )::Nothing
  #
  destMu[offset] = 0.0;
  destCov[offset] = 0.0;
  # suppress partials if no information available on this dimension
  checkpartials = (bv->bv[dim]).(glb.partialDimMask)
  # incorporate skip
  skip <= 0 ? nothing : (checkpartials[skip]=false)
  if !any(checkpartials)
    # zero introduces no artificial information
    return nothing
  end
  # Compute mean and variances (product) of selected particles

    # on manifold development
    @inbounds @fastmath @simd for j in 1:glb.Ndens
      if (j!=skip) && glb.partialDimMask[j][dim]
        # who populated glb.particles?
        glb.calclambdas[j] = 1.0/glb.variance[dim,j]
        glb.calcmu[j] = glb.particles[dim,j]
      else
        # adding zeros does not influence the result because `lambda_i = 0`
        glb.calclambdas[j] = 0.0
        glb.calcmu[j] = 0.0
      end
    end
    @inbounds destCov[offset] = getLambda(glb.calclambdas)
    @inbounds destCov[offset] = 1.0/destCov[offset]
    # μ = 1/Λ * Λμ  ## i.e. already scaled to mean only
    @inbounds destMu[offset] = getMu(glb.calcmu, glb.calclambdas, destCov[offset])

  nothing
end


"""
    $SIGNATURES

Calculate the likelihoods of left out kernel means (label_i's) against Normal distribution 
with mean `glb.Malmost` and variance `glb.Calmost` (to be verified)

Assumes
-------

Temporary product gaussian on-manifold has already been computed and stored in:
- `glb.Calmost` and `glb.Malmost`.

Notes
-----

- Incoming temporary product mean is `glb.Malmost` and covariance `glb.Calmost`.
- `j` is the left out density -- i.e. the density for which we want to select a new label (kernel)
- `cmo.tmpM` is the manifold difference from all dimensions prior to evaluating Gaussian pdf, and relates to levels of the k-d tree...
- `z` is index of leave one out density kernels -- i.e. evaluate likelihood of each kernel mean `glb.p[z] = p(μ_z)`
- New kernel will be selected randomly from likelihood weights in `glb.p[z:z]`.

Potential Concerns
------------------

- Make sure tmpC can happen on Euclidean vs user manifold (99% positive)

Development Notes
-----------------
- TODO, this function is one of the bottlenecks in computation
- Easy PARALLELs overhead here is much slower, already tried -- rather search for BLAS optimizations.  Could be related to memory allocation from multiple threads, worth retrying as Julia's automatic (and modern) threading model is enhanced.
"""
function makeFasterSampleIndex!(j::Int,
                                cmo::MSCompOpt,
                                glb::GbGlb,
                                muValue::Vector{Float64},
                                covValue::Vector{Float64},
                                offset::Int=0,
                                diffop=(-,),
                                doCalmost::Bool=true  )
  #
  # internal helper to for IR optimization to remove
  
  #
  cmo.tmpC = 0.0
  cmo.tmpM = 0.0
  cmo.pT = 0.0

  # zz=zz0
  zz=glb.levelList[j,1]

  # select active LOO dims from (not-j)
  dimmask = 0 .== ones(Int, glb.Ndim)
  for __j in 1:glb.Ndens
    __j == j ? continue : nothing
    dimmask .|= glb.partialDimMask[__j]
  end

  # iterate over kernels in the left out density
  for z in 1:(glb.dNpts[j])
    glb.p[z] = 0.0
    # compute mean `cmo.tmpM` and covariance `cmo.tmpC` across all KDE dimensions.
    for i in 1:glb.Ndim
      # skip both inactive dim on j, or skip inactive dim from all others in LOO (not-j)
      if !glb.partialDimMask[j][i] || !dimmask[i] ## tmpC tmpM only from partial LOO
        # Skip inactive dim (partial)
        continue
      end
      # tmpC is calculated on linear (Euclidean) manifold
      cmo.tmpC = bw(glb.trees[j], zz, i)
      doCalmost ? (cmo.tmpC += covValue[i]) : nothing
      # tmpM on-manifold differencing operation
      cmo.tmpM = diffop[i](mean(glb.trees[j], zz, i), muValue[i+offset])
      distr = (cmo.tmpM*cmo.tmpM)/cmo.tmpC
      if !isnan(distr)
        # This is the slowest piece
        glb.p[z] += distr
        glb.p[z] += log(cmo.tmpC)
      end
    end
    # final stage in Gaussian kernel evaluation
    @inbounds @fastmath glb.p[z] = exp( -0.5 * glb.p[z] ) * weight(glb.trees[j].bt, zz) # slowest piece
    
    # suppress NaNs
    if isnan(glb.p[z])
      glb.p[z] = 0.0
    end
    #
    cmo.pT += glb.p[z]
    z < glb.dNpts[j] ? (zz = glb.levelList[j,(z+1)]) : nothing
  end
  
  # stick with selection of others to preserve correlations
  if cmo.pT < 1e-99
    p_ = view(glb.p, 1:(glb.dNpts[j]))
    p_ .= weight(glb.trees[j].bt, zz)
    cmo.pT = sum(p_)
  end
  
  # Normalize the new probabilty for selecting a new kernel
  @simd for z in 1:glb.dNpts[j]
    glb.p[z] /= cmo.pT
  end
  
  # construct CDF for sampling a new kernel
  @simd for z in 2:glb.dNpts[j]
    glb.p[z] += glb.p[z-1]
  end

  nothing
end

function selectLabelOnLevel(glb::GbGlb, j::Int)

  dNp = glb.dNpts[j]    #trees[j].Npts();
  z = 1
  zz=glb.levelList[j,z]
  while z <= (dNp-1)
    # Selected a new kernel (`z`) from the `j`th density
    if (glb.randU[glb.ruptr] <= glb.p[z])  break;  end
    z+=1
    if z<=dNp
      zz=glb.levelList[j,z]
    end
  end

  # the new label in bt.means index
  glb.ind[j] = zz

  # increase randU counter
  glb.ruptr += 1 

  nothing
end


"""
    $SIGNATURES

??

Notes
-----
- This function does kernel evaluation internally.
- Does not have a loo skip step.
"""
function sampleIndices!(X::Array{Float64,1},
                        cmoi::MSCompOpt,
                        glb::GbGlb,
                        offset::Int,
                        diffop=(-,)  )
  #


  # calculate likelihood of all kernel means
  # zz=0
  for j in 1:glb.Ndens
    # how many kernels in this density
    makeFasterSampleIndex!(j, cmoi, glb, X, Vector{Float64}(undef, 0), offset, diffop, false)
    
    # sample a new kernel from CDF
    selectLabelOnLevel(glb, j)
  end

  # recompute particles, variance
  calcIndices!(glb);
  return nothing
end


"""
    $SIGNATURES

Sample new kernel in leave out density according to multiscale Gibbs sampling.

- calculate temporary product of leave in density components (labels previously selected)
- evaluate the likelihoods of the kernel means from the left out kernel on temporary product
- randomly select a new label (kernel_i) for left out density according to temporarily evaluated likelihoods

Notes
-----
- `j` is the left out density.
- Needs on-manifold getMu for product of two leave in density kernels.
- Needs diffop (on manifold difference for kernel evaluation).
- Sudderth PhD, p.139, Fig. 3.3, top-left operation
"""
function sampleIndex( j::Int,
                      cmo::MSCompOpt,
                      glb::GbGlb,
                      addop=(+,), 
                      diffop=(-,),
                      getMu=(getEuclidMu,),
                      getLambda=(getEuclidLambda,),
                      idx::Int=0,
                      level::Int=0  )
  #
  # determine product of selected kernel-labels from all but jth density (leave out)
  for i in 1:glb.Ndim
    gaussianProductMeanCov!(glb, i, glb.Malmost, glb.Calmost, i, j, addop[i], getMu[i], getLambda[i] )
    # indexMeanCovDens!(glb, i, j)
  end

  # evaluates the likelihoods of the left out density for each kernel mean, and stores in `glb.p[z]`.
  makeFasterSampleIndex!(j, cmo, glb, glb.Malmost, glb.Calmost, 0, diffop, true)

  selectLabelOnLevel(glb, j)

  # prep new particles and variances for calculation
  updateGlbParticlesVariance!(glb, j, idx, level, glb.recordChoosen)

  return nothing
end


"""
    $SIGNATURES

Sampling a point from the product of kernels (density components) listed in `glb.variance` and 
`glb.particles` without skipping a kernel (not a leave-one-out case).

Manifold defined by `addop`, `getMu`, and `getLambda`.
"""
function samplePoint!(X::Array{Float64,1},
                      glb::GbGlb,
                      idx::Int,
                      addop::Tuple=(+,),
                      getMu::Tuple=(getEuclidMu,),
                      getLambda::Tuple=(getEuclidLambda,),
                      addEntropy::Bool=true  )
  #
  for dim in 1:glb.Ndim
    # Calculate on-manifold mean and covariance.  Does not skip a density here -- i.e. skip = -1;  see `sampleIndex(...)`
    gaussianProductMeanCov!(glb, dim, glb.mn, glb.vn, 1, -1, addop[dim], getMu[dim], getLambda[dim] ) # getMeanCovDens!
    # then draw a sample from it
    glb.rnptr += 1
    
    # using .mn[1] is in conjunction with `gaussianProductMeanCov!(.., 1,..)` above
    X[dim+idx] = if addEntropy
      addop[dim](glb.mn[1], sqrt(glb.vn[1]) * glb.randN[glb.rnptr] )
    else
      glb.mn[1]
    end
    # @show addEntropy, dim, idx, X[dim+idx] 
  end
  return nothing
end

## Level and sampling operations

function levelInit!( glb::GbGlb, smplIdx::Int )
  for j in 1:glb.Ndens
    glb.dNpts[j] = 1
    glb.levelList[j,1] = root(glb.trees[j])
    if glb.recordChoosen
      glb.labelsChoosen[smplIdx][j] = Dict{Int,Int}()
    end
  end
end

function initIndices!( glb::GbGlb )
  # j is the number of incoming densities
  for j in 1:glb.Ndens
    dNp = glb.dNpts[j]
    zz=glb.levelList[j,1]
    z = 1
    while z <= dNp # ++z, used to be a for loop
      glb.p[z] = weight(glb.trees[j], zz);  # init by sampling from weights
        z+=1
        if z<=dNp
          zz=glb.levelList[j,z];
        end
    end

    for z in 2:dNp
        glb.p[z] += glb.p[z-1];
    end

    selectLabelOnLevel(glb, j)
  end
end


function levelDown!( glb::GbGlb )
  for j in 1:glb.Ndens
    z = 1
    for y in 1:(glb.dNpts[j])
      if validIndex(glb.trees[j], left(glb.trees[j], glb.levelList[j,y]))
        glb.levelListNew[j,z] = left(glb.trees[j], glb.levelList[j,y])
        z+=1
      end
      if validIndex(glb.trees[j], right(glb.trees[j], glb.levelList[j,y]))
        glb.levelListNew[j,z] = right(glb.trees[j], glb.levelList[j,y]);
        z+=1
      end
      if (glb.ind[j] == glb.levelList[j,y])                      # make sure ind points to
        glb.ind[j] = glb.levelListNew[j,z-1]                     #  a child of the old ind
      end
    end
    # ? lower down in the tree has more points
    glb.dNpts[j] = z-1
  end
  tmp = glb.levelList                            # make new list the current
  glb.levelList = glb.levelListNew               #   list and recycle the old
  glb.levelListNew=tmp
  nothing
end



function gibbs1(Ndens::Int, trees::Array{BallTreeDensity,1},
                Np::Int, Niter::Int,
                pts::Array{Float64,1}, ind::Array{Int},
                randU::Array{Float64,1}, randN::Array{Float64,1};
                addop=(+,), diffop=(-,),
                getMu=(getEuclidMu,),
                getLambda=(getEuclidLambda,),
                glbs = makeEmptyGbGlb(),
                addEntropy::Bool=true,
                ndims::Int=maximum(Ndim.(trees)),
                partialDimMask::AbstractVector{<:BitVector} = [ones(Int,ndims) .== 1 for i in 1:Ndens]  )
  #
  
  # number densities in the product
  glbs.Ndens = Ndens
  # the actual BTD data
  glbs.trees = trees
  # location for final posterior product samples
  glbs.newPoints = pts
  # REMEMBER -- all multi-dim arrays are column-indexing!  [i,j] => [j*N+i]
  glbs.newIndices = ind # reshape(ind, Ndens, :)
  # preemptive entropy for use during calculation
  glbs.randU = randU
  glbs.randN = randN
  # the dimension of the incoming densities
  glbs.Ndim = maximum(Ndim.(trees))
  glbs.partialDimMask = partialDimMask

  maxNp = 0                         # largest # of particles we deal with
  for tree in trees
      if (maxNp < Npts(tree))
          maxNp = Npts(tree)
      end
  end

    glbs.ind = ones(Int,Ndens)
    glbs.p = zeros(maxNp)
    glbs.Malmost = zeros(glbs.Ndim)
    glbs.Calmost = zeros(glbs.Ndim)
    glbs.calcmu = zeros(glbs.Ndens)
    glbs.calclambdas = zeros(glbs.Ndens)
    glbs.Nlevels = floor(Int,((log(maxNp)/log(2))+1))
    # working memory where kernels are multiplied together
    glbs.particles = zeros(glbs.Ndim, Ndens)
    glbs.variance  = zeros(glbs.Ndim, Ndens)
    glbs.dNpts = zeros(Int,Ndens)
    glbs.levelList = ones(Int,Ndens,maxNp)
    glbs.levelListNew = ones(Int,Ndens,maxNp)
    glbs.labelsChoosen = Dict{Int,Dict{Int,Dict{Int,Int}}}()

  cmo = MSCompOpt(0.0, 0.0, 0.0)
  cmoi = MSCompOpt(0.0, 0.0, 0.0)

    # loop for all output kernels in product (how many samples do you want from product)
    for s in 1:Np
        glbs.labelsChoosen[s] = Dict{Int, Dict{Int,Int}}()
        # index of where to put new sampled point in final posterior product
        frm = ((s-1)*glbs.Ndim)

        # initial assignments
        levelInit!(glbs, s)
        initIndices!(glbs)
        calcIndices!(glbs)

        # iterate down multi-scales of the Ball (k-d) tree (of the posterior belief?)
        for l in 1:glbs.Nlevels
          # multiply selected kernels from incoming densities, and sample a new point from the product.
          samplePoint!(glbs.newPoints, glbs, frm, addop, getMu, getLambda )

          # step a level down in the tree
          levelDown!(glbs);

          # ??
          sampleIndices!(glbs.newPoints, cmoi, glbs, frm, diffop);

          ## Sequential Gibbs logic where LOO element is rotated through all densities Niter times
          # After T iters, selected a kernel from each density
          @inbounds @fastmath for i in 1:Niter
            for j in 1:glbs.Ndens
              # pick a new label (kernel_i) from the LOO density (j) -- assumed leave in Gaussian product already computed
              sampleIndex(j, cmo, glbs, addop, diffop, getMu, getLambda, s, l );
            end
          end
        end

        for j in 1:glbs.Ndens
          # return particle label
          # REMEMBER -- all multi-dim arrays are column-indexing!  [i,j] => [j*N+i]
          glbs.newIndices[j,s] = getIndexOf(glbs.trees[j], glbs.ind[j])+1;
          # glbs.newIndices[(s-1)*glbs.Ndens+j] = getIndexOf(glbs.trees[j], glbs.ind[j])+1;
          
          # if glbs.recordChoosen
          #   glbs.labelsChoosen[s][j][glbs.Nlevels+1] = glbs.newIndices[j,s]
          # end
        end

        # take Gaussian product from a kernel component in all densities (inside samplePoint ??)
        # and then sample a value from new posterior kernel
        samplePoint!(glbs.newPoints, glbs, frm, addop, getMu, getLambda, addEntropy);
    end
    glbs = 0
    nothing
end


function prodAppxMSGibbsS(npd0::BallTreeDensity,
                          trees::Array{BallTreeDensity,1},
                          anFcns,
                          anParams,
                          Niter::Int )
  @warn "prodApproxMSGibbs has new keyword interface, use (..; Niter::Int=5 ) instead"
  prodAppxMSGibbsS( npd0,
                    trees,
                    anFcns,
                    anParams;
                    Niter=Niter )
end

function prodAppxMSGibbsS(npd0::BallTreeDensity,
                          trees::Array{BallTreeDensity,1},
                          anFcns,
                          anParams;
                          Niter::Int=3,
                          addop::T1=(+,),
                          diffop::T2=(-,),
                          getMu::T3=(getEuclidMu,),
                          getLambda::T4=(getEuclidLambda,),
                          glbs = makeEmptyGbGlb(),
                          addEntropy::Bool=true,
                          ndims::Integer=maximum(Ndim.(trees)),
                          Ndens = length(trees),              # of densities
                          Np    = Npts(npd0),                 # of points to sample
                          maxNp = maximum([Np; Npts.(trees)]),
                          Nlevels = floor(Int,(log(Float64(maxNp))/log(2.0))+1.0), # how many levels to a balanced binary tree?
                          randU = rand(Int(Np*Ndens*(Niter+2)*Nlevels)),
                          randN = randn(Int(ndims*Np*(Nlevels+1))),                
                          partialDimMask::AbstractVector{<:BitVector} = [ones(Int,ndims) .== 1 for i in 1:length(trees)],
                        ) where {T1<:Tuple,T2<:Tuple,T3<:Tuple,T4<:Tuple}
  #
  #
  
  # See  Ihler,Sudderth,Freeman,&Willsky, "Efficient multiscale sampling from products
  #         of Gaussian mixtures", in Proc. Neural Information Processing Systems 2003

  # prepare stack manifold add and diff operations functions (manifolds must match dimension)
  addopT = length(addop)!=ndims ? ([ (addop[1]) for i in 1:ndims]...,) : addop
  diffopT = length(diffop)!=ndims ? ([ (diffop[1]) for i in 1:ndims]...,) : diffop
  getMuT = length(getMu)!=ndims ? ([ getMu[1] for i in 1:ndims]...,) : getMu
  getLambdaT = length(getLambda)!=ndims ? ([ getLambda[1] for i in 1:ndims]...,) : getLambda

  # skipping analytic functions for now TODO ??
  UseAn = false
  #??pointsM = zeros(ndims, Np)
  points = zeros(ndims*Np)
  #??plhs[1] = mxCreateNumericMatrix(Ndens, Np, mxUINT32_CLASS, mxREAL);
  indices=ones(Int,Ndens, Np)
  # for tree in trees
  #   if (maxNp < Npts(tree))
  #     maxNp = Npts(tree)
  #   end
  # end


  # Generate enough random numbers to get us through the rest of this
  if false
    # FIXME using DelimitedFiles
    randU = vec(readdlm("randU.csv"))
    randN = vec(readdlm("randN.csv"))
  end

  gibbs1( Ndens, trees, Np, Niter, points, indices, 
          randU, randN, addop=addopT, diffop=diffopT, ndims=ndims,
          getMu=getMuT, getLambda=getLambdaT, glbs=glbs, 
          addEntropy=addEntropy, partialDimMask=partialDimMask );
  #
  return reshape(points, ndims, Np), indices
end



function *( trees::Vector{BallTreeDensity};
            glbs = makeEmptyGbGlb(),
            addEntropy::Bool=true )
  #

  # hack fix for #70
  if length(trees) == 1 && !addEntropy
    pts = getPoints(trees[1]) |> deepcopy
    return kde!(pts)
  end

  numpts = round(Int, Statistics.mean(Npts.(trees)))
  d = maximum(Ndim.(trees))
  for p in trees
    d != Ndim(p) ? error("kdes must have same dimension") : nothing
  end
  dummy = kde!(rand(d,numpts),[1.0]);
  pGM, = prodAppxMSGibbsS(dummy, trees, nothing, nothing, Niter=5, glbs=glbs, addEntropy=addEntropy)
  return kde!(pGM)
end

function *(p1::BallTreeDensity, p2::BallTreeDensity)
  return *([p1;p2])
  # numpts = round(Int,(Npts(p1)+Npts(p2))/2)
  # d = Ndim(p1)
  # d != Ndim(p2) ? error("kdes must have same dimension") : nothing
  # dummy = kde!(rand(d,numpts),[1.0]);
  # pGM, = prodAppxMSGibbsS(dummy, [p1;p2], nothing, nothing, Niter=5)
  # return kde!(pGM)
end
