mutable struct GbGlb
    # REMEMBER -- all multi-dim arrays are column-indexing!  [i,j] => [j*N+i]
    particles::Array{Float64,1} # [Ndim x Ndens]      // means of selected particles
    variance::Array{Float64,1}  # [Ndim x Ndens]      //   variance of selected particles
    p::Array{Float64,1}         # [Np]                // probability of ith kernel
    ind::Array{Int,1}                    # current indexes of MCMC step
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
    newIndices::Array{Int,1}
    trees::Array{BallTreeDensity,1}
    levelList::Array{Int,2}
    levelListNew::Array{Int,2}
    dNpts::Array{Int,1}
    ruptr::Int
    rnptr::Int
    mn::Vector{Float64}
    vn::Vector{Float64}
    calclambdas::Vector{Float64}
    calcmu::Vector{Float64}
end

function makeEmptyGbGlb()
   return GbGlb(zeros(0),
                zeros(0),
                zeros(0),
                zeros(Int,0),
                zeros(0),
                zeros(0),
                zeros(0),
                zeros(0),
                0,0,0,0,
                zeros(0),
                zeros(0),
                ones(Int,0),
                Vector{BallTreeDensity}(undef, 1),
                ones(Int,1,0),
                ones(Int,1,0),
                zeros(Int,0),
                0, 0,
                Float64[0.0;], Float64[0.0;],
                zeros(0), zeros(0)  )
end


mutable struct MSCompOpt
  pT::Float64
  tmpC::Float64
  tmpM::Float64
end


function updateGlbParticlesVariance!(glb::GbGlb, j::Int)::Nothing
  for dim in 1:glb.Ndim
    glb.particles[dim+glb.Ndim*(j-1)] = mean(glb.trees[j], glb.ind[j], dim)
    glb.variance[dim+glb.Ndim*(j-1)]  = bw(glb.trees[j], glb.ind[j], dim)
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
getEuclidLambda(lambdas::Vector{Float64})::Float64 = sum(lambdas)

"""
    $SIGNATURES

Returns `Λμ`.

```math
Λμ = Σ_i (Λ_i*μ_i)
```
"""
function getEuclidMu(mus::Vector{Float64}, lambdas::Vector{Float64}, scale::Float64=1.0)
  lambdamu = 0.0
  for z in 1:length(mus)
    lambdamu += mus[z]*lambdas[z]
  end
  return scale*lambdamu
end


"""
    $SIGNATURES

Multiplication of Gaussians using leave in densities (if skip `skip` > 0).  For on-manifold operations, set `getMu` and `getLambda` operations accordingly.

Notes
-----
- Used twice in samplePoint! (won't skip) and sampleIndex (will skip LOO).
- Assumes manifold `diffop` baked into `getMu`.
- use `j`th dimension
"""
function gaussianProductMeanCov!(glb::GbGlb,
                                 dim::Int,
                                 destMu::Vector{Float64},
                                 destCov::Vector{Float64},
                                 idx::Int,
                                 skip::Int,
                                 addop=+,  # currently not required -- baked into getMu
                                 getMu::Function=getEuclidMu,
                                 getLambda::Function=getEuclidLambda  )::Nothing
  destMu[idx] = 0.0;
  destCov[idx] = 0.0;
  # Compute mean and variances (product) of selected particles
  if false
    @inbounds @fastmath @simd for z in 1:glb.Ndens
      if (z!=skip)
        # TODO: change to on-manifold operation
        destCov[idx] += 1.0/glb.variance[dim+glb.Ndim*(z-1)]
        destMu[idx] = addop(destMu[idx], glb.particles[dim+glb.Ndim*(z-1)]/glb.variance[dim+glb.Ndim*(z-1)])
      end
    end
    destCov[idx] = 1.0/destCov[idx];
    destMu[idx] *= destCov[idx];
  else
    # on manifold development
    @inbounds @fastmath @simd for z in 1:glb.Ndens
      if (z!=skip)
        glb.calclambdas[z] = 1.0/glb.variance[dim+glb.Ndim*(z-1)]
        glb.calcmu[z] = glb.particles[dim+glb.Ndim*(z-1)]
      else
        # adding zeros does not influence the result because `lambda_i = 0`
        glb.calclambdas[z] = 0.0
        glb.calcmu[z] = 0.0
      end
    end
    @inbounds destCov[idx] = getLambda(glb.calclambdas)
    @inbounds destCov[idx] = 1.0/destCov[idx]
    # μ = 1/Λ * Λμ  ## i.e. already scaled to mean only
    @inbounds destMu[idx] = getMu(glb.calcmu, glb.calclambdas, destCov[idx])
    # destMu[idx] = destCov[idx]*getMu(glb.calcmu, glb.calclambdas)
    # destMu[idx] = getMu(glb.calcmu, glb.calclambdas)
  end
  nothing
end


"""
    $SIGNATURES

Calculate the likelihoods of left out kernel means (label_i's) against Normal distribution with mean `glb.Malmost` and variance `glb.Calmost` -- we think ...

Assumes
-------

Temporary product gaussian on manifold has already been computed and stored in `glb.Calmost` and `glb.Malmost`.

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
                                doCalmost::Bool=true  )::Nothing
  #
  cmo.tmpC = 0.0
  cmo.tmpM = 0.0
  cmo.pT = 0.0

  # zz=zz0
  zz=glb.levelList[j,1]

  # iterate over kernels in the left out density
  for z in 1:(glb.dNpts[j])
    glb.p[z] = 0.0
    # compute mean `cmo.tmpM` and covariance `cmo.tmpC` across all KDE dimensions.
    for i in 1:glb.Ndim
      # tmpC is calculated on linear (Euclidean) manifold
      cmo.tmpC = bw(glb.trees[j], zz, i)
      doCalmost ? (cmo.tmpC += covValue[i]) : nothing
      # tmpM on-manifold differencing operation
      cmo.tmpM = diffop[i](mean(glb.trees[j], zz, i), muValue[i+offset])
      # This is the slowest piece
      glb.p[z] += (cmo.tmpM*cmo.tmpM)/cmo.tmpC
      glb.p[z] += log(cmo.tmpC)
    end
    # final stage in Gaussian kernel evaluation
    @inbounds @fastmath glb.p[z] = exp( -0.5 * glb.p[z] ) * weight(glb.trees[j].bt, zz) # slowest piece
    #
    cmo.pT += glb.p[z]
    z < glb.dNpts[j] ? (zz = glb.levelList[j,(z+1)]) : nothing
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
                        diffop=(-,)  )::Nothing #pT::Array{Float64,1}
  #


  # calculate likelihood of all kernel means
  # zz=0
  for j in 1:glb.Ndens
    # how many kernels in this density
    dNp = glb.dNpts[j]    #trees[j].Npts();

    makeFasterSampleIndex!(j, cmoi, glb, X, Vector{Float64}(undef, 0), offset, diffop, false)

    # # ?? which level of the tree are you?
    # zz=glb.levelList[j,1]
    # for z in 1:dNp
    #   glb.p[z] = 0.0
    #
    #   for i in 1:glb.Ndim
    #     tmpC = bw(glb.trees[j], zz, i)
    #     tmpM = diffop[i]( mean(glb.trees[j], zz, i), X[i+offset] )
    #     glb.p[z] += (tmpM*tmpM) / tmpC
    #     glb.p[z] += log(bw(glb.trees[j], zz, i)) # Base.Math.JuliaLibm.log
    #   end
    #   # final step in calculating Gaussian kernel
    #   glb.p[z] = exp( -0.5 * glb.p[z] ) * weight(glb.trees[j].bt, zz)
    #   cmoi.pT += glb.p[z]
    #   z < dNp ? (zz = glb.levelList[j,(z+1)]) : nothing
    # end
    #
    # # Normalize the new probabilty for selecting a new kernel
    # @simd for z in 1:dNp
    #     glb.p[z] /= cmoi.pT
    # end
    #
    # # construct CDF for sampling a new kernel
    # @simd for z in 2:dNp
    #     glb.p[z] += glb.p[z-1]
    # end

    # sample a new kernel from CDF
    counter=1
    z=1
    zz=glb.levelList[j,z]
    while z<=(dNp-1)
      if (glb.randU[glb.ruptr] <= glb.p[z])
        # Selected a new kernel (`z`) from the `j`th density
        break;
      end
      z+=1
      if z<=dNp
        zz=glb.levelList[j,z]
      else
        error("This should never happen due to -1")
      end
    end
    glb.ind[j] = zz
    counter+=1
    glb.ruptr += 1 # increase randU counter
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

Sudderth PhD, p.139, Fig. 3.3, top-left operation
"""
function sampleIndex(j::Int,
                     cmo::MSCompOpt,
                     glb::GbGlb,
                     addop=(+,), diffop=(-,),
                     getMu=(getEuclidMu,),
                     getLambda=(getEuclidLambda,)  )::Nothing
  # determine product of selected kernel-labels from all but jth density (leave out)
  for i in 1:glb.Ndim
    gaussianProductMeanCov!(glb, i, glb.Malmost, glb.Calmost, i, j, addop[i], getMu[i], getLambda[i] )
    # indexMeanCovDens!(glb, i, j)
  end

  # evaluates the likelihoods of the left out density for each kernel mean, and stores in `glb.p[z]`.
  makeFasterSampleIndex!(j, cmo, glb, glb.Malmost, glb.Calmost, 0, diffop, true)

  zz=glb.levelList[j,1]
  z=1
  while z<=(glb.dNpts[j]-1)
    if (glb.randU[glb.ruptr] <= glb.p[z]) break;  end   #1  #   a new kernel from the jth
    z+=1
    if z<=glb.dNpts[j]
      zz=glb.levelList[j,z]
    end
  end
  glb.ind[j] = zz;
  glb.ruptr += 1

  # prep new particles and variances for calculation
  updateGlbParticlesVariance!(glb, j)
  # @simd for dim in 1:glb.Ndim
  #   glb.particles[dim+glb.Ndim*(j-1)] = mean(glb.trees[j], glb.ind[j], dim)
  #   glb.variance[dim+glb.Ndim*(j-1)]  = bw(glb.trees[j], glb.ind[j], dim)
  # end
  return nothing
end

## Level and sampling operations

function levelInit!(glb::GbGlb)
  for j in 1:glb.Ndens
    glb.dNpts[j] = 1
    glb.levelList[j,1] = root(glb.trees[j])
  end
end

function initIndices!(glb::GbGlb)
  count=1
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
    zz=glb.levelList[j,1]
    z = 1
    while z <= (dNp-1)
      if (glb.randU[glb.ruptr] <= glb.p[z]) break; end #count
        z+=1
        if z<dNp zz=glb.levelList[j,z]; end
    end
    glb.ind[j] = zz;
    count+=1
    glb.ruptr += 1
  end
end

"""
    $SIGNATURES

Sampling a point from the product of kernels (density components) listed in `glb.variance` and `glb.particles` without skipping a kernel (not a leave-one-out case).

Manifold defined by `addop`, `getMu`, and `getLambda`.
"""
function samplePoint!(X::Array{Float64,1},
                      glb::GbGlb,
                      frm::Int,
                      addop::T1=(+,),
                      getMu::T2=(getEuclidMu,),
                      getLambda::T3=(getEuclidLambda,) )::Nothing  where {T1<:Tuple, T2<:Tuple, T3<:Tuple}
  #counter = 1
  for j in 1:glb.Ndim
    # Calculate on-manifold mean and covariance.  Does not skip a density here -- i.e. skip = -1;  see `sampleIndex(...)`
    gaussianProductMeanCov!(glb, j, glb.mn, glb.vn, 1, -1, addop[j], getMu[j], getLambda[j] ) # getMeanCovDens!
    # then draw a sample from it
    glb.rnptr += 1
    X[j+frm] = addop[j](glb.mn[1], sqrt(glb.vn[1]) * glb.randN[glb.rnptr] ) #counter
  end
  return nothing
end

function levelDown!(glb::GbGlb)::Nothing
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
    glb.dNpts[j] = z-1
  end
  tmp = glb.levelList                            # make new list the current
  glb.levelList = glb.levelListNew               #   list and recycle the old
  glb.levelListNew=tmp
  nothing
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
end

function gibbs1(Ndens::Int, trees::Array{BallTreeDensity,1},
                Np::Int, Niter::Int,
                pts::Array{Float64,1}, ind::Array{Int,1},
                randU::Array{Float64,1}, randN::Array{Float64,1};
                addop=(+,), diffop=(-,),
                getMu=(getEuclidMu,),
                getLambda=(getEuclidLambda,) )::Nothing
    #

    glbs = makeEmptyGbGlb()
    glbs.Ndens = Ndens
    glbs.trees = trees
    # location for final posterior product samples
    glbs.newPoints = pts
    glbs.newIndices = ind
    glbs.randU = randU
    glbs.randN = randN
    glbs.Ndim = trees[1].bt.dims

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
    glbs.particles = zeros(glbs.Ndim*Ndens)
    glbs.variance  = zeros(glbs.Ndim*Ndens)
    glbs.dNpts = zeros(Int,Ndens)
    glbs.levelList = ones(Int,Ndens,maxNp)
    glbs.levelListNew = ones(Int,Ndens,maxNp)
    cmo = MSCompOpt(0.0, 0.0, 0.0)
    cmoi = MSCompOpt(0.0, 0.0, 0.0)

    # loop for all output kernels in product (how many samples do you want from product)
    for s in 1:Np
        # index of where to put new sampled point in final posterior product
        frm = ((s-1)*glbs.Ndim)

        # initial assignments
        levelInit!(glbs)
        initIndices!(glbs)
        calcIndices!(glbs)

        # iterate down multi-scales of the Ball (k-d) tree (of the posterior belief?)
        for l in 1:glbs.Nlevels
          # multiply selected kernels from incoming densities, and sample a new point from the product.
          samplePoint!(glbs.newPoints, glbs, frm, addop, getMu, getLambda )

          # step a level down in the tree`
          levelDown!(glbs);

          # ??
          sampleIndices!(glbs.newPoints, cmoi, glbs, frm, diffop);

          # After T iters, selected a kernel from each density
          @inbounds @fastmath for i in 1:Niter
            for j in 1:glbs.Ndens
              # pick a new label (kernel_i) from the LOO density (j) -- assumed leave in Gaussian product already computed
              sampleIndex(j, cmo, glbs, addop, diffop, getMu, getLambda );
            end
          end
        end

        for j in 1:glbs.Ndens
          glbs.newIndices[(s-1)*glbs.Ndens+j] = getIndexOf(glbs.trees[j], glbs.ind[j])+1;  # return particle label
        end

        # take Gaussian product from a kernel component in all densities (inside samplePoint ??)
        # and then sample a value from new posterior kernel
        samplePoint!(glbs.newPoints, glbs, frm, addop, getMu, getLambda );
    end
    glbs = 0
    nothing
end


function prodAppxMSGibbsS(npd0::BallTreeDensity,
                          npds::Array{BallTreeDensity,1},
                          anFcns,
                          anParams,
                          Niter::Int )
  @warn "prodApproxMSGibbs has new keyword interface, use (..; Niter::Int=5 ) instead"
  prodAppxMSGibbsS(npd0,
                   npds,
                   anFcns,
                   anParams;
                   Niter=Niter )
end

function prodAppxMSGibbsS(npd0::BallTreeDensity,
                          npds::Array{BallTreeDensity,1},
                          anFcns,
                          anParams;
                          Niter::Int=3,
                          addop::T1=(+,),
                          diffop::T2=(-,),
                          getMu::T3=(getEuclidMu,),
                          getLambda::T4=(getEuclidLambda,)  ) where {T1<:Tuple,T2<:Tuple,T3<:Tuple,T4<:Tuple}
    # See  Ihler,Sudderth,Freeman,&Willsky, "Efficient multiscale sampling from products
    #         of Gaussian mixtures", in Proc. Neural Information Processing Systems 2003
    Ndens = length(npds)              # of densities
    Ndim  = npds[1].bt.dims           # of dimensions
    Np    = Npts(npd0)                # of points to sample

    # prepare stack manifold add and diff operations functions (manifolds must match dimension)
    addopT = length(addop)!=Ndim ? ([ (addop[1]) for i in 1:Ndim]...,) : addop
    diffopT = length(diffop)!=Ndim ? ([ (diffop[1]) for i in 1:Ndim]...,) : diffop
    getMuT = length(getMu)!=Ndim ? ([ getMu[1] for i in 1:Ndim]...,) : getMu
    getLambdaT = length(getLambda)!=Ndim ? ([ getLambda[1] for i in 1:Ndim]...,) : getLambda

    # skipping analytic functions for now TODO ??
    UseAn = false
    #??pointsM = zeros(Ndim, Np)
    points = zeros(Ndim*Np)
    #??plhs[1] = mxCreateNumericMatrix(Ndens, Np, mxUINT32_CLASS, mxREAL);
    indices=ones(Int,Ndens*Np)
    maxNp = Np
    for tree in npds
        if (maxNp < Npts(tree))
            maxNp = Npts(tree)
        end
    end

    # how many levels to a balanced binary tree?
    Nlevels = floor(Int,(log(Float64(maxNp))/log(2.0))+1.0)

    # Generate enough random numbers to get us through the rest of this
    if true
      randU = rand(Int(Np*Ndens*(Niter+2)*Nlevels))
      randN = randn(Int(Ndim*Np*(Nlevels+1)))
    else
      randU = vec(readdlm("randU.csv"))
      randN = vec(readdlm("randN.csv"))
    end

    gibbs1(Ndens, npds, Np, Niter, points, indices, randU, randN, addop=addopT, diffop=diffopT, getMu=getMuT, getLambda=getLambdaT );
    return reshape(points, Ndim, Np), reshape(indices, Ndens, Np)
end



function *(pp::Vector{BallTreeDensity})
  numpts = round(Int, Statistics.mean(Npts.(pp)))
  d = Ndim(pp[1])
  for p in pp
    d != Ndim(p) ? error("kdes must have same dimension") : nothing
  end
  dummy = kde!(rand(d,numpts),[1.0]);
  pGM, = prodAppxMSGibbsS(dummy, pp, nothing, nothing, Niter=5)
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
