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
                0, 0)
end


mutable struct MSCompOpt
  pT::Float64
  tmpC::Float64
  tmpM::Float64
end

function levelInit!(glb::GbGlb)
  for j in 1:glb.Ndens
    glb.dNpts[j] = 1
    glb.levelList[j,1] = root(glb.trees[j])
  end
end

# this is a weird function -- not tested at all at this point
function initIndices!(glb::GbGlb)
  count=1
  for j in 1:glb.Ndens
    dNp = glb.dNpts[j]
    zz=glb.levelList[j,1]
    z = 1
    while z <= dNp # ++z, used to be a for loop
      glb.p[z] = weight(glb.trees[j], zz);  # init by sampling from weights
        z+=1
        if z<=dNp zz=glb.levelList[j,z]; end
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

function calcIndices!(glb::GbGlb)
  @fastmath @inbounds begin
    for j in 1:glb.Ndens
      for z in 1:glb.Ndim
        glb.particles[z+glb.Ndim*(j-1)] = mean(glb.trees[j],glb.ind[j], z)#[z];
        glb.variance[z+glb.Ndim*(j-1)]  = bw(glb.trees[j],glb.ind[j], z);
      end
    end
  end
end

function samplePoint!(X::Array{Float64,1}, glb::GbGlb, frm::Int)
  #counter = 1
  for j in 1:glb.Ndim
    mn=0.0; vn=0.0;
    for z in 1:glb.Ndens                             # Compute mean and variances of
      vn += 1.0/glb.variance[j+glb.Ndim*(z-1)]       #   product of selected particles
      mn += glb.particles[j+glb.Ndim*(z-1)]/glb.variance[j+glb.Ndim*(z-1)]
    end
    vn = 1.0/vn; mn *= vn;
    glb.rnptr += 1
    X[j+frm] = mn + sqrt(vn) * glb.randN[glb.rnptr] #counter       # then draw a sample from it
  end
  return Union{}
end

function levelDown!(glb::GbGlb)
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
  glb.levelList = glb.levelListNew              #   list and recycle the old
  glb.levelListNew=tmp
end

function sampleIndices!(X::Array{Float64,1}, cmoi::MSCompOpt, glb::GbGlb, frm::Int)#pT::Array{Float64,1}
  counter=1
  zz=0
  for j in 1:glb.Ndens
    dNp = glb.dNpts[j]    #trees[j].Npts();
    cmoi.pT = 0.0

    zz=glb.levelList[j,1]
    for z in 1:dNp
      glb.p[z] = 0.0
      for i in 1:glb.Ndim
        tmp = X[i+frm] - mean(glb.trees[j], zz, i)#[i]
        glb.p[z] += (tmp*tmp) / bw(glb.trees[j], zz, i)#[i]
        glb.p[z] += Base.log(bw(glb.trees[j], zz, i)) # Base.Math.JuliaLibm.log
      end
      glb.p[z] = exp( -0.5 * glb.p[z] ) * weight(glb.trees[j], zz)
      cmoi.pT += glb.p[z]
      zz = z<dNp ? glb.levelList[j,z+1] : zz
    end
    @simd for z in 1:dNp
        glb.p[z] /= cmoi.pT
    end
    @simd for z in 2:dNp
        glb.p[z] += glb.p[z-1]              # construct CDF and sample a
    end

    z=1
    zz=glb.levelList[j,z]
    while z<=(dNp-1)
      if (glb.randU[glb.ruptr] <= glb.p[z]) # counter
        break;                              # new kernel from jth density
      end
      z+=1
      if z<=dNp
        zz=glb.levelList[j,z]
      else
        error("This should never happend due to -1 MSGibbs01.jl")
      end
    end
    glb.ind[j] = zz
    counter+=1
    glb.ruptr += 1
  end
  calcIndices!(glb);                         # recompute particles, variance
end

## SLOWEST PIECE OF THE COMPUTATION -- TODO
# easy PARALLELs overhead here is much slower, already tried -- rather search for BLAS optimizations...
function makeFasterSampleIndex!(j::Int, cmo::MSCompOpt, glb::GbGlb)
  cmo.tmpC = 0.0
  cmo.tmpM = 0.0

  zz=glb.levelList[j,1]
  for z in 1:(glb.dNpts[j])
    glb.p[z] = 0.0
    for i in 1:glb.Ndim
      cmo.tmpC = bw(glb.trees[j], zz, i) + glb.Calmost[i]
      cmo.tmpM = mean(glb.trees[j], zz, i) - glb.Malmost[i]
      glb.p[z] += abs2(cmo.tmpM)/cmo.tmpC + log(cmo.tmpC) # This is the slowest piece
    end
    glb.p[z] = exp( -0.5 * glb.p[z] ) * weight(glb.trees[j].bt, zz) # slowest piece
    z < glb.dNpts[j] ? zz = glb.levelList[j,(z+1)] : nothing
  end

  nothing
end

function sampleIndex(j::Int, cmo::MSCompOpt, glb::GbGlb)
  cmo.pT = 0.0
  # determine product of selected particles from all but jth density
  for i in 1:glb.Ndim
    iCalmost = 0.0; iMalmost = 0.0;
    for k in 1:glb.Ndens
      if (k!=j) iCalmost += 1.0/glb.variance[i+glb.Ndim*(k-1)]; end
      if (k!=j) iMalmost += glb.particles[i+glb.Ndim*(k-1)]/glb.variance[i+glb.Ndim*(k-1)]; end
    end
    glb.Calmost[i] = 1/iCalmost;
    glb.Malmost[i] = iMalmost * glb.Calmost[i];
  end

  makeFasterSampleIndex!(j, cmo, glb)

  @simd for k in 1:glb.dNpts[j]
    cmo.pT += glb.p[k]
  end
  @simd for z in 1:glb.dNpts[j]
    glb.p[z] /= cmo.pT
  end
  @simd for z in 2:glb.dNpts[j]
    glb.p[z] += glb.p[z-1]
  end
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
  @simd for i in 1:glb.Ndim
    glb.particles[i+glb.Ndim*(j-1)] = mean(glb.trees[j], glb.ind[j], i)#[i]
    glb.variance[i+glb.Ndim*(j-1)]  = bw(glb.trees[j], glb.ind[j], i)#[i];
  end
end

function printGlbs(g::GbGlb, tag=Union{})
    if tag==Union{}
        println("=========================================================================")
    else
        println(string(tag,"================================================================"))
    end
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
                randU::Array{Float64,1}, randN::Array{Float64,1})

    glbs = makeEmptyGbGlb()
    glbs.Ndens = Ndens
    glbs.trees = trees
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
    glbs.Nlevels = floor(Int,((log(maxNp)/log(2))+1))
    glbs.particles = zeros(glbs.Ndim*Ndens)
    glbs.variance  = zeros(glbs.Ndim*Ndens)
    glbs.dNpts = zeros(Int,Ndens)
    glbs.levelList = ones(Int,Ndens,maxNp)
    glbs.levelListNew = ones(Int,Ndens,maxNp)
    cmo = MSCompOpt(0.0, 0.0, 0.0)
    cmoi = MSCompOpt(0.0, 0.0, 0.0)

    for s in 1:Np
        frm = ((s-1)*glbs.Ndim)

        levelInit!(glbs)
        initIndices!(glbs)
        calcIndices!(glbs)

        for l in 1:glbs.Nlevels
          samplePoint!(glbs.newPoints, glbs, frm)
          levelDown!(glbs);
          sampleIndices!(glbs.newPoints, cmoi, glbs, frm);

          for i in 1:Niter
            for j in 1:glbs.Ndens
              @fastmath @inbounds sampleIndex(j, cmo, glbs);
            end
          end
        end

        for j in 1:glbs.Ndens
          glbs.newIndices[(s-1)*glbs.Ndens+j] = getIndexOf(glbs.trees[j], glbs.ind[j])+1;  # return particle label
        end
        samplePoint!(glbs.newPoints, glbs, frm);
    end
    glbs = 0
    nothing
end

function prodAppxMSGibbsS(npd0::BallTreeDensity,
                          npds::Array{BallTreeDensity,1}, anFcns, anParams,
                          Niter::Int=5)
    # See  Ihler,Sudderth,Freeman,&Willsky, "Efficient multiscale sampling from products
    #         of Gaussian mixtures", in Proc. Neural Information Processing Systems 2003
    Ndens = length(npds)              # of densities
    Ndim  = npds[1].bt.dims           # of dimensions
    Np    = Npts(npd0)                # of points to sample

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
    Nlevels = floor(Int,(log(Float64(maxNp))/log(2.0))+1.0)  # how many levels to a balanced binary tree?

    # Generate enough random numbers to get us through the rest of this
    if true
      randU = rand(Int(Np*Ndens*(Niter+2)*Nlevels))
      randN = randn(Int(Ndim*Np*(Nlevels+1)))
    else
        randU = vec(readdlm("randU.csv"))
        randN = vec(readdlm("randN.csv"))
    end

    gibbs1(Ndens, npds, Np, Niter, points, indices, randU, randN);
    return reshape(points, Ndim, Np), reshape(indices, Ndens, Np)
end

function *(p1::BallTreeDensity, p2::BallTreeDensity)
  numpts = round(Int,(Npts(p1)+Npts(p2))/2)
  d = Ndim(p1)
  d != Ndim(p2) ? error("kdes must have same dimension") : nothing
  dummy = kde!(rand(d,numpts),[1.0]);
  pGM, = prodAppxMSGibbsS(dummy, [p1;p2], Union{}, Union{}, 5)
  return kde!(pGM)
end


function *(pp::Vector{BallTreeDensity})
  numpts = round(Int, Statistics.mean(Npts.(pp)))
  d = Ndim(pp[1])
  for p in pp
    d != Ndim(p) ? error("kdes must have same dimension") : nothing
  end
  dummy = kde!(rand(d,numpts),[1.0]);
  pGM, = prodAppxMSGibbsS(dummy, pp, Union{}, Union{}, 5)
  return kde!(pGM)
end
