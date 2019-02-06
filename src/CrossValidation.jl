# conventional leave one out likelihood cross validation for bandwidth selection



function updateBandwidth!(bd::BallTreeDensity, bw::Array{Float64, 1})
    if (bd.multibandwidth==0)
        bd.bandwidth = bw
        bd.bandwidthMax = bd.bandwidthMin = bd.bandwidth[(bd.bt.num_points*bd.bt.dims+1):end]
    else
        error("updateBandwidth! -- multibandwidth==0 ELSE not implemented yet")
    end
end


function nLOO_LL(alpha::Float64, bd::BallTreeDensity, addop, diffop)
  # assume always a Gaussian
  alpha = alpha.^2

  updateBandwidth!(bd,bd.bandwidth*alpha)
  H = entropy(bd, addop, diffop)
  updateBandwidth!(bd,bd.bandwidth/alpha)

  return H
end

"""
    $(SIGNATURES)

GOLDEN   Minimize the nLOO_LL function for KDE bandwidth selection
of one variable using golden section search.

xmin, fmin = golden(npd, f, ax, bx, cx, tol) computes a local minimum
of f = KDE.nLOO_LL. xmin is the computed local minimizer of f and fmin is
f(xmin). xmin is computed to an relative accuracy of TOL.

The parameters ax, bx and cx must satisfy the following conditions:
ax < bx < cx, f(bx) < f(ax) and f(bx) < f(cx).

xmin satisfies ax < xmin < cx. golden is guaranteed to succeed if f
is continuous between ax and cx

Roman Geus, ETH Zuerich, 9.12.97
"""
function golden(bd::BallTreeDensity, ax::Float64, bx::Float64, cx::Float64, tol::Float64, addop, diffop)

    C = (3.0-sqrt(5.0))/2.0
    R = 1.0-C

    x0 = ax
    x3 = cx
    if (abs(cx-bx) > abs(bx-ax))
      x1 = bx
      x2 = bx + C*(cx-bx)
    else
      x1 = bx - C*(bx-ax)
      x2 = bx
    end

    #@show x1, x2
    f1 = nLOO_LL(x1,bd, addop, diffop)
    f2 = nLOO_LL(x2,bd, addop, diffop)

    k = 1;
    #tic()
    while abs(x3-x0) > tol*(abs(x1)+abs(x2))
    #  fprintf(1,'k=%4d, |a-b|=%e\n', k, abs(x3-x0));
      if f2 < f1
        x0 = x1
        x1 = x2
        x2 = R*x1 + C*x3   # x2 = x1+c*(x3-x1)
        f1 = f2
        f2 = nLOO_LL(x2,bd, addop, diffop)
      else
        x3 = x2
        x2 = x1
        x1 = R*x2 + C*x0   # x1 = x2+c*(x0-x2)
        f2 = f1
        f1 = nLOO_LL(x1,bd, addop, diffop)
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
    ssumt = sqrt.(sum( (2.0*rang[:,1:(bd.bt.num_points-1)]).^2 , dims=1))
    minm = minimum(ssumt)
    minm = max(minm, 1e-6)
    return minm, maxm
end

function ksize(bd::BallTreeDensity, addop, diffop )
    Nd = bd.bt.dims;
    Np = bd.bt.num_points;

    minm,maxm = neighborMinMax(bd)
    p = kde!(getPoints(bd),[(minm+maxm)/2.0],getWeights(bd))
    ks, dummy = golden(p,2.0*minm/(minm+maxm),1.0,2.0*maxm/(minm+maxm), 1e-2, addop, diffop)
    ks = ks * ( minm + maxm )/2.0
    npd = kde!(getPoints(p),[ks],getWeights(p), addop, diffop)
    return npd
end
