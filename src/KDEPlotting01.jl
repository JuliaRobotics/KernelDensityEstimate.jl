# Guide.manual_color_key("Legend", ["Points", "Line"], ["green", "deepskyblue"])

DOYTICKS=true

function draw1D!(bd::BallTreeDensity, bins::Union{Array{Float64,1},LinSpace{Float64}}, e, c::ASCIIString="blue", myStyle::ASCIIString="";
  xlbl="X")
  global DOYTICKS
  yV = evaluateDualTree(bd,bins)
  # clamp max y values
  yV[3.0 .< yV] = 3.0
  if e == Union{}
    if DOYTICKS
      e=Gadfly.plot(x=bins,y=yV,Geom.line, Gadfly.Theme(default_color=parse(Colorant,c)),Guide.xlabel(xlbl))
    else
      e=Gadfly.plot(x=bins,y=yV,Geom.line, Gadfly.Theme(default_color=parse(Colorant,c)),Guide.xlabel(xlbl),Guide.ylabel(""),Guide.yticks(ticks=nothing))
    end
  else
    push!(e.layers, layer(x=bins, y=yV, Geom.line, Gadfly.Theme(default_color=parse(Colorant,c)))[1])
  end

  #mx = maximum(y)
  #wts = getWeights(bd)
  return e
end

function plotKDE(bd::BallTreeDensity;
          N=200, c::Array{ASCIIString,1}=["black"],rmax=-Inf,rmin=Inf,
          xlbl="X")
    plotKDE([bd],N=N,c=c,rmax=rmax,rmin=rmin,xlbl=xlbl)
end

function getKDERange(bd::BallTreeDensity; extend::Float64=0.1)
  if (bd.bt.dims == 1)
    pts = getPoints(bd)
    rangeV = [minimum(pts),maximum(pts)]
  else
      return error("getKDERange(::BTD) -- multidimensional not implemented yet")
  end

  dr = extend*(rangeV[2]-rangeV[1])
  rangeV[1] = rangeV[1] - dr;
  rangeV[2] = rangeV[2] + dr;
  return rangeV
end

function plotKDE(darr::Array{BallTreeDensity,1};
      N::Int=200, c::Array{ASCIIString,1}=["black"],rmax=-Inf,rmin=Inf,
      xlbl="X")
    if (length(c)<2)
        c = repmat(c,length(darr))
    end
    H = Union{}
    i = 0
    for bd in darr
        i+=1
        dim = 1:bd.bt.dims
        if (bd.bt.dims == 1)
          rangeV = getKDERange(bd)
          if rangeV[1] > rmin  rangeV[1] = rmin end
          if rmax > rangeV[2]  rangeV[2] = rmax end
          H=draw1D!(bd,linspace(rangeV[1],rangeV[2],N), H, c[i],xlbl=xlbl) #,argsPlot,argsKDE
        else
            error("plotKDE(::BTD) -- multidimensional plotting not implemented yet")
        end
    end
    return H
end


# function drawProdElement!(bd::BallTreeDensity, bins::Union{Array{Float64,1},LinSpace{Float64}}, H, offs, height, mcmc; c::ASCIIString="blue", myStyle::ASCIIString="")
#   a = 0.4
#   if offs==0.0
#     c = "red"
#     a = 1.0
#   end
#
#   ZZ = evaluateDualTree(bd,bins) .+ height
#   YY = zeros(length(ZZ)) .+ (offs - 0.0*height) .+ mcmc
#
#   H=PyPlot.plot3D(YY,bins,ZZ, color=c, alpha=a)
#
#   if offs == 0.0
#     offs -= 0.5
#   end
#   offs -= 1.0
#   return H, offs
# end
#
# function plotKDEProd!(bdArr::Array{BallTreeDensity,1},minmax; h=0.0, mcmc=0.0, N=200)
#   offs = 0.0
#   H = []
#   for bd in bdArr
#         dim = 1:bd.bt.dims
#         if (bd.bt.dims == 1)
#           pts = getPoints(bd)
#           rangeV = [minimum(pts),maximum(pts)]
#           rangeV[1] = rangeV[1] - .15*(rangeV[2]-rangeV[1]);
#           rangeV[2] = rangeV[2] + .15*(rangeV[2]-rangeV[1]);
#           bins = linspace(rangeV[1],rangeV[2],N)
#           if (offs==0)
#             minval = rangeV[1]
#             #if minval > 0.0
#             #  minval = 0
#             #end
#             bins = linspace(minval,rangeV[2],N)
#           end
#           H,offs =drawProdElement!(bd, bins, H, offs, h, mcmc) #,argsPlot,argsKDE
#         else
#             error("plotKDEProd(::BTD) -- multidimensional plotting not implemented yet")
#         end
#
#         if rangeV[1] < minmax[1]
#           minmax[1] = rangeV[1]
#         end
#         if minmax[2]<rangeV[2]
#           minmax[2] = rangeV[2]
#         end
#     end
# end


macro ifund(exp)
    local e = :($exp)
    isdefined(e.args[1]) ? :($(e.args[1])) : :($(esc(exp)))
end

function drawHorDens(pDens::Array{BallTreeDensity,1}; N::Int=200,
                    gt=Union{}, lbls=ASCIIString[], extend::Float64=0.1)
    len = length(pDens)
    h = Array{Gadfly.Plot,1}(len)
    r = Array{RemoteRef,1}(len)

    if gt!=Union{}
      for i in 1:length(pDens)
        g = kde!(gt[i][1,:],[gt[i][2,1]])
        v = getKDERange(pDens[i],extend=extend)
        p = plotKDE([pDens[i];g],c=["black";"red"],rmin=v[1],rmax=v[2],xlbl=lbls[i])
        h[i] =p;
      end
    else
      for i in 1:length(pDens)
         h[i]=plotKDE([pDens[i]],xlbl=lbls[i])
      end
    end
    #println("draw $(len) plots")
    # [r[i] = @spawn plotKDE(pDens[i], N=N) for i in 1:len]
    # [h[i] = fetch(r[i]) for i in 1:len]

    hstack(h)
end

function stackMarginals(P::Array{BallTreeDensity,1}, m::Int64)
  ret = Array{BallTreeDensity,1}(length(P))
  for i in 1:length(P)
    ret[i] =  marginal(P[i],[m])#push!()
  end
  return ret
end


function vstackedPlots(plots::Array{Compose.Context,1})
    evalstr = ""
    for i in 1:length(plots)
        evalstr = string(evalstr, ",plots[$(i)]")
    end
    evalstr = string("vstack(",evalstr[2:end],")")
    eval(parse(evalstr))
end

# function vArrPotentials(potens::Array{BallTreeDensity,1})
#   vv = Array{Compose.Context,1}(length(potens))
#   i = 0
#   for p in potens
#       i+=1
#       vv[i] = investigatePoseKDE(p)
#   end
#   return vv
# end
