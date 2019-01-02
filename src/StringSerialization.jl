function string(d::KernelDensityEstimate.BallTreeDensity)
  # TODO only supports single bandwidth per dimension at this point
  pts = getPoints(d)
  return "KDE:$(size(pts,2)):$(getBW(d)[:,1]):$(pts)"
end

function parsestringvector(str::AS; dlim=',') where {AS <: AbstractString}
  sstr = split(split(strip(str),'[')[end],']')[1]
  ssstr = strip.(split(sstr,dlim))
  parse.(Float64, ssstr)
end

function convert(::Type{BallTreeDensity}, str::AS) where {AS <: AbstractString}
  @assert occursin(r"KDE:", str) # ismatch
  sstr = strip.(split(str, ':'))
  N = parse(Int, sstr[2])
  bw = parsestringvector(sstr[3])
  dims = length(bw)
  ptrowstrs = split(sstr[4],';')
  @assert dims == length(ptrowstrs)
  pts = zeros(dims, N)
  for i in 1:dims
    pts[i,:] = parsestringvector(ptrowstrs[i], dlim=' ')
  end
  kde!(pts, bw)
end
# psubs = split(psubs, '[')[end]
# psubsub = split(psubs, ']')[1]
# pw = split(psubsub, ',')
