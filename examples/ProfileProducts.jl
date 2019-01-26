# profile products

using KernelDensityEstimate
using KernelDensityEstimatePlotting

function main()
  p = kde!(randn(100));
  q = kde!(randn(100));
  r = kde!(randn(100));

  numpts = Npts(p)
  dummy = kde!(rand(numpts),[1.0]);
  pGM = zeros(1,numpts)
  @time pGM[:,:], = prodAppxMSGibbsS(dummy, [p;q], nothing, nothing, Niter=5)

  pq = kde!(pGM)
  plotKDE([p;q;pq], c=["red";"green";"blue"])
end

main()

# using Profile
# Profile.clear_malloc_data()

main()
# exit()
# 0
