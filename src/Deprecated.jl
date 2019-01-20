

function getMeanCovDens(glb::GbGlb, j::Int, skip::Int=-1)::Tuple{Float64, Float64}
  mn = 0.0;
  vn = 0.0;
  # Compute mean and variances (product) of selected particles
  @inbounds @fastmath @simd for z in 1:glb.Ndens
    if (z!=skip)
      # TODO: change to on-manifold operation
      vn += 1.0/glb.variance[j+glb.Ndim*(z-1)]
      mn += glb.particles[j+glb.Ndim*(z-1)]/glb.variance[j+glb.Ndim*(z-1)]
    end
  end
  vn = 1.0/vn;
  mn *= vn;
  return mn, vn
end


function indexMeanCovDens!(glb::GbGlb, i::Int, skip::Int=-1)
  getMeanCovDens!(glb, i, glb.Malmo,t, glb.Calmost, i, skip)
  # second
  # glb.Malmost[i], glb.Calmost[i] = getMeanCovDens(glb, i, skip)
  # first
  # iCalmost = 0.0;
  # iMalmost = 0.0;
  # @inbounds @fastmath @simd for z in 1:glb.Ndens
  #   # TODO change to on-manifold operation
  #   if (z!=skip)
  #     iCalmost += 1.0/glb.variance[i+glb.Ndim*(z-1)];
  #     iMalmost += glb.particles[i+glb.Ndim*(z-1)]/glb.variance[i+glb.Ndim*(z-1)];
  #   end
  # end
  # glb.Calmost[i] = 1.0/iCalmost;
  # glb.Malmost[i] = iMalmost * glb.Calmost[i];
  nothing
end
