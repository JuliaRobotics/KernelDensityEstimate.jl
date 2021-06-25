

using KernelDensityEstimate


##

@testset "test partial product" begin

##

pts1 = rand(2,100) .+ 10.0
mask1 = 1 .== ones(2)
mask1[2] = false

pts2 = rand(2,100)
mask2 = 1 .== ones(2)

pts3 = rand(2,100) .- 10.0
mask3 = 1 .== ones(2)
mask3[1] = false


##

P1 = kde!(pts1)
P2 = kde!(pts2)
P3 = kde!(pts3)

bw1 = getBW(P1)[:,1]
bw2 = getBW(P2)[:,1]
bw3 = getBW(P3)[:,1]

pts1[2,:] .= 9999999.0
pts3[1,:] .= 9999999.0

P1 = kde!(pts1, bw1)
P3 = kde!(pts3, bw3)


mask = [mask1, mask2, mask3]

##

dummy = kde!(rand(2,100));

pGM, = prodAppxMSGibbsS(dummy, [P1;P2;P3], nothing, nothing, partialDimMask=mask )

# P = kde!(pGM)

@test 80 < sum(0 .< pGM[1,:] .< 10 )

@test 80 < sum(-10 .< pGM[2,:] .< 0)


##

end


#