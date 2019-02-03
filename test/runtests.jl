using KernelDensityEstimate
using Test
using LinearAlgebra, Statistics, DelimitedFiles

include("testKnownConstructions.jl")

# parse the output from matlab process
function parseMatPrintKDE(filename::String)
  fid = open(filename,"r")
  dict = Dict{String, Array{Float64,1}}()
  for line in eachline(fid)
    twos = split(split(line,"\n")[1],"=")
    s = (split(split(twos[2],"[")[2],"]")[1])
    dict[twos[1]] = vec(readdlm(IOBuffer(s),','))
  end
  close(fid)
  return dict
end

function constructBTD(dict::Dict{String, Array{Float64,1}})
  refbt = BallTree(dict["dims"][1],dict["num_points"][1],
  dict["centers"], #
  dict["ranges"], # =
  dict["weights"], #
  # can use . syntax when Julia 0.5 support is dropped
  map(t -> round(Int(t)),dict["left_child"]), # =
  map(t -> round(Int(t)),dict["right_child"]), # =
  map(t -> round(Int(t)),dict["lowest_leaf"]), # =
  map(t -> round(Int(t)),dict["highest_leaf"]), # =
  map(t -> round(Int(t)),dict["permutation"]), # =
  0, +, +, nothing )

  refbtd = BallTreeDensity(refbt, nothing, 0,
  dict["means"], # =
  dict["bandwidth"], # =
  dict["bwMin"],
  dict["bwMax"],
  +, +)
  return refbtd
end

function testSubtract(b1::BallTreeDensity, b2::BallTreeDensity, tol=1e-10)
    #@show b1.bt.dims, b2.bt.dims, b1.bt.num_points, b2.bt.num_points
    if (b1.bt.dims - b2.bt.dims) != 0 || (b1.bt.num_points - b2.bt.num_points) != 0
        error("testSubtract -- b1 and b2 are not the same size.")
    end
    function testDiffs(a1::Array{Float64},a2::Array{Float64}; tol=1e-10, name="")
        n = norm(a1 - a2)
        if n > tol
          @show name, n
          error("testSubtract -- b1 and b2 are not the same.")
        end
        nothing
    end
    function testInds(a1::Array{Int},a2::Array{Int})
        for i in 1:length(a1)
            res = 0
            if (a1[i] < 0)
                res = (a2[i]+1)
            else
                res = (a1[i]+1-a2[i])
            end
            if res != 0
              error("testSubtract -- b1 and b2 have different indices.")
            end
        end
        nothing
    end
    testDiffs(b1.bt.centers, b2.bt.centers, name="centers", tol=tol)
    testDiffs(b1.bt.ranges, b2.bt.ranges, name="ranges", tol=tol)
    testDiffs(b1.bt.weights, b2.bt.weights, name="weights", tol=tol)
    testInds(b1.bt.left_child, b2.bt.left_child)
    testInds(b1.bt.right_child, b2.bt.right_child)
    testInds(b1.bt.lowest_leaf, b2.bt.lowest_leaf)
    testInds(b1.bt.highest_leaf, b2.bt.highest_leaf)
    testInds(b1.bt.permutation[(b1.bt.num_points+1):end], b2.bt.permutation[(b2.bt.num_points+1):end])
    testDiffs(b1.means, b2.means, name="means", tol=tol)
    testDiffs(b1.bandwidth, b2.bandwidth, tol=tol)
    testDiffs(b1.bandwidthMin, b2.bandwidthMin, tol=tol)
    testDiffs(b1.bandwidthMax, b2.bandwidthMax, tol=tol)
    println("Success")
    true
end

function testdatapath(ss::String)
  # joinpath(Pkg.dir("KernelDensityEstimate"),"test","testdata",ss)
  joinpath(dirname(pathof(KernelDensityEstimate)),"..","test","testdata",ss)
end

function UnitTest1D01()

print("Running UnitTest1D01...")
p = kde!([.1,.45,.55,3.8],[0.08])

d = parseMatPrintKDE(testdatapath("test1DResult.txt"))
refbtd = constructBTD(d)

# printBallTree(p)
testSubtract(refbtd, p, 1e-5)

end


function UnitTest1Dlcv01()

  print("Running UnitTest1Dlcv01...")
  x = vec(readdlm(testdatapath("test1Dlcv100.txt"))')
  p = kde!(collect(x) )

  d = parseMatPrintKDE(testdatapath("test1Dlcv100Result.txt"))
  refbtd = constructBTD(d)

  # printBallTree(p)
  testSubtract(refbtd, p, 1e-4)

end

function UnitTest2D01()
  print("Running UnitTest2D01...")
  pts = [[0.5172, 0.7169, 0.4049]';
         [0.0312, 1.0094, 2.0204]']
  p = kde!(pts,[0.1])

  d = parseMatPrintKDE(testdatapath("test2DResult.txt"))
  refbtd = constructBTD(d)

  # printBallTree(p)
  testSubtract(refbtd, p, 1e-5)
end

function UnitTest2Dlcv01()
  print("Running UnitTest2Dlcv01...")
  x = readdlm(testdatapath("test2Dlcv100.txt"))'
  p = kde!(collect(x) )

  d = parseMatPrintKDE(testdatapath("test2Dlcv100Result.txt"))
  refbtd = constructBTD(d)

  # printBallTree(p)
  testSubtract(refbtd, p, 1e-4)
end

function UnitTest2Dvar01()
  print("Running UnitTest2Dvar01...")
  pts = [[0.5172, 7.169, 4.049]';
         [0.0312, 10.0094, -2.0204]']
  p=kde!(pts,[0.1; 1.0]);
  d = parseMatPrintKDE(testdatapath("test2DvarResult.txt"))
  refbtd = constructBTD(d)

  # printBallTree(p)
  testSubtract(refbtd, p, 1e-4)
end

function UnitTest2Dvarlcv01()
  print("Running UnitTest2Dvarlcv01...")
  x = readdlm(testdatapath("test2Dvarlcv100.txt"))'
  p = kde!(x )

  d = parseMatPrintKDE(testdatapath("test2Dvarlcv100Result.txt"))
  refbtd = constructBTD(d)

  # printBallTree(p)
  testSubtract(refbtd, p, 2e-3)
end

function testProds(;D=3,M=6,N=100,n=100, dev=1.0, MCMC=5)
  P = BallTreeDensity[];
  [push!(P, kde!(dev*randn(D,N))) for i in 1:M];
  dummy = kde!(randn(D,n),[1.0]);
  pGM, = prodAppxMSGibbsS(dummy, P, nothing, nothing, Niter=MCMC);
  sum(abs, pGM) < 1e-14 ? error("testProds -- prodAppxMSGibbsS, nothing in pGM, len $(length(P))") : nothing
  prodDev = sqrt(dev^(2*M)/(M*(dev^2)))
  T1 = norm(Statistics.mean(pGM,dims=2)) < 1.0*prodDev
  T2 = true
  for i in 1:D
    tv = Statistics.std(pGM[i,:])
    lv = (0.66*prodDev < tv < 1.33*prodDev)
    T2 = T2 && lv
  end
  T1 && T2
end

function rangeTestProds(;D=3,M=6,N=100,n=100, dev=1.0, MCMC=5)
  @show v = [testProds(D=D,M=M,N=N,n=n, dev=dev, MCMC=MCMC) for i in 1:10]
  sum(map(Int,v)) >= 5
end

function rangeUnitTests()
  passt = true
  @show passt = passt && rangeTestProds(D=2,M=2);
  @show passt = passt && rangeTestProds(D=2,M=4);
  @show passt = passt && rangeTestProds(D=2,M=6);
  @show passt = passt && rangeTestProds(D=3,M=6, MCMC=10);
  @show passt = passt && rangeTestProds(D=4,M=6, n=200, MCMC=10);
  # @show passt = passt && rangeTestProds(D=3,M=10, MCMC=10);
  @show passt = passt && rangeTestProds(D=3,M=5,N=300);
  @show passt = passt && rangeTestProds(D=2,M=7,n=300);
  @show passt = passt && rangeTestProds(D=3,M=2,MCMC=100);
  return passt
end

function intgAppxGaussianOffs(;offs::Float64=0.0, N::Int=201, dim::Int=1)
  p = kde!(randn(dim,100));
  pts = randn(dim,150)
  pts[1,:] .+= offs
  q = kde!(pts);
  return intersIntgAppxIS(p,q, N=N)
end

function integralAppxUnitTests()
  testflag = true
  @show a = intgAppxGaussianOffs(offs=0.0,dim=1)
  @show testflag = testflag && 0.2 < a < 0.35
  @show a = intgAppxGaussianOffs(offs=1.0,dim=1,N=1000)
  @show testflag = testflag && 0.1 < a < 0.3
  @show a = intgAppxGaussianOffs(offs=-2.0,dim=1,N=1000)
  @show testflag = testflag && 0.01 < a < 0.17
  @show a = intgAppxGaussianOffs(offs=0.0,dim=2)
  @show testflag = testflag && 0.05 < a < 0.15

  return testflag
end

function testRand()
  println("testing rand functionality")
  p = kde!(rand(2,100))
  pts = rand(p,100)
  return true
end


@test UnitTest1D01()
@test UnitTest1Dlcv01()
@test UnitTest2D01()
# UnitTest2Dlcv01()
@test UnitTest2Dvar01()
#UnitTest2Dvarlcv01()
@test rangeUnitTests()
@test integralAppxUnitTests()

@test testRand()



@testset "test string and back conversion" begin

p = kde!(randn(2,3))
ps = string(p)
pp = convert(BallTreeDensity, ps)

@test norm(getPoints(pp)-getPoints(p)) < 1e-4
@test norm(getBW(pp)-getBW(p)) < 1e-4

end
