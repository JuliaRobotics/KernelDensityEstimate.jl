# capturing MSGibbs' final label selection per sample

##

using KernelDensityEstimate

const KDE = KernelDensityEstimate


##

X1 = kde!([1;2;3.0],[1.0;]);
X2 = kde!([0.5;1.5;2.5],[1.0;]);
X3 = kde!([4;5;6.0],[1.0;]);

##

glbs = KDE.makeEmptyGbGlb();
glbs.recordChoosen = true

p123 = *([X1;X2;X3], glbs=glbs, addEntropy=false );


## Analyze what happened


lc = glbs.labelsChoosen

for i in 1:3
  μ1 = getPoints(X1)[lc[i][1][2] - 3]
  μ2 = getPoints(X2)[lc[i][2][2] - 3]
  μ3 = getPoints(X3)[lc[i][3][2] - 3]

  Λμ = μ1 + μ2 + μ3

  @show μ = 1/3*Λμ
end


## check with the product we have...

p123 |> getPoints




#