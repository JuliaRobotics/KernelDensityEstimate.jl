using KernelDensityEstimate, Gadfly

toggleYTicks()

p = resample(kde!([0.0],[10.0]),300);
q = resample(kde!([-8.0; 13.0],[1.5]),300);

# easy product operator example, use prodAppxMSGibbsS(...) for more options and multiple terms
# dummy = kde!(rand(1,300),[1.0]);
# pGM, = prodAppxMSGibbsS(dummy, [p;q], nothing, nothing, 5)
pq = p*q


pl = plotKDE([p;q;pq],c=["blue";"red";"deepskyblue"],legend=["p(Θ | Y)";"p(Θ | Z)";"p(Θ | Y, Z)"],xlbl="p(Θ | Y, Z) ∝ p(Θ | Y) p(Θ | Z)")



r = resample(kde!([-35;-11.0; 26.0],[2.5]),300);

dummy = kde!(rand(1,300),[1.0]);
pGM, = prodAppxMSGibbsS(dummy, [p;q;r], nothing, nothing, Niter=5)
pqr = kde!(pGM);

pl2 = plotKDE([p;q;r;pqr],c=["blue";"red";"green";"deepskyblue"],legend=["p(Θ | Y)";"p(Θ | Z1)";"p(Θ | Z2)";"p(Θ | Y, Z:)"],xlbl="p(Θ | Y, Z:) ∝ p(Θ | Y) p(Θ | Z1) p(Θ | Z2)")


draw(PDF("consensus.pdf",12cm,12cm),vstack(pl,pl2))
