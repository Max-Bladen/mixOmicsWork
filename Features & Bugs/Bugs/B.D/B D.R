library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
toxicity.spls <- pls(X, Y, ncomp = 3)

plotVar(toxicity.spls, cutoff = 0.99)
plotVar(toxicity.spls, cutoff = 0.0)
