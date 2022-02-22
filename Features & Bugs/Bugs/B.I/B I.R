library(mixOmics)

data(liver.toxicity)

X <- liver.toxicity$gene
X <- rbind(X,X)[1:72, 1:48]

X[sample(nrow(X), 10), sample(ncol(X), 10)] <- NA

grid.keepX<-c(seq(5,30,5))
tune.spca.result<-tune.spca(X, ncomp=3, folds=4, test.keepX=grid.keepX, nrepeat=2)
