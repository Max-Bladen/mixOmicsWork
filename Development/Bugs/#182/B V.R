library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene

MyResult.pca <- pca(X) 
plotLoadings(MyResult.pca)
