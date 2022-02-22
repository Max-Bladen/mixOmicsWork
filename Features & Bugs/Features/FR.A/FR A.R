# FR.A - project new data onto existing PCA() object components

library(mixOmics)

data("liver.toxicity")
X <- liver.toxicity$gene
test <- sample(1:nrow(X), 16)
X.test <- X[test, ]
X.train <- X[-test, ]

pca.res <- prcomp(X.train, scale. = TRUE, center = TRUE)
res.new <- as.matrix(scale(X.test, scale = pca.res$scale, 
                          center = pca.res$center)) %*% pca.res$rotation
