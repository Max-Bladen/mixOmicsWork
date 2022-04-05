library(devtools)
library(mixOmics)
setwd("C:/Users/Work/Desktop/mixOmics/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================ REPRODUCING ERROR ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)
data(liver.toxicity)
liver.pca <- pca(liver.toxicity$gene, ncomp=2) 
cim(liver.pca) 

# ---------------------------------------------------------------------------- #

# control case
liver.spca <- spca(liver.toxicity$gene, 
                   ncomp=2)#, keepX = c(10,10))
# test case
liver.pca <- pca(liver.toxicity$gene, 
                 ncomp=2) 

# ---------------------------------------------------------------------------- #

# control cases
cim(liver.spca) 

# test case
cim(liver.pca) 

# ---------------------------------------------------------------------------- #
# ================================== CHECKS ================================== #
# ---------------------------------------------------------------------------- #

plotIndiv(liver.pca)
plotIndiv(liver.spca)

plotVar(liver.pca)
plotVar(liver.spca)

cim(liver.pca)
cim(liver.spca)

biplot(liver.pca)
biplot(liver.spca)
