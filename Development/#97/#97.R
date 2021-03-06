# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)

# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")

load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

data(multidrug)
X <- multidrug$ABC.trans
dim(X) # check dimension of data

# ---------------------------------------------------------------------------- #
# =========================== CASES WITH NO ERROR ============================ #
# ---------------------------------------------------------------------------- #

pca.result <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)

# ---------------------------------------------------------------------------- #
# ============================= CASES WITH ERROR ============================= #
# ---------------------------------------------------------------------------- #

cell.line = as.factor(multidrug$cell.line$Class)
biplot(pca.result, group = cell.line, cutoff = c(0.99, 0.1))

