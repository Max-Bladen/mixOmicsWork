# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)

# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/bladen-devel/R")

load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic



# ---------------------------------------------------------------------------- #
# ================================== TESTS =================================== #
# ---------------------------------------------------------------------------- #



liver.pls <- pls(X, Y, ncomp = 5, mode = "canonical")
liver.pls.val <- perf(liver.pls, validation = "Mfold", folds = 5, nrepeat = 3)

#names(liver.pls.val)
#liver.pls.val$features
names(liver.pls.val$measures)

# ---------------------------------------------------------------------------- #

liver.spls <- spls(X, Y, ncomp = 5, keepX = rep(5,5), keepY = rep(10, 5))
liver.spls.val <- perf(liver.spls, validation = "Mfold", folds = 5, nrepeat = 3)

names(liver.spls.val$features)
names(liver.spls.val$features$stability.X)
head(liver.spls.val$features$stability.X$comp1)
names(liver.pls.val$measures)
