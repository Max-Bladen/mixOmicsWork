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

rm(list = c("MCVfold.spls"))

# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# =========================== CASES WITH NO ERROR ============================ #
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# ============================= CASES WITH ERROR ============================= #
# ---------------------------------------------------------------------------- #

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$treatment$Dose.Group
# create a class with one sample only
Y[c(1,2)] <- 'FOO'

liver.plsda <- plsda(X, Y, ncomp = 2)
set.seed(12)
liver.val <- perf(liver.plsda, validation = "Mfold", folds = 3, nrepeat = 3)

# ---------------------------------------------------------------------------- #

data(nutrimouse)
Y = nutrimouse$diet
# create a class with one sample only
Y = as.character(Y)
Y[1] <- 'FOO'
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
nutrimouse.sgccda <- block.splsda(X=data,
                                  Y = Y,
                                  keepX = list(gene=c(10,10), lipid=c(15,15)),
                                  ncomp = 2,
                                  scheme = "horst")

perf.nutri = perf(nutrimouse.sgccda)
