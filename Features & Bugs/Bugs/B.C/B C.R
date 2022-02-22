library(devtools)
setwd("~/GitHub/mixOmics/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================ REPRODUCING ERROR ============================= #
# ---------------------------------------------------------------------------- #

data(breast.TCGA)

breast.data1 <- list(mirna = breast.TCGA$data.train$mirna,
                    mrna = breast.TCGA$data.train$mrna,
                    protein = breast.TCGA$data.train$protein)

breast.block.spls.comb <- block.spls(breast.data1, indY = 3,
                                     keepX = list(mirna = c(10,10),
                                                  mrna = c(10,10)))

breast.data2 <- list(mirna = breast.TCGA$data.train$mirna,
                    mrna = breast.TCGA$data.train$mrna)

breast.block.spls.sep <- block.spls(breast.data2, breast.TCGA$data.train$protein,
                                    keepX = list(mirna = c(10,10),
                                                 mrna = c(10,10)))

breast.group = breast.TCGA$data.train$subtype

circosPlot(breast.block.spls.sep, group = breast.group, 
           cutoff = 0.8, Y.name = "protein") # works, can find object$X$Y

# ---------------------------------------------------------------------------- #

load_all()
breast.block.spls.comb <- block.spls(breast.data1, indY = 3, ncomp = 2,
                                     keepX = list(mirna = c(10,10),
                                                  mrna = c(10,10)))
circosPlot(breast.block.spls.comb, group = breast.group, 
           cutoff = 0.8) # doesnt work


# ---------------------------------------------------------------------------- #
# ================================== CHECKS ================================== #
# ---------------------------------------------------------------------------- #


plotIndiv(breast.block.spls.comb)
plotVar(breast.block.spls.comb)
network(breast.block.spls.comb)
cim(breast.block.spls.comb)
