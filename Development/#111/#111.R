# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)

# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")


rm(list = c("Check.entry.pls", "Check.entry.rgcca", "Check.entry.sgcca", "Check.entry.single", "Check.entry.wrapper.mint.block", "get.keepA", "block.plsda"))
devtools::load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

data(srbct)
X <- srbct$gene
Y <- srbct$class

data("breast.TCGA")
X.b <- list(mirna = breast.TCGA$data.train$mirna,
            mrna = breast.TCGA$data.train$mrna)
Y.b <- breast.TCGA$data.train$subtype
Y.b.2 <- breast.TCGA$data.train$protein

data("stemcells")
X.m <- stemcells$gene
Y.m <- stemcells$celltype
S.m <- stemcells$study


Y[c(1,2,3)] <- NA
Y.b[c(1,2,3)] <- NA
Y.m[c(1,2,3)] <- NA

# ---------------------------------------------------------------------------- #
# ================================== TESTS =================================== #
# ---------------------------------------------------------------------------- #

# plsda
obj <- plsda(X, Y)

# splsda
obj <- splsda(X, Y)

# block.plsda
obj <- block.plsda(X.b, Y.b)

# block.splsda
obj <- block.splsda(X.b, Y.b)

# mint.plsda
obj <- mint.plsda(X.m, Y.m, study = S.m)

# mint.splsda
obj <- mint.splsda(X.m, Y.m, study = S.m)

# block.pls
obj <- block.pls(X.b, Y.b.2)

# block.spls
obj <- block.spls(X.b, Y.b.2)

# mint.pls
obj <- mint.pls(X.m, X.m, study = S.m)

# mint.spls
obj <- mint.spls(X.m, X.m, study = S.m)


#### 
data(nutrimouse)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = nutrimouse$diet)
data$Y[c(1,2,3)] <- NA
res = block.plsda(X = data, indY = 3)

MyResult.splsda <- splsda(X, Y, keepX = c(50,50))


