library(devtools)
library(mixOmics)
setwd("C:/Users/Work/Desktop/Dev Build/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================== LIVER TOXICITY ============================== #
# ---------------------------------------------------------------------------- #

data(liver.toxicity)

grid.keepX<-c(seq(5,30,5))

BPPARAM <- BiocParallel::SnowParam(workers = 15)

# ---------------------------------------------------------------------------- #
# Run with no NAs
X <- liver.toxicity$gene

tune.spca.result.full <- tune.spca(X, ncomp=3, folds=3, nrepeat=50,
                                   test.keepX = grid.keepX,
                                   BPPARAM = BPPARAM)

# ---------------------------------------------------------------------------- #
# Run with NAs
X <- liver.toxicity$gene

for (r in sample(nrow(X), floor(0.15*nrow(X)))) {
  X[r, sample(ncol(X), 1)] <- NA
}

tune.spca.result.NAs <- tune.spca(X, ncomp=3, folds=3, nrepeat=50,
                                  test.keepX = grid.keepX,
                                  BPPARAM = BPPARAM)

# ---------------------------------------------------------------------------- #

tune.spca.result.full$cor.comp
tune.spca.result.NAs$cor.comp

plot(tune.spca.result.full)
plot(tune.spca.result.NAs)

# ---------------------------------------------------------------------------- #
# ================================= MULTIDRUG ================================ #
# ---------------------------------------------------------------------------- #

data(multidrug)

grid.keepX<-c(seq(5,30,5))

# ---------------------------------------------------------------------------- #
# Run with no NAs
X <- multidrug$ABC.trans

tune.spca.result.full <- tune.spca(X, ncomp=3, folds=3, nrepeat=50,
                                   test.keepX = grid.keepX,
                                   BPPARAM = BPPARAM)

# ---------------------------------------------------------------------------- #
# Run with NAs
X <- multidrug$ABC.trans

for (r in sample(nrow(X), floor(0.15*nrow(X)))) {
  X[r, sample(ncol(X), 1)] <- NA
}

tune.spca.result.NAs <- tune.spca(X, ncomp=3, folds=3, nrepeat=50,
                                  test.keepX = grid.keepX,
                                  BPPARAM = BPPARAM)

# ---------------------------------------------------------------------------- #

tune.spca.result.full$cor.comp
tune.spca.result.NAs$cor.comp

plot(tune.spca.result.full)
plot(tune.spca.result.NAs)

# ---------------------------------------------------------------------------- #
# ================================= NUTRIMOUSE =============================== #
# ---------------------------------------------------------------------------- #

suppressMessages(library(mixOmics))
data(nutrimouse)

grid.keepX<-c(seq(5,30,5)) # recreate grid from issue #161

## run tuning with no NAs
X <- nutrimouse$gene
tune.spca.result.full <- tune.spca(X, ncomp=3, folds=3, nrepeat=3,
                                   test.keepX = grid.keepX)

## introduce a single NA to ~15% of rows in the data frame
for (r in sample(nrow(X), floor(0.15*nrow(X)))) {
  X[r, sample(ncol(X), 1)] <- NA
}

## run tuning with NAs
tune.spca.result.NAs <- tune.spca(X, ncomp=3, folds=3, nrepeat=3,
                                  test.keepX = grid.keepX)

# ---------------------------------------------------------------------------- #

tune.spca.result.full$cor.comp
tune.spca.result.NAs$cor.comp

plot(tune.spca.result.full)
plot(tune.spca.result.NAs)

# ---------------------------------------------------------------------------- #
# =============================== TESTING NIPALS ============================= #
# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #

data(liver.toxicity)

grid.keepX <- c(seq(5,30,5))

BPPARAM <- BiocParallel::SnowParam(workers = 15)

# ---------------------------------------------------------------------------- #
# Run with no NAs
X <- liver.toxicity$gene

tune.spca.result.full <- tune.spca.nipals(X, ncomp=3, folds=3, nrepeat=50,
                                   test.keepX = grid.keepX,
                                   BPPARAM = BPPARAM,
                                   impute.NAs = T)

# ---------------------------------------------------------------------------- #
# Run with NAs
X <- liver.toxicity$gene#[, 1:100]

for (r in sample(nrow(X), floor(0.15*nrow(X)))) {
  X[r, sample(ncol(X), 1)] <- NA
}

tune.spca.result.NAs <- tune.spca.nipals(X, ncomp=3, folds=3, nrepeat=3,
                                  test.keepX = grid.keepX,
                                  BPPARAM = BPPARAM,
                                  impute.NAs = T)

# ---------------------------------------------------------------------------- #

tune.spca.result.full$cor.comp
tune.spca.result.NAs$cor.comp

plot(tune.spca.result.full)
plot(tune.spca.result.NAs)
