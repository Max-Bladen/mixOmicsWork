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

data(stemcells)
X = stemcells$gene
Y = stemcells$celltype
study <- stemcells$study
# data("srbct")
# X <- srbct$gene
# Y <- srbct$class
# study <- sample(stemcells$study)[1:63]


# ---------------------------------------------------------------------------- #
# =========================== CASES WITH NO ERROR ============================ #
# ---------------------------------------------------------------------------- #

# rm(list = c("perf.mint.pls", "perf.mint.plsda", "perf.mint.spls", "perf.mint.splsda", "LOGOCV"))
# load_all()

# ---------------------------------------------------------------------------- #
# ============================= CASES WITH ERROR ============================= #
# ---------------------------------------------------------------------------- #


# c('max.dist', 'centroids.dist', 'mahalanobis.dist')
metricsC1 <- matrix(0, nrow = 2, ncol = 3)
colnames(metricsC1) <- c('max.dist', 'centroids.dist', 'mahalanobis.dist')
rownames(metricsC1) <- c("overall", "BER")

metricsC2 <- metricsC1

for (dist in c('max.dist', 'centroids.dist', 'mahalanobis.dist')) {
  for (measure in c("overall", "BER")) {
    tune.mint = tune.mint.splsda(X = X, Y = Y, study = study, ncomp = 2, test.keepX = seq(1, 51, 5),
                                 dist = dist, progressBar = FALSE, measure = measure)
    
    mint.splsda.res = mint.splsda(X = X, Y = Y, study = study, ncomp = 2,
                                  keepX = tune.mint$choice.keepX)
    
    perf.mint = perf(mint.splsda.res, progressBar = FALSE, dist = dist)
    
    metricsC1[measure, dist] <- tune.mint$error.rate[which(rownames(tune.mint$error.rate) == tune.mint$choice.keepX[1]), "comp1"]
    metricsC2[measure, dist] <- tune.mint$error.rate[which(rownames(tune.mint$error.rate) == tune.mint$choice.keepX[2]), "comp2"]
    
    if (round(perf.mint$global.error[[measure]]["comp1",], 4) != round(metricsC1[[measure, dist]], 4)) {
      cat("DIFF ON COMP 1, ", measure, dist, "\n")
      cat("perf: ", perf.mint$global.error[[measure]]["comp1",], "\n")
      cat("tune: ", metricsC1[measure, dist], "\n")
      }
    if (round(perf.mint$global.error[[measure]]["comp2",], 4) != round(metricsC2[[measure, dist]], 4)) {
      cat("DIFF ON COMP 2, ", measure, dist, "\n")
      cat("perf: ", perf.mint$global.error[[measure]]["comp2",], "\n")
      cat("tune: ", metricsC2[measure, dist], "\n")
      }
  }
}

