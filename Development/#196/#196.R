# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)
library(readxl)

# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")

rm(list = c(".perf.mixo_pls_cv", "perf", "perf.mixo_pls", "perf.mixo_plsda", "perf.mixo_spls", "perf.mixo_splsda"))
load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

interleuch <- read_xlsx("C:/Users/Work/Desktop/mO Work/Support Files/Leandro/Interleuch.xlsx")
SCFA <- read_xlsx("C:/Users/Work/Desktop/mO Work/Support Files/Leandro/SCFA.xlsx")

data("liver.toxicity")

# ---------------------------------------------------------------------------- #
# =========================== CASES WITH NO ERROR ============================ #
# ---------------------------------------------------------------------------- #

X <- liver.toxicity$gene#[1:24, 1:24]
Y <- liver.toxicity$clinic#[1:24, 1:9]

MyResult.pls1 <- pls(Y, X, mode = "canonical", ncomp = 4)
plotIndiv(MyResult.pls1)
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
perf.pls <- perf(MyResult.pls1, validation = "Mfold", folds = 4, progressBar = T, nrepeat = 3)
plot(perf.pls, criterion = "Q2.total")

# ---------------------------------------------------------------------------- #
# ============================= CASES WITH ERROR ============================= #
# ---------------------------------------------------------------------------- #

X <- data.frame(interleuch)
rownames(X) <- X[,1]
X <- X[, -1]

Y <- data.frame(SCFA)
rownames(Y) <- Y[,1]
Y <- Y[, -1]

MyResult.pls1 <- pls(Y, X, mode = "classic", ncomp = 2)
#plotIndiv(MyResult.pls1)
set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
perf.pls <- perf(MyResult.pls1, validation = "Mfold", folds = 4, progressBar = FALSE, nrepeat = 3)
plot(perf.pls, criterion = "Q2.total")


Yc <- c(rep(1,12), rep(2,12))
MyResult.plsda <- plsda(X, Yc, ncomp=2)
set.seed(30)
perf.plsda <- perf(MyResult.plsda, validation = "Mfold", folds = 4, progressBar = FALSE, nrepeat = 3)
plot(perf.plsda)

# ---------------------------------------------------------------------------- #
# ================================== REPREX ================================== #
# ---------------------------------------------------------------------------- #

library(mixOmics)

data("liver.toxicity")

# reducing number of features to reduce run time
X <- liver.toxicity$gene[, 1:1000]
Y <- liver.toxicity$clinic

# to reproduce error, we need to induce some features to have near zero variance
X[, c(1, 23, 62, 234, 789)] <- 0

pls.obg <- pls(Y, X, ncomp = 4)
pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 4, 
                     progressBar = F, 
                     nrepeat = 3)