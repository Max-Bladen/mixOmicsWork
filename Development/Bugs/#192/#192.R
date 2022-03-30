# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")

load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

suppressMessages(library(mixOmics))

data(breast.TCGA) # load in the data

# extract data
X.train = list(mirna = breast.TCGA$data.train$mirna,
               mrna = breast.TCGA$data.train$mrna)

X.test = list(mirna = breast.TCGA$data.test$mirna,
              mrna = breast.TCGA$data.test$mrna)

Y.train = breast.TCGA$data.train$subtype

# use optimal values from the case study on mixOmics.org
optimal.ncomp = 2
optimal.keepX = list(mirna = c(10,5),
                     mrna = c(26, 16))

# set design matrix
design = matrix(0.1, ncol = length(X.train), nrow = length(X.train),
                dimnames = list(names(X.train), names(X.train)))
diag(design) = 0

# generate model
final.diablo.model = block.splsda(X = X.train, Y = Y.train, ncomp = optimal.ncomp, # set the optimised DIABLO model
                                  keepX = optimal.keepX, design = design)


# create new test data with one dataframe being reordered
new.var.order = sample(1:dim(X.test$mirna)[2])
X.test.dup <- X.test
X.test.dup$mirna <- X.test.dup$mirna[, new.var.order]

predict.diablo = predict(final.diablo.model, newdata = X.test)

predict.diablo.reordered = predict(final.diablo.model, newdata = X.test.dup)

# ---------------------------------------------------------------------------- #
# =========================== CASES WITH NO ERROR ============================ #
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# ============================= CASES WITH ERROR ============================= #
# ---------------------------------------------------------------------------- #

# create new test data with one dataframe being reordered
new.var.order = sample(1:dim(X.test$mirna)[2])

X.test.dup <- X.test
X.test.dup$mirna <- X.test.dup$mirna[, new.var.order]

predict.diablo = predict(final.diablo.model, newdata = X.test)
predict.diablo.reordered = predict(final.diablo.model, newdata = X.test.dup)

homogenity <- matrix(NA, nrow = 2, ncol = 3)
colnames(homogenity) <- c("max.dist", "centroids.dist", "mahalanobis.dist")
rownames(homogenity) <- c("mirna", "mrna")

for (dist in colnames(homogenity)) {
  for (block in rownames(homogenity)) {
    homogenity[block, dist] = all(predict.diablo$class[[dist]][[block]] == predict.diablo.reordered$class[[dist]][[block]])
  }
}
homogenity