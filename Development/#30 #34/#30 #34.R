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

data("nutrimouse")
data1 <- nutrimouse$gene
data2 <- nutrimouse$lipid

data("stemcells")
mint.X <- stemcells$gene[1:62,]
mint.Y <- stemcells$gene[63:124,]
study <- stemcells$study[1:62]
mint.class <- stemcells$celltype[1:62]

# ---------------------------------------------------------------------------- #
# ================================== TESTS =================================== #
# ---------------------------------------------------------------------------- #

optimal.ncomp <- 5

obj <- pca(X=data1, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp

obj <- spca(X=data1, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp

obj <- pls(X=data1, Y=data2, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp

obj <- spls(X=data1, Y=data2, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp

obj <- rcc(X=data1, Y=data2, ncomp = optimal.ncomp, verbose.call=T, method = "shrinkage")
obj$call$simple.call
obj$call$ncomp


obj <- block.pls(X=list(data1 = data1,
                  data2 = data2), indY=2, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp


obj <- block.plsda(X=list(data1 = data1,
                        data2 = data2), Y = c(rep(1,20), rep(2,20)), ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp


obj <- block.spls(X=list(data1 = data1,
                        data2 = data2), indY=2, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp


obj <- block.splsda(X=list(data1 = data1,
                          data2 = data2), Y = c(rep(1,20), rep(2,20)), ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp


obj <- mint.pca(X=mint.X, study = study, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp

obj <- mint.pls(X=mint.X, Y = mint.Y, study=study, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp


# ---------------------------------------------------------------------------- #

data("nutrimouse")
data1 <- nutrimouse$gene
data2 <- nutrimouse$lipid

optimal.ncomp <- 5

obj <- pca(X=data1, ncomp = optimal.ncomp, verbose.call=T)
obj$call$simple.call
obj$call$ncomp

# ---------------------------------------------------------------------------- #

"

@template arg/verbose.call

\item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
  #' If \code{verbose.call = TRUE} then all the inputted values are accessable via
  #' this component}


,
verbose.call = FALSE

"

if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }