library(devtools)
library(mixOmics)
library(BiocParallel)
setwd("C:/Users/Work/Desktop/Dev Build/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================ REPRODUCING ERROR ============================= #
# ---------------------------------------------------------------------------- #
# Set up data
data("breast.TCGA")

# set a list of all the X dataframes
X = list(miRNA = breast.TCGA$data.train$mirna, 
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype

design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

# ---------------------------------------------------------------------------- #

BPPARAM <- BiocParallel::SnowParam(workers = 4)


noParaTimes <- c()
for (i in 1:1) {
  start = Sys.time()
  
  diablo.tuning <- tune.block.splsda(X, Y, design = design, ncomp = 1, nrepeat = 1,
                                     test.keepX = list(miRNA = 1:3, mRNA = 1:3, proteomics = 1:3),
                                     progressBar = TRUE) # 
  end = Sys.time()
  
  noParaTimes <- c(noParaTimes, end-start)
}

fullParaTimes <- c()
for (i in 1:1) {
  start = Sys.time()
  
  diablo.tuning <- tune.block.splsda(X, Y, design = design, ncomp = 1, nrepeat = 1,
                                     test.keepX = list(miRNA = 1:3, mRNA = 1:3, proteomics = 1:3),
                                     progressBar = TRUE, BPPARAM = BPPARAM) # 
  end = Sys.time()
  
  fullParaTimes <- c(fullParaTimes, end-start)
}



# ---------------------------------------------------------------------------- #
# ================================== CHECKS ================================== #
# ---------------------------------------------------------------------------- #