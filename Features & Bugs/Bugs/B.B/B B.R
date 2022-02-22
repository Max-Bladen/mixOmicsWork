library(devtools)
setwd("~/GitHub/mixOmics/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================ REPRODUCING ERROR ============================= #
# ---------------------------------------------------------------------------- #

# load data in
data(nutrimouse)
data(breast.TCGA)

# set dataframes to be analysed
nutri.data <- list(gene = nutrimouse$gene,
                   lpipd = nutrimouse$lipid)

breast.data <- list(mirna = breast.TCGA$data.train$mirna,
                    mrna = breast.TCGA$data.train$mrna,
                    protein = breast.TCGA$data.train$protein)
breast.outcome <- breast.TCGA$data.train$subtype

# generate various models
nutri.pls <- pls(nutrimouse$gene, nutrimouse$lipid, ncomp = 3)
breast.pls <- pls(breast.TCGA$data.train$mirna, breast.TCGA$data.train$mrna)

nutri.cca <- rcc(nutrimouse$gene, nutrimouse$lipid, method = 'shrinkage')
breast.cca <- rcc(breast.TCGA$data.train$mirna, breast.TCGA$data.train$mrna, 
                  method = "ridge", lambda1 = 0.1, lambda2 = 0.1)

breast.block.pls <- block.pls(breast.data, indY = 3, ncomp = 3)
breast.diablo <- block.plsda(breast.data, breast.outcome)

# run plotArrow

plotArrow(nutri.pls, ind.names = TRUE, ind.names.position = "end", 
          group = breast.outcome)











