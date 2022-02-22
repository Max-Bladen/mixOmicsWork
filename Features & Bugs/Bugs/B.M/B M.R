library(devtools)
setwd("~/GitHub/mixOmics/R")

data(nutrimouse)
diablo <- block.splsda(X = list(gene = nutrimouse$gene, 
                                lipid = nutrimouse$lipid), 
                       Y = nutrimouse$diet,
                       ncomp = 3)

nutri.spls <- spls(nutrimouse$gene, nutrimouse$lipid)

# ---------------------------------------------------------------------------- #
# ============================================================================ #
# ---------------------------------------------------------------------------- #

load_all()
loadings1 <- plotLoadings(diablo, comp = 1)
loadings3 <- plotLoadings(diablo, comp = 3)
loadings1
loadings3

