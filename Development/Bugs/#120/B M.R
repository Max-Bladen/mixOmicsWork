library(devtools)
setwd("C:/Users/Work/Desktop/bladen-devel/R")

load_all()

data(nutrimouse)
diablo <- block.splsda(X = list(gene = nutrimouse$gene, 
                                lipid = nutrimouse$lipid), 
                       Y = nutrimouse$diet,
                       ncomp = 3)

nutri.spls <- spls(nutrimouse$gene, nutrimouse$lipid)
nutri.pca <- pca(nutrimouse$gene)
nutri.splsda <- splsda(nutrimouse$gene, nutrimouse$diet)

# ---------------------------------------------------------------------------- #
# ============================================================================ #
# ---------------------------------------------------------------------------- #

loadings.spls <- plotLoadings(nutri.spls)
loadings.pca <- plotLoadings(nutri.pca)
loadings.splsda <- plotLoadings(nutri.splsda)


loadings.diablo1 <- plotLoadings(diablo, comp = 1)
loadings.diablo3 <- plotLoadings(diablo, comp = 3)

