# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)

# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")

rm(list = c("background.predict", "internal_graphicModule", "mint.splsda", "internal_wrapper.mint"))
load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

data(srbct)
data(stemcells)

# ---------------------------------------------------------------------------- #
# =========================== CASES WITH NO ERROR ============================ #
# ---------------------------------------------------------------------------- #

stnd.res.srbct <- splsda(X = srbct$gene, Y = srbct$class, ncomp = 3)

bg.srbct = background.predict(stnd.res.srbct, comp.predicted = 2, 
                              dist = "centroids.dist")

plotIndiv(stnd.res.srbct, background = bg.srbct, legend = T, style = "ggplot2")

# ---------------------------------------------------------------------------- #
# DOESN'T PLOT THE POLYGONS PROPERLY
stnd.res.stem <- splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 4)

bg.stem = background.predict(stnd.res.stem, comp.predicted = 2, 
                             dist = "centroids.dist")

plotIndiv(stnd.res.stem, background = bg.stem, legend = T, style = "ggplot2")

# ---------------------------------------------------------------------------- #
# ============================= CASES WITH ERROR ============================= #
# ---------------------------------------------------------------------------- #

# mint.res.stem$ind.mat <- scale(mint.res.stem$ind.mat)

plot.new()
plot(c(mint.bg.stem.max$hiPSC[,1], mint.bg.stem.max$hESC[,1], mint.bg.stem.max$Fibroblast[,1]),
     c(mint.bg.stem.max$hiPSC[,2], mint.bg.stem.max$hESC[,2], mint.bg.stem.max$Fibroblast[,2]))



mint.res.stem <- mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 4,
                             study = stemcells$study)

mint.bg.stem.max = background.predict(mint.res.stem, comp.predicted = 2, # raises error found in #24 on GitHub
                                      dist = "max.dist")

mint.bg.stem.cen = background.predict(mint.res.stem, comp.predicted = 2, # no error but the background does not display at all
                             dist = "centroids.dist")

mint.bg.stem.mah = background.predict(mint.res.stem, comp.predicted = 2, # no error but the background does not display at all
                                      dist = "mahalanobis.dist")

plotIndiv(mint.res.stem, background = mint.bg.stem.mah, legend = T, style = "ggplot2")

# ---------------------------------------------------------------------------- #
# note when passing in an object of one type (eg mint.splsda) and then using the
# background object which inherits from another tpye (eg. splsda) then no error
# will be raised as well as no background will be shown
bg = background.predict(stnd.res.srbct, comp.predicted = 2, dist = "max.dist")

plotIndiv(stnd.res.stem, background = bg, legend = TRUE)

# ---------------------------------------------------------------------------- #
