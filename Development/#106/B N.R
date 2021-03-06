library(devtools)
library(mixOmics)
setwd("~/GitHub/mixOmics/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================ REPRODUCING ERROR ============================= #
# ---------------------------------------------------------------------------- #

# load data in
data("nutrimouse")
data("liver.toxicity")
data("breast.TCGA")
data("stemcells")

# ---------------------------------------------------------------------------- #

nutri.pls <- pls(nutrimouse$gene, nutrimouse$lipid)
cutoff_00 <- network(nutri.pls)
cutoff_05 <- network(nutri.pls, cutoff = 0.5)

cutoff_00$cutoff

# ---------------------------------------------------------------------------- #
# ================================== CHECKS ================================== #
# ---------------------------------------------------------------------------- #


