library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")
load_all()
# ---------------------------------------------------------------------------- #

library(mixOmics)

# ---------------------------------------------------------------------------- #

data("stemcells")
X <- stemcells$gene # use genetic expression levels as predictors
Y <- stemcells$celltype # use stem cell type as response
study <- stemcells$study # extract study allocation of samples

# ---------------------------------------------------------------------------- #

stem.mint.tune <- tune.mint.splsda(X = X, Y = Y, study = study,  # tune the number of features
                       ncomp = 2,
                       test.keepX = seq(1, 100, 10),
                       measure = 'BER', # balanced error rate
                       dist = "centroids.dist")
plot(stem.mint.tune)

# ---------------------------------------------------------------------------- #

stem.mint.tune <- tune.mint.splsda(X = X, Y = Y, study = study,  # tune the number of features
                                   ncomp = 2,
                                   test.keepX = seq(1, 100, 10),
                                   measure = 'BER', # balanced error rate
                                   dist = "centroids.dist")
plot(stem.mint.tune)


stem.tune <- tune.splsda(X = X, Y = Y,  # tune the number of features
                                   ncomp = 2,
                                   test.keepX = seq(1, 100, 10),
                                   measure = 'BER', # balanced error rate
                                   dist = "centroids.dist") # 
plot(stem.tune) # 
