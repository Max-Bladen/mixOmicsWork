# FR.K - background predcit for MINT sPLS-DA

library(mixOmics) # import the mixOmics library


## -------------------------------------------------------------------------------------------------------------------
data(stemcells) # extract stem cells data

X <- stemcells$gene # use genetic expression levels as predictors
Y <- stemcells$celltype # use stem cell type as response
study <- stemcells$study # extract study allocation of samples

basic.plsda.model <- mint.plsda(X, Y, study = study, ncomp = 2) # generate basic MINT pls-da model

plotIndiv(basic.plsda.model, legend = TRUE, title = ' ', subtitle = ' ', ellipse = TRUE) # plot the samples

bgp <- background.predict(basic.splsda.model, comp.predicted = 2, dist = "max.dist")

## -------------------------------------------------------------------------------------------------------------------

data(stemcells)
res = mint.splsda(X = stemcells$gene, Y = stemcells$celltype, ncomp = 3,
                  study = stemcells$study)
## before background
plotIndiv(res)

background = background.predict(res, comp.predicted = 2, dist = "mahalanobis.dist")
## after background
plotIndiv(res, background = background, legend = TRUE )