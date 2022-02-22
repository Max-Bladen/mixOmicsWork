library(mixOmics)

data("stemcells")

X <- stemcells$gene
Y <- stemcells$celltype
S <- stemcells$study

optimal.ncomp <- 2
optimal.keepX <- c(24, 45)

splsda.stemcells <- mint.splsda(X = X, Y = Y, 
                                  study = S, 
                                  ncomp = optimal.ncomp, 
                                  keepX = optimal.keepX)
plotArrow(splsda.stemcells)



data(srbct)
splsda.srbct <- splsda(srbct$gene, srbct$class, ncomp = 3, keepX = c(9,280,30))
plotArrow(splsda.srbct)