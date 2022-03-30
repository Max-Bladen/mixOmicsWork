library(mixOmics)
library("BiocParallel")

#Read RDS file
X_new <- readRDS("X_new.RData")
X_new1 <- readRDS("X_new1.RData")

Y <- read.csv("Y.csv")
#Converting Outcome variable to a Factor
Y3 <- as.factor(Y$Y3) #PBF at 1-year [8]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******Parameter choice: Design*******
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

design <- matrix(0.1, ncol = length(X_new1), nrow = length(X_new1), 
                 dimnames = list(names(X_new1), names(X_new1)))
diag(design) <- 0

design 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******Tuning the number of components*******
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sgccda.res <- block.splsda(X = X_new1, Y = Y3, ncomp = 3, 
                           design = design)

set.seed(155) # for reproducibility, only when the `cpus' argument is not used
perf.diablo <- perf(sgccda.res, validation = 'Mfold', folds = 5, nrepeat = 50,
                    cpus = 5, progressBar = TRUE)

#perf.diablo  # lists the different outputs
plot(perf.diablo) 

#Considering this distance and the BER, the output $choice.ncomp indicates an optimal number of components.
perf.diablo$choice.ncomp$WeightedVote

ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******Tuning keepX - number of variables*******
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#test.keepX <- list (asv_clr = c(5:9, seq(10, 18, 2), seq(20,30,5)),
#                    metablog = c(5:9, seq(10, 18, 2), seq(20,30,5)))

test.keepX = list (asv_clr = seq(5, 25, 2),
                   metablog =  seq(5, 25, 2))

tune.TCGA_Y3 <- tune.block.splsda(X = X_new, Y = Y3, ncomp = ncomp, 
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = 5, nrepeat = 20,
                                      dist = "centroids.dist", scale = FALSE, progressBar = TRUE,
                                  BPPARAM = MulticoreParam())

list.keepX <- tune.TCGA_Y3$choice.keepX
list.keepX

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*******Final model: The final DIABLO model is run as:*******
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sgccda.res <- block.splsda(X = X_new1, Y = Y3, ncomp = ncomp, 
                           keepX = list.keepX, design = design)

#sgccda.res   # list the different functions of interest related to that object
sgccda.res$design

selectVar(sgccda.res, block = 'asv_clr', comp = 1)$asv_clr$name 
selectVar(sgccda.res, block = 'metablog', comp = 1)$metablog$name 

plotDiablo(sgccda.res, ncomp = 1)

plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(0.5,0.5), col = c('red', 'grey3'))

circosPlot(sgccda.res, cutoff = 0.5, line = TRUE, 
           color.blocks= c('brown1', 'lightgreen'),
           color.cor = c("blue","grey20"), size.variables = 0.75, size.labels = 1.5, comp = 1)

network.res <- network(sgccda.res, save='jpeg', name.save = 'mynetwork1', cutoff = 0.5, cex.node.name=0.6)

plotLoadings(sgccda.res, comp = 1, contrib = 'max', method = 'median', size.name = 0.85, size.subtitle = rel(1))

cim.res <- cimDiablo(sgccda.res, margins=c(12,12), size.legend = 0.5, legend.position = "topright", save = 'jpeg', name.save = 'PLS_CIM_image1')