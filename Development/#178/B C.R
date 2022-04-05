library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ============================ REPRODUCING ERROR ============================= #
# ---------------------------------------------------------------------------- #

### --- SET UP DATA --- ###

suppressMessages(library(mixOmics)) # load package and data
data(breast.TCGA)

## class membership
breast.group = breast.TCGA$data.train$subtype

## data1 has all dataframes
breast.data1 <- list(mirna = breast.TCGA$data.train$mirna,
                    mrna = breast.TCGA$data.train$mrna,
                    protein = breast.TCGA$data.train$protein)

## data2 has all but protein data
breast.data2 <- list(mirna = breast.TCGA$data.train$mirna,
                    mrna = breast.TCGA$data.train$mrna)

### --- RUN BLOCK.SPLS --- ###

## form block.spls with combined dataframes
breast.block.spls.comb <- block.spls(breast.data1, indY = 3,
                                     keepX = list(mirna = c(10,10),
                                                  mrna = c(10,10)))

## form block.spls with separated dataframes
breast.block.spls.sep <- block.spls(breast.data2, breast.TCGA$data.train$protein,
                                    keepX = list(mirna = c(10,10),
                                                 mrna = c(10,10)))

### --- PASS TO CIRCOS.PLOT --- ###

## works, can find object$X$Y
circosPlot(breast.block.spls.sep, group = breast.group, 
           cutoff = 0.8) 

## ERROR RAISED !
circosPlot(breast.block.spls.comb, group = breast.group, 
           cutoff = 0.8) # doesnt work


# ---------------------------------------------------------------------------- #
# ================================== CHECKS ================================== #
# ---------------------------------------------------------------------------- #


plotIndiv(breast.block.spls.comb)
plotVar(breast.block.spls.comb)
network(breast.block.spls.comb)
cim(breast.block.spls.comb)
