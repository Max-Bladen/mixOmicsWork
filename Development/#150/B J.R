library(devtools)
library(mixOmics)
setwd("C:/Users/Work/Desktop/mixOmics/R")

# ---------------------------------------------------------------------------- #

load_all()

# ---------------------------------------------------------------------------- #
# ================================ STEM CELLS ================================ #
# ---------------------------------------------------------------------------- #

data("stemcells")

X <- stemcells$gene
Y <- stemcells$celltype
S <- stemcells$study

optimal.ncomp <- 4
optimal.keepX <- c(24, 45, 20, 30)

# ----------------------------- Standard Framework --------------------------- #

splsda.stemcells <- splsda(X = X, Y = Y, 
                                ncomp = optimal.ncomp, 
                                keepX = optimal.keepX)
plotArrow(splsda.stemcells) # produces error

# ------------------------------- MINT Framework ----------------------------- #

mint.splsda.stemcells <- mint.splsda(X = X, Y = Y, 
                                study = S, 
                                ncomp = optimal.ncomp, 
                                keepX = optimal.keepX)

#load_all()
plotArrow(mint.splsda.stemcells) # produces error

# ---------------------------------------------------------------------------- #
# =============================== BREAST SRBCT =============================== #
# ---------------------------------------------------------------------------- #

data(srbct)

X <- srbct$gene
Y <- srbct$class

optimal.ncomp <- 3
optimal.keepX = c(9,280,30)

# ----------------------------- Standard Framework --------------------------- #

splsda.srbct <- splsda(X = X, Y = Y,
                            ncomp = optimal.ncomp, 
                            keepX = optimal.keepX)
plotArrow(splsda.srbct) # produces error




# ---------------------------------------------------------------------------- #
# =============================== BREAST TCGA ================================ #
# ---------------------------------------------------------------------------- #

data(breast.TCGA) # load in the data

# set a list of all the X dataframes
X <- list(miRNA = breast.TCGA$data.train$mirna, 
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype # set the response variable as the Y df

optimal.ncomp <- 2
optimal.keepX <-  list(miRNA = c(10, 5), 
                       mRNA = c(25,16),
                       proteomics = c(8,5))

# ------------------------------ DIABLO Framework ---------------------------- #

tcga.diablo <- block.splsda(X, Y,
                            ncomp = optimal.ncomp,
                            keepX = optimal.keepX)

#load_all()
plotArrow(tcga.diablo) # DOESN'T produce error




# ---------------------------------------------------------------------------- #
# ============================== LIVER TOXICITY ============================== #
# ---------------------------------------------------------------------------- #

data(liver.toxicity) # extract the liver toxicity data
X <- liver.toxicity$gene # use the gene expression data as the X matrix
Y <- liver.toxicity$clinic # use the clinical data as the Y matrix

optimal.ncomp <- 2
optimal.keepX = c(35,45)
optimal.keepY = c(3,3)

# ------------------------------- sPLS Framework ----------------------------- #

liver.spls <- spls(X, Y, 
                   ncomp = optimal.ncomp, 
                   keepX = optimal.keepX, 
                   keepY = optimal.keepY)

plotArrow(liver.spls) # DOESN'T produce error

# ------------------------------- rCCA Framework ----------------------------- #

liver.rcca <- rcc(X, Y, 
                  ncomp = optimal.ncomp, 
                  method = "shrinkage")

plotArrow(liver.rcca) # DOESN'T produce error
