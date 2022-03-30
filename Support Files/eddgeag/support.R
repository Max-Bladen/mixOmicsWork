library(mixOmics)
data(breast.TCGA) # load in the data

# set a list of all the X dataframes
data = list(miRNA = breast.TCGA$data.train$mirna, 
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)

Y = breast.TCGA$data.train$subtype # set the response variable as the Y df

BPPARAM <- BiocParallel::SnowParam(workers = parallel::detectCores()-1)

design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

test.keepX <- list(mrna = c(10, 30), mirna = c(15, 25), protein = c(4, 8))


tune.comp12 <- tune.block.splsda(X = data, Y = Y,
                                 ncomp = 2,
                                 test.keepX = test.keepX,
                                 design = design,
                                 nrepeat = 2,
                                 BPPARAM = BPPARAM
)
tune.comp12$choice.keepX


# Now tuning a new component given previous tuned keepX
already.tested.X = tune.comp12$choice.keepX
tune.comp13 <- tune.comp12 <- tune.block.splsda(X = data, Y = Y,
                                                ncomp = 3,
                                                already.tested.X = tune.comp12$choice.keepX,
                                                test.keepX = test.keepX,
                                                design = design,
                                                nrepeat = 2,
                                                BPPARAM = BPPARAM
)

tune.comp34$choice.keepX
