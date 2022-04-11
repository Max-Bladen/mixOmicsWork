
setwd("C:/Users/Work/Desktop/mixOmics/R")

devtools::load_all()

# ---------------------------------------------------------------------------- #
# ================================== (s)PCA ================================== #
# ---------------------------------------------------------------------------- #

data(multidrug)
X <- multidrug$ABC.trans

# ---------------------------------------------------------------------------- #
# ================================ (s)PLS-DA ================================= #
# ---------------------------------------------------------------------------- #

data(srbct)
X <- srbct$gene
Y <- srbct$class

# ---------------------------------------------------------------------------- #
# ================================== (s)PLS ================================== #
# ---------------------------------------------------------------------------- #

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

# ---------------------------------------------------------------------------- #
# ================================== (r)CCA ================================== #
# ---------------------------------------------------------------------------- #

data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene

# ---------------------------------------------------------------------------- #
# =============================== Block.(s)PLS =============================== #
# ---------------------------------------------------------------------------- #

data(breast.TCGA)
X = list(miRNA = breast.TCGA$data.train$mirna,
         mRNA = breast.TCGA$data.train$mrna,
         proteomics = breast.TCGA$data.train$protein)

# ---------------------------------------------------------------------------- #

data.GH.URL <- "https://github.com/Max-Bladen/mixOmicsWork/blob/main/Development/Data/nmt_data_processed.RData?raw=true"
.load_url(data.GH.URL)
X1 <- data$rna
X2 <- data$met_genebody
X3 <- data$acc_genebody
X <- list(rna = X1, methylation = X2, accessibility = X3) 

# ---------------------------------------------------------------------------- #
# ================================== DIABLO ================================== #
# ---------------------------------------------------------------------------- #

data(breast.TCGA)
X = list(miRNA = breast.TCGA$data.train$mirna,
         mRNA = breast.TCGA$data.train$mrna,
         proteomics = breast.TCGA$data.train$protein)
Y = breast.TCGA$data.train$subtype

# ---------------------------------------------------------------------------- #

data(breast.TCGA)
X = list(miRNA = breast.TCGA$data.test$mirna,
         mRNA = breast.TCGA$data.test$mrna)
Y = breast.TCGA$data.test$subtype

# ---------------------------------------------------------------------------- #
# ============================== MINT.(s)PLSDA =============================== #
# ---------------------------------------------------------------------------- #

data(stemcells)
X <- stemcells$gene
Y <- stemcells$celltype
s <- stemcells$study

# ---------------------------------------------------------------------------- #
# =========================== MINT.Block.(s)PLSDA ============================ #
# ---------------------------------------------------------------------------- #

data(breast.TCGA)
mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
X = list(mrna = mrna, mirna = mirna)
Y = c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)

study = c(rep("study1",150), rep("study2",70))


# ---------------------------------------------------------------------------- #
# ================================ Multilevel ================================ #
# ---------------------------------------------------------------------------- #

data(vac18)
X <- vac18$genes
Y <- vac18$stimulation
ml <- vac18$sample