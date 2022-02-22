library(mixOmics)

data(breast.TCGA) # load in the data

data = list(a = breast.TCGA$data.train$mirna, 
            b = breast.TCGA$data.train$mrna,
            c = breast.TCGA$data.train$protein,
            d = breast.TCGA$data.train$mirna, 
            e = breast.TCGA$data.train$mrna,
            f = breast.TCGA$data.train$protein,
            g = breast.TCGA$data.train$protein, 
            h = breast.TCGA$data.train$mrna,
            i = breast.TCGA$data.train$protein,
            j = breast.TCGA$data.train$protein,
            k = breast.TCGA$data.train$protein)

Y = breast.TCGA$data.train$subtype

design = matrix(0.1, ncol = length(data), nrow = length(data), # for square matrix filled with 0.1s
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

multi.block.plsda = block.splsda(X = data, Y = Y, ncomp = 5, design = design) # form basic DIABLO model
my.colors <- c("blue", "red", "green")
groups <- breast.TCGA$data.train$subtype

# code from user
plotIndiv(multi.block.plsda,
          ind.names = FALSE,
          ellipse = TRUE,
          group = groups,
          col.per.group = my.colors,
          blocks = c(4:10),
          pch = c(rep(16,8)),
          cex = 1,
          X.label =c(expression(paste("PLS-DA component 1"))),
          Y.label =c(expression(paste("PLS-DA component 2"))),
          size.xlabel = rel(1), 
          size.ylabel = rel(1),
          size.axis = rel(1),
          size.legend = rel(1),
          layout = c(4,2),
          legend.title = "", legend.box = NULL, legend.position = "right",
          
          
          subtitle = c("Biochemical indices in jejunal tissue",
                       ""mRNA targets in jejunal tissue"",
                       ""Morphology characteristics of jejunal tissue"",
                       ""Cell types in jejunal tissue"",
                       ""Free amino acids in jejunal digesta"",
                       ""Protein bound amino acids in jejunal digesta"",
                       ""Amino metabolites in jejunal digesta""
          ),
          style = "graphics",
          title = "")

