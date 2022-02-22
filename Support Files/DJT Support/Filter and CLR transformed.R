#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prevalence filer - Microbiome data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Only Keep 72 samples with both omics data before applying filter
to_keep <- metab$Sample_ID
my_subset <- prune_samples((sample_names(microbiome) %in% to_keep), microbiome)
#new_physeq <- as(otu_table(my_subset), "matrix")
#new_physeq1 <- as.data.frame(new_physeq)

prevelancedf <- apply(X = otu_table(my_subset),
                      MARGIN = 1,
                      FUN = function(x){sum(x > 0)})

prevelancedf <- data.frame(Prevalence = prevelancedf,
                           TotalAbundance = taxa_sums(my_subset),
                           tax_table(my_subset)) # Add taxonomy and total read counts to this data.frame
prevelancedf[1:10,]

prevalenceThreshold <- 0.004 * nsamples(my_subset) #0.004 is present in at least 5 samples
prevalenceThreshold

keepTaxa <- rownames(prevelancedf)[(prevelancedf$Prevalence >= prevalenceThreshold)]
length(keepTaxa)

microbiome_0.004 <- prune_taxa(keepTaxa, my_subset)
microbiome_0.004

#confirm1 <- as(otu_table(microbiome_0.004), "matrix")
#confirm1 <- as.data.frame(confirm1)

saveRDS(microbiome_0.004, "C:/Users/trafi/Desktop/Paper 3/Files/Remote Access Files/ASV_besthit_phyloseq_0.004.rds")

# glom at the best taxonomic assignment (READ ABOUT THIS)
besthit_0.004 <- tax_glom(microbiome_0.004, taxrank="Genus") 
besthit_0.004

#confirm <- as(otu_table(besthit_0.004), "matrix")
#confirm <- as.data.frame(confirm)

saveRDS(besthit_0.004, "besthit_0.004.rds")
head(tax_table(besthit_0.004))
otu_table(besthit_0.004)[1:10,1:10]

OTU_1 <- as(otu_table(besthit_0.004), "matrix")
OTU_df <- as.data.frame(OTU_1)
OTU_dfT <- t(OTU_df)
OTU_dfT <- as.data.frame(OTU_dfT)

asv_fin1 <- data.frame(lapply(OTU_dfT, function(x) as.numeric(as.character(x))))

#Centered log ratio transformation (AFTER MATCHING IDs PROPERLY BELOW)
#1) offset
data_offset <- asv_fin1+1
sum(which(data_offset == 0)) # ok

#2) Log-ratio transformation (package from mixomics)
# we input the data as a matrix, here no need to add an offset as it is already done
data_clr <- logratio.transfo(as.matrix(data_offset), logratio = 'CLR', offset = 0)
class(data_clr) <- "matrix"
asv_clr <- as.data.frame(data_clr)
