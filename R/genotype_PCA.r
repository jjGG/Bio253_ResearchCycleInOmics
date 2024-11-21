library(vcfR)
library(adegenet)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)

#Load in a vcf
vcf <- read.vcfR("../../../course_Research cycle in genomics/JE2_set_forIGV/2020-09-25_mod_SAJE2_filteredFromAncestor_noOutlier.1.vcf", verbose = FALSE)

#Convert a vcf to genind object and then to df
x <- vcfR2genind(vcf, return.alleles = TRUE, ploidy = 1, sep = "/")
xT <- tab(x, freq=TRUE, NA.method="zero")

#Create a metadata file
metadata = as.data.frame(rownames(xT)) %>% 
  mutate(ID = rownames(xT)) %>% 
  separate('rownames(xT)', into = c("clone", "Treatment", "Tube", "Clone")) %>% 
  column_to_rownames("ID")

#Filter "P" data from the sequence
xT = xT[rownames(metadata[!metadata$Clone == "P",]),]
metadata = metadata[!metadata$Clone == "P",]

#Run PCA
pca1 <- dudi.pca(df = xT, center = TRUE, scale = FALSE, scannf = FALSE, nf = 15)

#Create PCA with labels
s.label(pca1$li)

#Plot a PCA
plot(pca1$eig)

#Plot loadings
loadingplot(pca1$c1, axis=1, thres=.07, lab.jitter=1)
loadingplot(pca1$c1, axis=2, thres=.07, lab.jitter=1)

#Phenotype
phenotype = readxl::read_xlsx("../../../course_Research cycle in genomics/resources/SA_evolution_phenotypes_v2_noCowan.xlsx")
phenotype = phenotype[phenotype$SA_strain_ID == "JE2",]
phenotype = phenotype %>% 
  mutate(sample = if_else(evolution_condition == "supernatant", paste0(SA_strain_ID,"_", "S"), 
                          paste0(SA_strain_ID,"_", "T"))) 
phenotype$sample_ID = if_else(nchar(sapply(str_split(phenotype$strain_ID, "\\."), .subset,1)) == 2, 
                              paste0(sapply(str_split(sapply(str_split(phenotype$strain_ID, "\\."), .subset,1), ""), .subset, 1), "0", sapply(str_split(sapply(str_split(phenotype$strain_ID, "\\."), .subset,1), ""), .subset, 2)), sapply(str_split(phenotype$strain_ID, "\\."), .subset,1))
phenotype$sample_ID_2 = sapply(str_split(phenotype$strain_ID, "\\."), .subset,2)
phenotype$sample =  paste0(phenotype$sample, "_", phenotype$sample_ID, "_", phenotype$sample_ID_2)
phenotype = phenotype %>% column_to_rownames("sample")
temp = merge(pca1$li, phenotype, by = "row.names")
melt(temp, id.vars = colnames(temp)[c(1,17:33)]) %>% ggplot(aes(value, average_H2O2_survival, color = variable)) + geom_line()

##Calculate the correlation of a PC with a particular phenotype
correlations <- sapply(temp[, 2:16], function(pc) cor(as.numeric(pc), temp$average_H2O2_survival, use = "complete.obs"))

# Sort and view top correlations
sort(abs(correlations), decreasing = TRUE)[1:10]

#Create a plot using ggplot, labelling points according to metadata
merge(pca1$li, metadata, by = "row.names") %>% 
  ggplot(aes(Axis3, Axis4, color = Treatment)) + geom_point()


merge(pca1$li, phenotype, by = "row.names") %>% data.frame() %>%
  ggplot(aes(Axis3, Axis4, color = average_H2O2_survival)) + geom_point()
