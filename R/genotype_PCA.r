library(vcfR)
library(adegenet)
library(dplyr)
library(tidyverse)
library(ggplot2)

#Load in a vcf
vcf <- read.vcfR("../course_RCIG/JE2_set_forIGV/2020-09-25_mod_SAJE2_filteredFromAncestor_noOutlier.1.vcf", verbose = FALSE)

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

s.label(pca1$li)

#Plot a PCA
plot(pca1$eig)

#Plot loadings
loadingplot(pca1$c1, axis=1, thres=.07, lab.jitter=1)
loadingplot(pca1$c1, axis=2, thres=.07, lab.jitter=1)

#Create a plot using ggplot, labelling points according to metadata
merge(pca1$li, metadata, by = "row.names") %>% 
  ggplot(aes(Axis3, Axis4, color = Treatment)) + geom_point()
