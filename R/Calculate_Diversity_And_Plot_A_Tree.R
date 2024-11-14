library(vcfR)
library(ape)
library(adegenet)
library(poppr)
library(reshape2)
library(dplyr)
library(tidyverse)

#Load vcf
vcf <- read.vcfR("../../course_Research cycle in genomics/JE2_set_forIGV/2020-09-25_mod_SAJE2_filteredFromAncestor_noOutlier.1.vcf", verbose = FALSE)

#Convert to Genind
my_genind <- vcfR2genind(vcf, check.ploidy=TRUE)
ploidy(my_genind) <-  rep(1,50) # Correct ploidy
remove_indices <- c(27,38,44,50) # Remove "P"s
my_genind <- my_genind[-remove_indices, ]

#Create metadata
metadata <- data.frame(
  ID = indNames(my_genind),  # Make sure IDs match those in your genind object
  Population = sapply(str_split(indNames(my_genind), "_"), .subset, 2))

#Add pop to gening object
my_genind@pop <- as.factor(metadata$Population)
xT <- tab(my_genind, freq=TRUE, NA.method="zero")

#Calculate diversity
df = diversity_ci(my_genind, n = 100L, raw = FALSE)

#Plot tree
plot.phylo(upgma(dist(xT)))
