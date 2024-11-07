#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# This script attempts to parse the description lines of a fasta file into a data frame
# It is a bit more complex than the one in the bash script because it has to deal with the fact that the description lines are not uniform
# The script is written in R and uses the tidyverse packages as well as stringr


# read in all description lines from txt file and parse into frame
# get rid of the first 2
# take out MS contaminants
# treat pseudoOrframeshift proteins differently

library(readr)
library(dplyr)
library(stringr)

# this txt file with all description lines is generated with grep on a bash console (cat *myfastaofInterest.fasta | grep '>' > myFile.txt)
# the first 2 proteins are describing the database generation and are CRAPCRAP
txt <- read_tsv("../data/SA6850_db5_allDescriptionLines.txt", col_names = FALSE, skip = 2)

# contaminant bool
# mass spec related -> get rid of Contaminants
bool_Nocontaminants <- grepl(x = txt$X1, pattern = "FGCZCont") == FALSE
table(bool_Nocontaminants)

MyAcc <- data.frame(txt[which(bool_Nocontaminants),])
nrow(MyAcc)
colnames(MyAcc) <- c("originalFullLine")

# init df
df <- data.frame(matrix(nrow = nrow(MyAcc), ncol = 16))
colnames(df) <- c("originalFullLine", "Acc", "protType","gbkey", "locusTag", "proteinDesc", "SAprotID", "SAlocation", "SAGenedirection", "SALocationStart", "SALocationEnd", "OrthoFull", "OrthoName", "OrthoAcc", "OrthoUPname", "OrthoEvalue")

# Acc
df$originalFullLine <- MyAcc$originalFullLine
df$Acc <- gsub(x = sapply(strsplit(df$originalFullLine, split = " "), function(x)x[1]), pattern = ">", replacement = "")


# pseudo or not
df$protType <- "CDS"
df$protType[grep(x = df$Acc, pattern = "pseudoOrframeshift")] <- "pseudoOrframeshift"
table(df$protType)

# split here for further processing
df_pseudo <- df[df$protType == "pseudoOrframeshift", ]
df_cds <- df[df$protType == "CDS", ]

# work ncbi proteins and orthologue
df_cds$gbkey <- gsub(x = df_cds$originalFullLine, pattern = ".*\\[gbkey=(.*)\\].*", replacement = "\\1")
#df_cds$gbkey[grep(x = df$Acc, pattern = "pseudoOrframeshift")] <- "pseudoOrframeshift"

# regexxing for the rest
# add more meta info columns to output
df_cds$locusTag <- gsub(x = df_cds$originalFullLine, pattern = ".*\\[locus_tag=(.*)\\] \\[protein=.*", replacement = "\\1")
df_cds$proteinDesc <- gsub(x = df_cds$originalFullLine, pattern = ".*\\[protein=(.*)\\] \\[protein_id=.*", replacement = "\\1")
df_cds$SAprotID <- gsub(x = df_cds$originalFullLine, pattern = ".*\\[protein_id=(.*)\\] \\[location=.*", replacement = "\\1")
df_cds$SAlocation <- gsub(x = df_cds$originalFullLine, pattern = ".*\\[location=(.*)\\] \\[gbkey=.*", replacement = "\\1")
df_cds$SAGenedirection <- "NA"
df_cds$SAGenedirection[grep(x = df_cds$SAlocation, pattern = "complement")] <- "revStrand"
df_cds$SAGenedirection[!is.na(as.numeric(sapply(strsplit(df_cds$SAlocation, split = "\\.\\."), function(x)x[1])))] <- "fwStrand"

# normal genes
# fw strand
df_cds$SALocationStart[df_cds$SAGenedirection == "fwStrand"] <- as.numeric(sapply(strsplit(df_cds$SAlocation[df_cds$SAGenedirection == "fwStrand"], split = "\\.\\."), function(x)x[1]))
df_cds$SALocationEnd[df_cds$SAGenedirection == "fwStrand"] <- as.numeric(sapply(strsplit(df_cds$SAlocation[df_cds$SAGenedirection == "fwStrand"], split = "\\.\\."), function(x)x[2]))

# rev strand
df_cds$SALocationEnd[df_cds$SAGenedirection == "revStrand"] <- as.numeric(gsub(x = sapply(strsplit(df_cds$SAlocation[df_cds$SAGenedirection == "revStrand"], split = "\\.\\."), function(x)x[1]), pattern = "complement\\(",replacement=""))
df_cds$SALocationStart[df_cds$SAGenedirection == "revStrand"] <- as.numeric(gsub(x = sapply(strsplit(df_cds$SAlocation[df_cds$SAGenedirection == "revStrand"], split = "\\.\\."), function(x)x[2]), pattern = "\\)",replacement=""))


# blast related
# some do not have an orthologue .. we again split it and only do it for those that have it
orthoFound_bool <- grepl(x = df_cds$originalFullLine, pattern = "BLASTORTHO")
table(orthoFound_bool)
df_cds_wOrtho <- df_cds[orthoFound_bool,]
df_cds_noOrtho <- df_cds[orthoFound_bool==FALSE,]


df_cds_wOrtho$OrthoFull<- gsub(x = df_cds_wOrtho$originalFullLine, pattern = ".*BLASTORTHO (.*)", replacement = "\\1")
df_cds_wOrtho$OrthoName <- gsub(x = df_cds_wOrtho$OrthoFull, pattern = "(.*) evalue=.*", replacement = "\\1")
df_cds_wOrtho$OrthoAcc <- sapply(strsplit(df_cds_wOrtho$OrthoName, split = "_"), function(x)x[2])
df_cds_wOrtho$OrthoUPname <- sapply(strsplit(df_cds_wOrtho$OrthoName, split = "_"), function(x)x[3])
df_cds_wOrtho$OrthoEvalue <- as.numeric(gsub(x = df_cds_wOrtho$OrthoFull, pattern = ".* evalue=(.*)", replacement = "\\1"))

ncol(df_cds_wOrtho)
colnames(df_cds_wOrtho)

# work on pseudo
# locusTag <- split by ; is No 1 (get rid of >pseudoOr... )
# description <- split by ; is the last
df_pseudo$locusTag <- gsub(x = sapply(strsplit(df_pseudo$originalFullLine, split = "; "), function(x)x[1]), pattern = ">.* ", replacement = "")

# # pseudo genes
df_pseudo$SAGenedirection[grep(x = df_pseudo$Acc, pattern = "pseudo.*_fw_")] <- "fwStrand"
df_pseudo$SAGenedirection[grep(x = df_pseudo$Acc, pattern = "pseudo.*_rev_")] <- "revStrand"


# get the correct description out
listSplitted <- strsplit(df_pseudo$originalFullLine, split = "; ")

for (i in 1:length(listSplitted)) {
  df_pseudo$proteinDesc[i] <- sapply(strsplit(df_pseudo$originalFullLine[i], split = "; "), function(x)x[length(listSplitted[[i]])])
}


# combine all of them again
resDF <- rbind(df_pseudo, df_cds_wOrtho, df_cds_noOrtho)

# write out
write_tsv(x = resDF, file = "../data/p3404_db5_SA6850_fastaHeadersParsedIntoDF.tsv")

protLen <- abs(resDF$SALocationEnd - resDF$SALocationStart)
hist(protLen)
