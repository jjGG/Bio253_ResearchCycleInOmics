# R

# Functional Genomics Center 2024
#
#
# Using R with DNA & Proteomics Data
# written by
# @ jg@fgcz.ethz.ch
# modified and adapted by
# @ timon.leupp@uzh.ch
# @ luuk.vanheugten@uzh.ch
#
# Using R along with genomics resources to come up with new protein fasta sequences
#

library(seqinr)
library(tidyverse)
library(rtracklayer)

rm(list = ls())


# translate all regions that are annotated as type "pseudogene"
#
# FGCZ-Analysis
# adapted and modified by Timon Leupp and Luuk van Heugten in November 2024
#- Project: Block course research cycle in (gen)omics
#- Task: translate in all frames the pseudogenes of SA6850

# read in DNA sequence
myDNA <- read.fasta(file = "SA_6850_ncbi_sequence_corrected_for_ancestor_mutation.fasta", seqtype = "DNA")

myLeadingTag <- "pseudo_TLLVH_"

myLeadingSeq <- c("C","R","A", "P", "C","R","A", "P")
fastaFileName <- "SA6850_pseudogenes_translated_in_3_frames_TLLVH.fasta"
write.fasta(sequences = myLeadingSeq, 
            "MyCRAPproteinFirst this was generated to have a leading entry for the fasta file",
            file.out = fastaFileName, 
            open = "a")

## look up pseudogenes in gff

# approach for old 2013 Genbank GFF
oldGenbankGFF <- read_tsv("SA_6850_GCA_000462955.1_ASM46295v1_genomic_forR.gff", 
                          skip = 7)
colnames(oldGenbankGFF) <- c("chromosome", "source", "whatLevel", "start", 
                             "stop", "dot", "plusOrminus", "codingOrNot", "description")

# not only do pseudogenes but also frameshifts and truncated
stringOfInterest <- "pseudogene"
sum(str_count(string = oldGenbankGFF$description, pattern = stringOfInterest) > 0, na.rm = TRUE)
pseudo_idx <- which(str_count(string = oldGenbankGFF$description, pattern = stringOfInterest) > 0)
stringOfInterest <- "frameshift"
sum(str_count(string = oldGenbankGFF$description, pattern = stringOfInterest) > 0, na.rm = TRUE)
frameShift_idx <- which(str_count(string = oldGenbankGFF$description, pattern = stringOfInterest) > 0)
stringOfInterest <- "truncated"
sum(str_count(string = oldGenbankGFF$description, pattern = stringOfInterest) > 0, na.rm = TRUE)
truncated_idx <- which(str_count(string = oldGenbankGFF$description, pattern = stringOfInterest) > 0)

#combine  all idx of interest
all_idx <- unique(c(frameShift_idx, pseudo_idx, truncated_idx))
length(all_idx)

# extract from gff interesting regions
myInterstingRegionsGFF <- myGFF[all_idx, ]

# all in the same direction
sum((myInterstingRegionsGFF$end - myInterstingRegionsGFF$start) < 0) == 0

#approach for 2013 Genbank
for (i in 1:nrow(myInterstingRegionsGFF)) {
  myPseudogene <- myDNA$CP006706.1[myInterstingRegionsGFF$start[i]:myInterstingRegionsGFF$stop[i]]
  myRSAU <- gsub(
    gsub(x = myInterstingRegionsGFF$description[i], pattern = ";.*", replacement = ""), 
    pattern = "ID=gene-", replacement = "")
  #parse description
  myDesc <- gsub(
    gsub(x = myInterstingRegionsGFF$description[i], pattern = ".*description=", replacement = ""), 
    pattern = ";gbkey.*", replacement = "")
  myNote <- gsub(
    gsub(x = myInterstingRegionsGFF$description[i], pattern = ".*;Note=", replacement = ""), 
    pattern = ";description=.*", replacement = "")
  # go in frames here fw
  for (j in 1:3) {
    myprot <- getTrans(myPseudogene,frame =  j, sens = "F")
    myName <- paste(myLeadingTag, "", i,"_","fw","_frm",j, " ", myRSAU, "; ", myNote, "; ", myDesc, sep="")
    write.fasta(sequences = myprot, myName,file.out = fastaFileName, open = "a")
  }
  # go in frames here rev
  for (j in 1:3) {
    myprot <- getTrans(myPseudogene,frame =  j, sens = "R")
    myName <- paste(myLeadingTag, "", i,"_","rev","_frm",j, " ", myRSAU, "; ", myDesc, sep="")
    write.fasta(sequences = myprot, myName,file.out = fastaFileName, open = "a")
  }
}


# approach for new 2023 RefSeq annotated GFF
newGFF <- as.data.frame(import.gff("SA_6850_2023_REFSEQ_annotation.gff"))

idx_of_interest <- which(str_count(string = newGFF$type, pattern = 'pseudogene') > 0)
InterstingRegionsNewGFF <- newGFF[idx_of_interest, ]

# modified for 2023 RefSeq GFF
file_conn <- file("SA_6850_REFSEQ_pseudogenes_protein_description_lines.txt", open = "a")
for (i in 1:nrow(InterstingRegionsNewGFF)) {
  myPseudogene <- myDNA$NC_022222.1[InterstingRegionsNewGFF$start[i]:InterstingRegionsNewGFF$end[i]]
  myRSAU <- gsub(InterstingRegionsNewGFF$ID[i], pattern = "gene-", replacement = "")
  myDesc <- InterstingRegionsNewGFF$Name[i]
  myNote <- InterstingRegionsNewGFF$ID[i]
  myLocation <- paste0("[location=", InterstingRegionsNewGFF$start[i], "..", InterstingRegionsNewGFF$end[i], "]")
  
  # go in frames here forward
  for (j in 1:3) {
    myprot <- getTrans(myPseudogene, frame =  j, sens = "F")
    myName <- paste(myLeadingTag, "", i,"_","fw","_frm",j, " ", myRSAU, "; ", myLocation, "; ", myDesc, "; ", myNote, sep="")
    cat(paste0(">", myName), "\n", file = file_conn)
    write.fasta(sequences = myprot, myName,file.out = fastaFileName, open = "a")
  }
  # go in frames here reverse
  for (j in 1:3) {
    myprot <- getTrans(myPseudogene,frame =  j, sens = "R")
    myName <- paste(myLeadingTag, "", i,"_","rev","_frm",j, " ", myRSAU, "; ", myLocation, "; ", myDesc, "; ", myNote, sep="")
    cat(paste0(">", myName), myName, "\n", file = file_conn)
    write.fasta(sequences = myprot, myName,file.out = fastaFileName, open = "a")
  }
}
close(file_conn)

