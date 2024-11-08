# R

# Functional Genomics Center 2022
#
#
# Using R with DNA & Proteomics Data
# @ jg@fgcz.ethz.ch
#
# Using R along with genomics resources to come up with new fasta sequences
#
#
#

library(seqinr)
library(stringr)
library(readr)

myDNA <- read.fasta(file = "SA_6850_ncbi_sequence.fasta", seqtype = "DNA")

myStartPos <- 563123
myEndPos <- 567056
posDelPos <- 566275

myDNA$CP006706.1[(posDelPos-1):(posDelPos+2)]
(orfLength <- myEndPos - myStartPos)

length(myDNA$CP006706.1)
newChromosomeWithDeletion <- myDNA$CP006706.1[-posDelPos]
length(newChromosomeWithDeletion)

myDNA$CP006706.1 <- myDNA$CP006706.1[-posDelPos]



mysdrD <- myDNA$CP006706.1[myStartPos:myEndPos]

myP <- getTrans(mysdrD,frame =  0, sens = "F")
write.fasta(sequences = myP, "testName",file.out = "MyProteinInAAfasta_deletion.fasta", open = "a")


for(i in 1:3) {
  myPP <- getTrans(mysdrD,frame =  i, sens = "F")
  myName <- paste0("sdrDprotein_", "fw", "_frm_", i, sep="")
  print("number of stops: ")
  print(sum(str_count(paste(myPP,sep = ""), pattern = "\\*")))
  write.fasta(sequences = myPP, names = myName,file.out = "MyProteinInAAfasta_deletion.fasta", open = "a")
}


write.fasta(sequences = myPP, names = "myName",file.out = "MyProteinInAAfasta.fasta", open = "a")


# translate all regions that are "frame shifted"
#
# FGCZ-Analysis (jonas)
#- Project: p3404
#- Date: 2022-11-23
#- Task: translate in all frames the pseudogenes of SA6850
#- Duration:
#- What has been done:

# look up pseudogenes in gff
myGFF <- read_tsv("SA_6850_GCA_000462955.1_ASM46295v1_genomic_forR.gff", skip = 7)
colnames(myGFF) <- c("chromosome", "source", "whatLevel", "start", "stop", "dot", "plusOrminus", "codingOrNot", "description")

# not only do frameshift but also pseudogene and truncated
stringOfInterest <- "frameshift"
sum(str_count(string = myGFF$description, pattern = stringOfInterest) > 0, na.rm = TRUE)
frameShift_idx <- which(str_count(string = myGFF$description, pattern = stringOfInterest) > 0)
stringOfInterest <- "pseudogene"
sum(str_count(string = myGFF$description, pattern = stringOfInterest) > 0, na.rm = TRUE)
pseudo_idx <- which(str_count(string = myGFF$description, pattern = stringOfInterest) > 0)
stringOfInterest <- "truncated"
sum(str_count(string = myGFF$description, pattern = stringOfInterest) > 0, na.rm = TRUE)
truncated_idx <- which(str_count(string = myGFF$description, pattern = stringOfInterest) > 0)

#combine  all idx of interest
all_idx <- unique(c(frameShift_idx, pseudo_idx, truncated_idx))
length(all_idx)

# extract from gff interestin regions
myInterstingRegionsGFF <- myGFF[all_idx, ]
myInterstingRegionsGFF$description
# read in DNA chromosome
myDNA <- read.fasta(file = "SA_6850_ncbi_sequence.fasta", seqtype = "DNA")

# approach -> deletions not accounted! (to not mess up positions within chromosome)
#allDeletions <- c()
#newChromosomeWithDeletion <- myDNA$CP006706.1[-posDelPos]

myLeadingTag <- "pseudoOrframeshift"
myDNA$CP006706.1 <- myDNA$CP006706.1[-posDelPos]

# all in the same direction
sum((myInterstingRegionsGFF$stop - myInterstingRegionsGFF$start) < 0) == 0

myLeadingSeq <- c("C","R","A", "P", "C","R","A", "P")
fastaFileName <- "SA6850_allLociNotInAAfasta_Translatedin3Frames_nextApproach.fasta"
write.fasta(sequences = myLeadingSeq, "MyCRAPproteinFirst this was generated to have a leading entry for the fasta file",file.out = fastaFileName, open = "a")

for (i in 1:nrow(myInterstingRegionsGFF)) {
  myPseudogene <- myDNA$CP006706.1[myInterstingRegionsGFF$start[i]:myInterstingRegionsGFF$stop[i]]
  myRSAU <- gsub(gsub(x = myInterstingRegionsGFF$description[i], pattern = ";.*", replacement = ""), pattern = "ID=gene-", replacement = "")
  #parse description
  myDesc <- gsub(gsub(x = myInterstingRegionsGFF$description[i], pattern = ".*description=", replacement = ""), pattern = ";gbkey.*", replacement = "")
  myNote <- gsub(gsub(x = myInterstingRegionsGFF$description[i], pattern = ".*;Note=", replacement = ""), pattern = ";description=.*", replacement = "")
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


# read fasta back in and modify stop codons

myFasta <- read.fasta(file = "SA6850_allLociNotInAAfasta_Translatedin3Frames_nextApproach.fasta", seqtype = "AA")
length(myFasta)

myFasta[[5]][1:3] <- "X"
getSequence(myFasta[5])
getSequence(myFasta[5])[getSequence(myFasta[5]) == "*"]
stopBool <- getSequence(myFasta[[5]]) == "*"

myFasta[[5]][stopBool] <- "X"
getSequence(myFasta[5])

# start in Fasta w CRAPentry
myLeadingSeq <- c("C","R","A", "P", "C","R","A", "P")
fastaFileName <- "SA6850_allStopsReplacedByAllPossibleAAs.fasta"
write.fasta(sequences = myLeadingSeq, "MyCRAPproteinFirst this was generated to have a leading entry for the fasta file",file.out = fastaFileName, open = "a")

myAAs <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
for (i in 2:length(myFasta)) {
  for (j in 1:length(myAAs)) {
    mySeqHere <- myFasta[[i]]
    # prepare new names
    protAcc <- getName(mySeqHere)
    protDesc <- getAnnot(mySeqHere)
    protDesc <- gsub(x = protDesc, pattern = paste(">",protAcc, sep=""), replacement = "")
    # replace
    stopBool <- getSequence(mySeqHere) == "*"
    mySeqHere[stopBool] <- myAAs[j]
    # also put in desc where it was replaced
    # do.call()
    allStops <- do.call(paste, c(as.list(which(stopBool)), sep = "_"))
    protDesc <- paste("replaceStopsFor ", myAAs[j],";", allStops, ";" ,protDesc, sep="")
    protAcc <- paste(protAcc, "_",myAAs[j]," ",protDesc, sep="")
    write.fasta(sequences = getSequence(mySeqHere), protAcc,file.out = fastaFileName, open = "a")
  }

}



