# R

# Functional Genomics Center 2023
# 
#
# Using R with DNA & Proteomics Data
# @ jg@fgcz.ethz.ch
#


library(seqinr)
library(stringr)
library(readr)

# read in
myDNA <- read.fasta(file = "SA_6850_ncbi_sequence.fasta", seqtype = "DNA")
# look up pseudogenes in gff
myGFF <- read_tsv("SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff", skip = 7)
colnames(myGFF) <- c("chromosome", "source", "whatLevel", "start", "stop", "dot", "plusOrminus", "codingOrNot", "description")

# vcf
myVCF <- read_tsv("allS68Samples.vcf", skip = 28)
# for SA6850 we work in prX with B04_3 -> select this clone
cloneOfInterest <- "B04_3"

# find regions where we have SNVs
colID <- which(grepl(x = colnames(myVCF), pattern = cloneOfInterest))
SNVidx <- myVCF[,colID] != "./.:.:.:.:.:.:."
# get interesting stuff
SNVdf<- myVCF[SNVidx,1:7]


# on the myGFF
mGFFallGenes <- myGFF[myGFF$whatLevel == "gene", ] # only go for the gene tags and not the CDS
#gsub(x = sapply(strsplit(mGFFallGenes$description, split = ";"), function(x)x[6]), pattern = "locus_tag=", replacement = "") # ambigous becaues some have and additional gene= in
#gsub(x = sapply(strsplit(mGFFallGenes$description, split = ";"), function(x)x[5]), pattern = "locus_tag=", replacement = "") # see above
mGFFallGenes$geneName <- gsub(x = sapply(strsplit(mGFFallGenes$description, split = ";"), function(x)x[2]), pattern = "Name=", replacement = "")
mGFFallGenes$RSAU <- gsub(x = sapply(strsplit(mGFFallGenes$description, split = ";"), function(x)x[1]), pattern = "ID=gene-", replacement = "")
mGFFallGenes$locusTag <- gsub(x = sapply(strsplit(mGFFallGenes$description, split = ";"), function(x)x[1]), pattern = "ID=gene-", replacement = "")
mGFFallGenes$protAcc <- paste(mGFFallGenes$RSAU, sep = "_")
mGFFallGenes$RestOfFastaHeader <- paste("locusTag=",mGFFallGenes$RSAU, " ", mGFFallGenes$geneName, sep = "")

slimGFF <- mGFFallGenes[1:50,]

rm(fastaFileName)
# prepare fasta file
myLeadingSeq <- c("C","R","A", "P", "C","R","A", "P")
fastaFileName <- paste("SA6850_newFastaAA_correctedForAllChangesOf_v2_",cloneOfInterest,".fasta", sep = "")
write.fasta(sequences = myLeadingSeq, "MyCRAPproteinFirst this was generated to have a leading entry for the fasta file",file.out = fastaFileName, open = "a")

gffDF <- mGFFallGenes
tail(gffDF)

for (i in 1:nrow(gffDF)) {
  # get DNA stretch
  myORF <- myDNA$CP006706.1[(gffDF$start[i]-1):gffDF$stop[i]]
  
  # do we need to modify DNA stretch?
  isModified <- "notModified"
  offSetForPos <- 0
  for (j in 1:nrow(SNVdf)) {
    if (SNVdf$POS[j] > gffDF$start[i] && SNVdf$POS[j] < gffDF$stop[i]) {
      print(paste("my i is: ", i,"we should do something at j=", j))
      print(paste("length of ORF: ", length(myORF), "original: ", (gffDF$stop[i] - gffDF$start[i]), "offset: ", offSetForPos))
      # call the exact position split in 2 parts and put together again
      # seqPrefixLen <- SNVdf$POS[j] - gffDF$start[i]
      #firstStretch <- getFrag(object = myDNA, begin = gffDF$start[i], end = SNVdf$POS[j])
      #firstStretchMinus1 <- getFrag(object = myDNA, begin = gffDF$start[i], end = (SNVdf$POS[j] - 1)) # till Pos
      # get stretch from ORF
      snvRefLength <- nchar(SNVdf$REF[j])
      snvAltLength <- nchar(SNVdf$ALT[j])
      posInOrf <- SNVdf$POS[j] - gffDF$start[i] - snvRefLength
      firstStretchMinus1 <- myORF[1:(posInOrf + offSetForPos)][!is.na(myORF[1:(posInOrf + offSetForPos)])] # (object = myORF, begin = 0, end = (SNVdf$POS[j] - 1)) # till Pos
      lastStretchPlusOne <- myORF[(length(firstStretchMinus1) + snvAltLength):length(myORF)]
      indelSNPwhatever <- SNVdf$ALT[j]
      newOrfSequence <- s2c(tolower(paste(paste(firstStretchMinus1, collapse = ""), indelSNPwhatever, paste(lastStretchPlusOne, collapse = ""),sep="",collapse = " ")))
      myORF <- newOrfSequence
      isModified <- "isModified"
      offSetForPos <- snvAltLength - snvRefLength
    }
  }
  
  # translate
  # myprot <- getTrans(myPseudogene,frame =  j, sens = "R")
  if (gffDF$plusOrminus[i] == "+") myProtSeq <- getTrans(myORF, frame = 1, sens = "F")
  # if (gffDF$plusOrminus[i] == "-") myProtSeq <- getTrans(myORF, frame = 1, sens = "R")
  if (gffDF$plusOrminus[i] == "-") myProtSeq <- getTrans(myORF, frame = 3, sens = "R")
  
  # myprot <- getTrans(myPseudogene,frame =  j, sens = "R")
  # myName <- paste(myLeadingTag, "", i,"_","rev","_frm",j, " ", myRSAU, "; ", myDesc, sep="")
  # write.fasta(sequences = myprot, myName,file.out = fastaFileName, open = "a")
  myName <- paste(gffDF$protAcc[i], gffDF$RestOfFastaHeader[i], isModified, paste("location:",gffDF$start[i],"-",gffDF$stop[i]," strand:",gffDF$plusOrminus[i], sep=""),sep=" ")
  write.fasta(sequences = myProtSeq, myName,file.out = fastaFileName, open = "a")
}


#
# all wrong on the minus strand.. how do we need to translate
# protein sequence starts with METLNSINIP
# https://biocyc.org/gene?orgid=GCF_000462955&id=RSAU_RS00035#tab=FTRS
# getTrans(rev(myORF), frame = 1, sens = "F")
# getTrans(rev(myORF), frame = 2, sens = "F")
# getTrans(rev(myORF), frame = 3, sens = "F")
# 
# getTrans(rev(myORF), frame = 1, sens = "R")
# getTrans(rev(myORF), frame = 2, sens = "R")
# getTrans(rev(myORF), frame = 3, sens = "R")
# 
# 
# getTrans(myORF, frame = 1, sens = "F")
# getTrans(myORF, frame = 2, sens = "F")
# getTrans(myORF, frame = 3, sens = "F")
# 
# getTrans(myORF, frame = 1, sens = "R")
# getTrans(myORF, frame = 2, sens = "R")
# getTrans(myORF, frame = 3, sens = "R") # we have a winner!

i <- 7
myORF <- myDNA$CP006706.1[(gffDF$start[i]-1):(gffDF$stop[i]+0)]
getTrans(myORF, frame = 3, sens = "R", )

# to be more investigated.. cases where there are multiple snps in the same orf
# [1] "my i is:  2204 we should do something at j= 86"
# [1] "length of ORF:  418 original:  416 offset:  0"
# [1] "my i is:  2204 we should do something at j= 87"
# [1] "length of ORF:  419 original:  416 offset:  -1"
# [1] "my i is:  2204 we should do something at j= 88"
# [1] "length of ORF:  420 original:  416 offset:  -1"

i <- 2204
myORF <- myDNA$CP006706.1[(gffDF$start[i]-1):(gffDF$stop[i]+0)]

gffDF[2204,]
SNVdf[86:88,]


# 
# 
# gffDF <- mGFFallGenes
# for (i in 1:nrow(gffDF)) {
#   # get DNA stretch
#   myORF <- myDNA$CP006706.1[(gffDF$start[i]-1):gffDF$stop[i]]
#   
#   # do we need to modify DNA stretch?
#   isModified <- "notModified"
#   for (j in 1:nrow(SNVdf)) {
#     if (SNVdf$POS[j] > gffDF$start[i] && SNVdf$POS[j] < gffDF$stop[i]) {
#       print("we should do something")
#       # call the exact position split in 2 parts and put together again
#       # seqPrefixLen <- SNVdf$POS[j] - gffDF$start[i]
#       #firstStretch <- getFrag(object = myDNA, begin = gffDF$start[i], end = SNVdf$POS[j])
#       firstStretchMinus1 <- getFrag(object = myDNA, begin = gffDF$start[i], end = (SNVdf$POS[j] - 1)) # till Pos
#       lastStretchPlusOne <- getFrag(object = myDNA, begin = (SNVdf$POS[j]+1), end = gffDF$stop[i])
#       indelSNPwhatever <- SNVdf$ALT[j]
#       newOrfSequence <- s2c(tolower(paste(paste(firstStretchMinus1[[1]], collapse = ""), indelSNPwhatever, paste(lastStretchPlusOne[[1]], collapse = ""),sep="",collapse = " ")))
#       myORF <- newOrfSequence
#       isModified <- "isModified"
#     }
#   }
#   
#   # translate
#   # myprot <- getTrans(myPseudogene,frame =  j, sens = "R")
#   if (gffDF$plusOrminus[i] == "+") myProtSeq <- getTrans(myORF, frame = 1, sens = "F")
#   if (gffDF$plusOrminus[i] == "-") myProtSeq <- getTrans(myORF, frame = 1, sens = "R")
#   
#   # myprot <- getTrans(myPseudogene,frame =  j, sens = "R")
#   # myName <- paste(myLeadingTag, "", i,"_","rev","_frm",j, " ", myRSAU, "; ", myDesc, sep="")
#   # write.fasta(sequences = myprot, myName,file.out = fastaFileName, open = "a")
#   myName <- paste(gffDF$protAcc[i], gffDF$RestOfFastaHeader[i], isModified, paste("location:",gffDF$start[i],"-",gffDF$stop[i]," strand:",gffDF$plusOrminus[i], sep=""),sep=" ")
#   write.fasta(sequences = myProtSeq, myName,file.out = fastaFileName, open = "a")
# }
# 
# 
# 

