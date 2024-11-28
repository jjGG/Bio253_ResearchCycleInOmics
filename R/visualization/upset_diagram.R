# R

# Functional Genomics Center 2024
# 
#
# Displaying sets of peptides
# @ timon.leupp@uzh.ch
# @ luuk.vanheugten@uzh.ch
#

library(tidyverse)
library(UpSetR)

rm(list = ls())

# protein fasta files:
# refseq_fasta <- read.fasta('C:/Users/timon/OneDrive - Universität Zürich UZH/General - Pseudogeniuses/data/genomic data/ncbi_6850_dataset/ncbi_dataset/data/GCF_000462955.1/protein.faa')
# genbank_fasta <- read.fasta('C:/Users/timon/OneDrive - Universität Zürich UZH/General - Pseudogeniuses/data/genomic data/ncbi_6850_dataset/ncbi_dataset/data/GCA_000462955.1/protein.faa')


# import peptide reports
refseq <- read_csv("refseq.csv", skip = 8)
refseq[nrow(refseq), ]
# remove last "end of file" row
refseq <- refseq[-nrow(refseq), ]
paste("num peptides with REFSEQ:", length(refseq$`Peptide Sequence`))

genbank <- read_csv("genbank.csv", skip = 8)
genbank[nrow(genbank), ]
genbank <- genbank[-nrow(genbank), ]
paste("num peptides with GENBANK:", length(genbank$`Peptide Sequence`))

uniprot <- read_csv("uniprot.csv", skip = 8)
uniprot[nrow(uniprot), ]
uniprot <- uniprot[-nrow(uniprot), ]
paste("num peptides with UNIPROT:", length(uniprot$`Peptide Sequence`))

db2 <- read_csv("db2.csv", skip = 8)
db2[nrow(db2), ]
db2 <- db2[-nrow(db2), ]
paste("num peptides with db2:", length(db2$`Peptide Sequence`))

db5 <- read_csv("db5.csv", skip = 8)
db5[nrow(db5), ]
db5 <- db5[-nrow(db5), ]
paste("num peptides with db5:", length(db5$`Peptide Sequence`))

db2_pseudo_uncorrected <- read_csv("db2_pseudo_uncorrected.csv",
                                   skip = 8)
db2_pseudo_uncorrected[nrow(db2_pseudo_uncorrected), ]
db2_pseudo_uncorrected <- db2_pseudo_uncorrected[-nrow(db2_pseudo_uncorrected), ]
paste("num peptides with db2_pseudo_uncorr:", length(db2_pseudo_uncorrected$`Peptide Sequence`))

uniprot_pseudo_corrected <- read_csv("uniprot_pseudo_corrected.csv", 
                                     skip = 8)
uniprot_pseudo_corrected[nrow(uniprot_pseudo_corrected), ]
uniprot_pseudo_corrected <- uniprot_pseudo_corrected[-nrow(uniprot_pseudo_corrected), ]
paste("num peptides with uniprot_pseudo_corr:", length(uniprot_pseudo_corrected$`Peptide Sequence`))


listInput <- list(db2 = db2$`Peptide Sequence`, pseudo_uncorrected = db2_pseudo_uncorrected$`Peptide Sequence`,
                  db5 = db5$`Peptide Sequence`, refseq = refseq$`Peptide Sequence`,
                  uniprot = uniprot$`Peptide Sequence`, pseudo_corrected = uniprot_pseudo_corrected$`Peptide Sequence`
                  )

# upset venn diagram
upset(fromList(listInput), nsets = 7, nintersects = 25, order.by = "freq", 
      set_size.show = TRUE, set_size.scale_max = 30000)




