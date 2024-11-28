# R

# Functional Genomics Center 2024
# 
#
# Using R with DNA & Proteomics Data
# @ timon.leupp@uzh.ch
# @ luuk.vanheugten@uzh.ch
#

library(VariantAnnotation)
library(Biostrings)


rm(list = ls())

apply_vcf_to_fasta <- function(vcf_file, fasta_file, output_fasta, variant_id) {
  # Read the genome from the FASTA file (assumes a single chromosome)
  fasta <- readDNAStringSet(fasta_file)
  if (length(fasta) > 1) {
    stop("FASTA contains multiple sequences; this script assumes a single chromosome.")
  }
  genome_seq <- as.character(fasta[[1]])
  
  # Read the VCF file
  vcf <- readVcf(vcf_file, genome = names(fasta))
  # Check if the variant ID exists
  if (!variant_id %in% colnames(geno(vcf)$GT)) {
    stop(paste("Variant ID", variant_id, "not found in VCF."))
  }
  
  # Extract genotypes for the given variant ID
  gt <- geno(vcf)$GT[, variant_id]
  
  # Process each variant, sorted by position to handle sequential adjustments
  chroms <- rowRanges(vcf)
  variants <- data.frame(
    pos = start(chroms),
    ref = as.character(ref(vcf)),
    alt = as.character(unlist(alt(vcf))),
    gt = gt,
    stringsAsFactors = FALSE
  )
  variants <- variants[order(variants$pos), ]
  
  # Initialize position shift
  shift <- 0
  
  for (i in seq_along(gt)) {
    if (gt[i] == "0/1" || gt[i] == "1/1") {
      ref <- as.character(ref(vcf)[i])
      alt <- as.character(unlist(alt(vcf))[i])
      # Adjusted position accounting for previous shifts
      pos <- start(chroms[i]) + shift  # Adjusted position for previous shifts
      
      # Modify the sequence
      if (nchar(ref) == 1 && nchar(alt) == 1) {
        # SNP: Replace a single nucleotide
        substr(genome_seq, pos, pos) <- alt
      } else {
        # Indels: Update the sequence after insertion/deletion 
        # and calculate the positional shift
        genome_seq <- paste0(
          substr(genome_seq, 1, pos - 1),
          alt,
          substr(genome_seq, pos + nchar(ref), nchar(genome_seq))
        )
        # update shift
        shift <- shift + (nchar(alt) - nchar(ref))
      }
    }
  }
  
  updated_fasta <- DNAStringSet(genome_seq)
  names(updated_fasta) <- names(fasta)
  writeXStringSet(updated_fasta, output_fasta)
}

vcf_file <- "allS68Samples.vcf"
fasta_file <- "SA_6850_ncbi_sequence.fasta"

output_fasta <- "output.fasta"
variant_id <- "S68_ancestr"

apply_vcf_to_fasta(vcf_file, fasta_file, output_fasta, variant_id)
