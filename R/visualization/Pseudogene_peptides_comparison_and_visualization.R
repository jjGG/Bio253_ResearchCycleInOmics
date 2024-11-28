# Research cycle in omics script for comparing the fictional (pseudogene) peptides
# with the proteomics data
# By Luuk van Heugten and Timon Leupp
# both authors contributed equally to this script and deserve equal credit

# reset global environment

rm(list=ls())

#load packages
library(tidyverse)
library(skimr)
library(seqinr)
library(rtracklayer)
library(ggVennDiagram)


# import peptide report tables 
# skip the metadata at the beginning

db2 <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_db2.csv", 
                skip = 8)

db2_wPseudogenes <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_db2_wPseudogenes.csv",
                             skip = 8)

db5 <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_db5_20241111_v2.csv",
                skip = 8)

Genbank <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_Genbank.csv",
                    skip = 8)

UPreference <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_UPreference.csv",
                        skip = 8)

wREFSEQ <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_wREFSEQ.csv",
                    skip = 8)

newORF_wPseudo_v2 <- read_csv("Peptide Report of Spectronaut2SDIA_SA6850_newORF_wPseudo_v2.csv", 
                              skip = 8)

# This script does:
# compare peptide sequences across the different data sets.
# list the sequences that match across different datasets in = equal
# list the sequences that dont match across = special

# used for debugging
peptide_data1 <- db2_wPseudogenes
peptide_data2 <- db5

# function to create a venn diagram for two datasets
# also saves matched information to csv file
# usage: create_vd_peptides(
#                           dataframe of peptide report 1, 
#                           dataframe of peptide report 2,
#                           TRUE if peptides map to pseudogenes should be visible for dataset 1,
#                           TRUE if peptides map to pseudogenes should be visible for dataset 2)
create_vd_peptides <- function(peptide_data1, peptide_data2, show_pseudo_1, show_pseudo_2){
  
  # extract the variable names of the input arguments of the function
  dp1 <- deparse(substitute(peptide_data1))
  dp2 <- deparse(substitute(peptide_data2))
  
  # find matches for both datasets
  matches <- intersect(peptide_data1$`Peptide Sequence`, peptide_data2$`Peptide Sequence`)
  
  # extract protein names
  matched_data1 <- peptide_data1[peptide_data1$`Peptide Sequence` %in% matches, ]
  matched_data2 <- peptide_data2[peptide_data2$`Peptide Sequence` %in% matches, ]
  
  matched_df <- data.frame(Peptide = matches, 
                           Status = "Matched",
                           Protein_Name1 = matched_data1$`Protein Name`,
                           Protein_Name2 = matched_data2$`Protein Name`,
                           Protein_Accession1 = matched_data1$`Protein Accessions`,
                           Protein_Accession2 = matched_data2$`Protein Accessions`)
  # this means all the peptide sequences here are in both datasets
  
  # find peptides in dataset1 but not in dataset2 and vice versa
  only_data1 <- setdiff(peptide_data1$`Peptide Sequence`, peptide_data2$`Peptide Sequence`)
  
  only_data2 <- setdiff(peptide_data2$`Peptide Sequence`, peptide_data1$`Peptide Sequence`)
  
  # extract protein names
  unique_data1 <- peptide_data1[peptide_data1$`Peptide Sequence` %in% only_data1, ]
  unique_data2 <- peptide_data2[peptide_data2$`Peptide Sequence` %in% only_data2, ]
  
  not_matched_df <- data.frame(Peptide = c(unique_data1$`Peptide Sequence`, unique_data2$`Peptide Sequence`),
                               Status = c(rep( dp1, nrow(unique_data1)),
                                          rep( dp2, nrow(unique_data2))),
                               Protein_Name = c(unique_data1$`Protein Name`, unique_data2$`Protein Name`),
                               Protein_Accessions = c(unique_data1$`Protein Accessions`, unique_data2$`Protein Accessions`))
 
  # Rename columns for Protein Names dynamically
  colnames(matched_df)[colnames(matched_df) == "Protein_Name1"] <- paste0("Protein_Name_", dp1)
  colnames(matched_df)[colnames(matched_df) == "Protein_Name2"] <- paste0("Protein_Name_", dp2)
  colnames(matched_df)[colnames(matched_df) == "Protein_Accession1"] <- paste0("Protein_Accessions_", dp1)
  colnames(matched_df)[colnames(matched_df) == "Protein_Accession2"] <- paste0("Protein_Accessions_", dp2)
  
  
  # filter pseudogene-associated entries from the Protein Accessions column
  pseudogene_data1 <- peptide_data1[grepl("pseudo_", peptide_data1$`Protein Accessions`, ignore.case = TRUE), ]
  pseudogene_data2 <- peptide_data2[grepl("pseudo_", peptide_data2$`Protein Accessions`, ignore.case = TRUE), ]
  
  # Calculate pseudogene-specific overlaps and uniques
  pseudo_matches <- intersect(pseudogene_data1$`Peptide Sequence`, pseudogene_data2$`Peptide Sequence`)
  pseudo_only_data1 <- setdiff(pseudogene_data1$`Peptide Sequence`, pseudogene_data2$`Peptide Sequence`)
  pseudo_only_data2 <- setdiff(pseudogene_data2$`Peptide Sequence`, pseudogene_data1$`Peptide Sequence`)
  pseudo_data1_in_data2 <- intersect(pseudo_only_data1, peptide_data2$`Peptide Sequence`)
  pseudo_data2_in_data1 <- intersect(pseudo_only_data2, peptide_data1$`Peptide Sequence`)
  pseudo_exclusive_data_1 <- setdiff(pseudo_only_data1, peptide_data2$`Peptide Sequence`)
  pseudo_exclusive_data_2 <- setdiff(pseudo_only_data2, peptide_data1$`Peptide Sequence`)
  
  # Create a separate pseudogene DataFrame for matches
  # it can be that no pseudogenes were even annotated in the dataset
  # therefore if statement with in case an empty data frame
  # Create a separate pseudogene DataFrame for matches
  # it can be that no pseudogenes were even annotated in the dataset
  # therefore if statement with in case an empty data frame
  if (length(pseudo_matches) > 0) {
    pseudogene_matched_df <- data.frame(
      Peptide = pseudo_matches,
      Status = "Matched Pseudogenes",
      Protein_Name1 = pseudogene_data1$`Protein Name`[pseudogene_data1$`Peptide Sequence` %in% pseudo_matches],
      Protein_Name2 = pseudogene_data2$`Protein Name`[pseudogene_data2$`Peptide Sequence` %in% pseudo_matches],
      Protein_Accession1 = pseudogene_data1$`Protein Accessions`[pseudogene_data1$`Peptide Sequence` %in% pseudo_matches],
      Protein_Accession2 = pseudogene_data2$`Protein Accessions`[pseudogene_data2$`Peptide Sequence` %in% pseudo_matches]
    )
  } else {
    # Create an empty data frame if no matches are found
    pseudogene_matched_df <- data.frame(
      Peptide = character(0),
      Status = character(0),
      Protein_Name1 = character(0),
      Protein_Name2 = character(0),
      Protein_Accession1 = character(0),
      Protein_Accession2 = character(0)
    )
  }
  
  if (length(c(pseudo_only_data1, pseudo_only_data2)) > 0) {
    not_matched_pseudogenes_df <- data.frame(
      Peptide = c(pseudo_only_data1, pseudo_only_data2),
      Status = c(
        rep(paste0("Pseudogenes only in ", dp1), length(pseudo_only_data1)),
        rep(paste0("Pseudogenes only in ", dp2), length(pseudo_only_data2))
      ),
      Protein_Name = c(
        pseudogene_data1$`Protein Name`[pseudogene_data1$`Peptide Sequence` %in% pseudo_only_data1],
        pseudogene_data2$`Protein Name`[pseudogene_data2$`Peptide Sequence` %in% pseudo_only_data2]
      ),
      Protein_Accessions = c(
        pseudogene_data1$`Protein Accessions`[pseudogene_data1$`Peptide Sequence` %in% pseudo_only_data1],
        pseudogene_data2$`Protein Accessions`[pseudogene_data2$`Peptide Sequence` %in% pseudo_only_data2]
      )
    )
  } else {
    # Create an empty data frame if no unique pseudogenes are found
    not_matched_pseudogenes_df <- data.frame(
      Peptide = character(0),
      Status = character(0),
      Protein_Name = character(0),
      Protein_Accessions = character(0)
    )
  }
  
  # Create a seperate dataframe with pseudogenes of one data frame for which 
  # the peptides also occur in the other peptide dataframe
  # essentially pseudogenes that are regular genes in the other db
  if (length(c(pseudo_data1_in_data2, pseudo_data2_in_data1)) > 0) {
    pseudo_in_other_as_normal <- data.frame(
      Peptide = c(pseudo_data1_in_data2, pseudo_data2_in_data1),
      Status = c(
        rep(paste0("Pseudogene in ", dp1, " but normal in ", dp2), length(pseudo_data1_in_data2)),
        rep(paste0("Pseudogene in ", dp2, " but normal in ", dp1), length(pseudo_data2_in_data1))
      ),
      Protein_Name = c(
        pseudogene_data1$`Protein Name`[pseudogene_data1$`Peptide Sequence` %in% pseudo_data1_in_data2],
        pseudogene_data2$`Protein Name`[pseudogene_data2$`Peptide Sequence` %in% pseudo_data2_in_data1]
      ),
      Protein_Accessions = c(
        pseudogene_data1$`Protein Accessions`[pseudogene_data1$`Peptide Sequence` %in% pseudo_data1_in_data2],
        pseudogene_data2$`Protein Accessions`[pseudogene_data2$`Peptide Sequence` %in% pseudo_data2_in_data1]
      )
    )
  } else {
    pseudo_in_other_as_normal <- data.frame(
      Peptide = character(0),
      Status = character(0),
      Protein_Name = character(0),
      Protein_Accessions = character(0)
    )
  }
  
  # create dataframe with pseudogenes that are exclusive for one database
  # pseudogene peptides that do not occur in the complete database of the other
  if (length(c(pseudo_exclusive_data_1, pseudo_exclusive_data_2)) > 0) {
    pseudo_df_exclusive <- data.frame(
      Peptide = c(pseudo_exclusive_data_1, pseudo_exclusive_data_2),
      Status = c(
        rep(paste0("Pseudogene in ", dp1, " but not at all in ", dp2), length(pseudo_exclusive_data_1)),
        rep(paste0("Pseudogene in ", dp2, " but not at all in ", dp1), length(pseudo_exclusive_data_2))
      ),
      Protein_Name = c(
        pseudogene_data1$`Protein Name`[pseudogene_data1$`Peptide Sequence` %in% pseudo_exclusive_data_1],
        pseudogene_data2$`Protein Name`[pseudogene_data2$`Peptide Sequence` %in% pseudo_exclusive_data_2]
      ),
      Protein_Accessions = c(
        pseudogene_data1$`Protein Accessions`[pseudogene_data1$`Peptide Sequence` %in% pseudo_exclusive_data_1],
        pseudogene_data2$`Protein Accessions`[pseudogene_data2$`Peptide Sequence` %in% pseudo_exclusive_data_2]
      )
    )
  } else {
    pseudo_df_exclusive <- data.frame(
      Peptide = character(0),
      Status = character(0),
      Protein_Name = character(0),
      Protein_Accessions = character(0)
    )
  }
  
  # Save results to CSV files
  dir.create(paste0("csv_export_", dp1, "_", dp2), showWarnings = FALSE)
  write.csv(matched_df, paste0("csv_export_", dp1, "_", dp2, "/matched_peptides_with_proteins_",dp1,"and",dp2,".csv"), row.names = FALSE)
  write.csv(not_matched_df,paste0("csv_export_", dp1, "_", dp2, "/non_matched_peptides_with_proteins_",dp1,"and",dp2,".csv"), row.names = FALSE)
  write.csv(pseudogene_matched_df, paste0("csv_export_", dp1, "_", dp2, "/matched_pseudogenes_", dp1, "_and_", dp2, ".csv"), row.names = FALSE)
  write.csv(not_matched_pseudogenes_df, paste0("csv_export_", dp1, "_", dp2, "/not_matched_pseudogenes_", dp1, "_and_", dp2, ".csv"), row.names = FALSE)
  write.csv(pseudo_in_other_as_normal, paste0("csv_export_", dp1, "_", dp2, "/pseudo_in_one_as_normal_in_other_",dp1,"_and_",dp2,".csv"), row.names = FALSE)
  write.csv(pseudo_df_exclusive, paste0("csv_export_", dp1, "_", dp2, "/pseudo_exclusive_",dp1,"_and_",dp2,".csv"), row.names = FALSE)
  
  # make a venn diagramm
  # we need to count total length of the dataset
  total_data1 <- length(peptide_data1$`Peptide Sequence`)
  total_data2 <- length(peptide_data2$`Peptide Sequence`)
  overlapping <-length(matched_df$Peptide)
  
  # generate databases output
  
  # Prepare data for Venn diagram
  # if statement to decide whether to show peptide numbers of pseudogenes
  if (show_pseudo_1 && show_pseudo_2) { # show pseudo for both
    venn_data <- list(
      a = unique(peptide_data1$`Peptide Sequence`),
      b = unique(peptide_data2$`Peptide Sequence`),
      c = unique(pseudogene_data1$`Peptide Sequence`),
      d = unique(pseudogene_data2$`Peptide Sequence`)
    )
    # rename the venn_data
    venn_data <- setNames(
      venn_data,
      c(
        paste0("Peptides_", dp1),
        paste0("Peptides_", dp2),
        paste0("Pseudogenes_", dp1),
        paste0("Pseudogenes_", dp2)
      )
    )
  } else if (show_pseudo_1) { # show pseudo only for data 1
    venn_data <- list(
      a = unique(peptide_data1$`Peptide Sequence`),
      b = unique(peptide_data2$`Peptide Sequence`),
      c = unique(pseudogene_data1$`Peptide Sequence`)
    )
    # rename the venn_data
    venn_data <- setNames(
      venn_data,
      c(
        paste0("Peptides_", dp1),
        paste0("Peptides_", dp2),
        paste0("Pseudogenes_", dp1)
      )
    )
  } else if (show_pseudo_2) { # show pseudo only for data 2
    venn_data <- list(
      a = unique(peptide_data1$`Peptide Sequence`),
      b = unique(peptide_data2$`Peptide Sequence`),
      c = unique(pseudogene_data2$`Peptide Sequence`)
    )
    # rename the venn_data
    venn_data <- setNames(
      venn_data,
      c(
        paste0("Peptides_", dp1),
        paste0("Peptides_", dp2),
        paste0("Pseudogenes_", dp2)
      )
    )
  } else { # show pseudo for none
    venn_data <- list(
      a = unique(peptide_data1$`Peptide Sequence`),
      b = unique(peptide_data2$`Peptide Sequence`)
    )
    # rename the venn_data
    venn_data <- setNames(
      venn_data,
      c(
        paste0("Peptides_", dp1),
        paste0("Peptides_", dp2)
      )
    )
  }
  
  
  # Plot Venn diagram using ggVennDiagram
  # (add "+ coord_flip()" before theme_minimal() to rotate by 90°)
  vennplot <- ggVennDiagram(
    venn_data,
    label_alpha = 0.5,
    category.names = names(venn_data),
    fill_color = c("blue", "red", "green", "purple")
  ) + theme_minimal()
  
  print("done")
  return(vennplot)
}

# comparison between Genbank and REFSEQ
print(create_vd_peptides(Genbank, wREFSEQ, FALSE, FALSE))

# comparison between db2 (genbank) including refseq annotated pseudogenes and 
# db5 (base db2 with pseudogense annotated by genbank)
print(create_vd_peptides(db2_wPseudogenes, db5, TRUE, TRUE))

# Uniprot reference and genbank (no pseudo)
print(create_vd_peptides(UPreference, Genbank, FALSE, FALSE))

# Uniprot reference and db2 (genbank, no pseudo)
print(create_vd_peptides(UPreference, db2, FALSE, FALSE))

# Uniprot reference and refseq annotation (no pseudo)
print(create_vd_peptides(UPreference, wREFSEQ, FALSE, FALSE))

# Uniprot reference and db2 (genbank with refseq pseudo)
print(create_vd_peptides(UPreference, db2_wPseudogenes, FALSE, TRUE))

# Genbank (no pseudo) and corrected genome (pgap, with pseudo)
print(create_vd_peptides(Genbank, newORF_wPseudo_v2, FALSE, TRUE))

# db2 (Genbank, no pseudo) and corrected genome (pgap, with pseudo)
print(create_vd_peptides(db2, newORF_wPseudo_v2, FALSE, TRUE))

# Uniprot reference and corrected genome (pgab, with pseudo)
print(create_vd_peptides(UPreference, newORF_wPseudo_v2, FALSE, TRUE))

library(patchwork)
print(
  create_vd_peptides(UPreference, Genbank, FALSE, FALSE) +
  create_vd_peptides(Genbank, newORF_wPseudo_v2, FALSE, TRUE) +
  create_vd_peptides(UPreference, newORF_wPseudo_v2, FALSE, TRUE)
)


# import the fasta headers RSAU
GFF <- read_tsv("Desktop/BIO253 Block Research cycle in genomics/Week 2 data files and R scripts /R stuff from Luuk /db2AND5_pseudo_SA6850_fastaHeadersParsedIntoDF (1).tsv")

# import the fasta headers pgaptmp
GFF <- read_tsv("Desktop/BIO253 Block Research cycle in genomics/Week 2 data files and R scripts /R stuff from Luuk /correctedFastaHeadersParsedIntoDF.tsv")

# Warning! use the right GFF for the right datasets. otherwise you get NAs!
# !! also watch out with the filter for RSAU or other calls in the function !!


# import the written csv files we want to investigate
pseudo_exclusive_db2_wPseudogenes_and_UPreference <- read_csv("pseudo_exclusive_db2_wPseudogenes_and_UPreference.csv")

# make a new function
# now we want to list the unique proteins from the pseudogenes into a histogram

# first make a useful function for GFF matches
match_locus_to_gff <- function(rsau_list, GFF) {
  # Ensure the necessary columns exist in the GFF data frame
  if (!all(c("locusTag", "SALocationStart") %in% colnames(GFF))) {
    stop("The GFF data frame must contain columns 'locusTag' and 'SALocationStart'.")
  }
  df_rsau_location <- data.frame(ProteinNames = as.character(rsau_list))
  df_rsau_location$startLocation <- numeric(nrow(df_rsau_location))
  
  for (i in seq_along(rsau_list)) {
    match_row <- which(GFF$locusTag == rsau_list[i])
    if (length(match_row) > 0) {
      df_rsau_location$startLocation[i] <- GFF$SALocationStart[match_row]
    } else {
      df_rsau_location$startLocation[i] <- NA
    }
  }
  
  # now we want to create a summary of the collected ProteinName data
  # how many peptides for each specific ProteinName
  # list the top 5 most expressed ProteinNames
  # Sum the start locations by ProteinNames
  df_rsau_location <<- df_rsau_location
  summed_ProteinNames_ <<- df_rsau_location %>%
    group_by(ProteinNames) %>%
    summarize(Peptide_count = n())
  
  # Sort by Peptide_count in descending order and select the top 5 rows
  top_5_ProteinNames <- summed_ProteinNames_ %>%
    arrange(desc(Peptide_count)) %>%
    head(5)
  
  # Print the result
  
  print(summary(summed_ProteinNames_))
  print(top_5_ProteinNames)
  
}

# Then make function for histogram. use GFF function in it
make_locus_histogram <- function(not_matched_data){
  # Ensure Protein_Name column exists in the input data
  if (!"Protein_Name" %in% colnames(not_matched_data)) {
    stop("The input data must contain a column named 'Protein_Name'")
  }
  
  # CHANGE VALUE "pgaptmp_" TO YOUR NEED!
  # Extract values that start with "pgaptmp_" followed by digits
  filtered_proteins <- na.omit(str_extract(not_matched_data$Protein_Name, "pgaptmp_\\d+"))
  # warning! not every headers has "pgaptmp_" annotation. Change in case of different annotations
  
  # Convert to a list if necessary
  protein_list <<- as.list(filtered_proteins) 
  # double arrow gives the result in the global environment
  print("list done")
  
  # use list now in GFF with the match_locus function
  
  check_locus <<- match_locus_to_gff(protein_list, GFF)
  print("dataframe done")
  
  # now ggplot the whole stuff
  
  plot <- ggplot(df_rsau_location, aes(x = startLocation, fill = ProteinNames))+
    geom_histogram(bins = 250) +
    xlab("")+
    ylab("")+
    ggtitle("")+
    theme(plot.title = element_text(size = 20))+  # Increase title size
    theme_minimal()# check if the plotter works properly, does not apply when function is used somehow
  print(plot)
  print("done")
}

make_locus_histogram(pseudo_exclusive_db2_wPseudogenes_and_UPreference)

# find out the stop and start sequence from the protein in the GFF

# mission:
# find the protein names of the unique pseudogene peptides, can be done with written csv file
# how many proteins are they actually? make a histogram and list these proteins/ peptides, 
# quantify unique peptides and also show start and stop of the protein


# make new analysis:

create_vd_peptides(newORF_wPseudo_v2, UPreference)

# import created csv for exclusive pseudogenes:
pseudo_exclusive_newORF_wPseudo_v2_and_UPreference <- read_csv("pseudo_exclusive_newORF_wPseudo_v2_and_UPreference.csv")

# important note! the protein names are not RSAU in this file, 
# change the function line: extract values that start with RSAU
# to values that start with pgaptmp_


# or add in the histogram function an if statement looking which ProteinNames method 
# you have and then apply directly the appropriate selection method

make_locus_histogram(pseudo_exclusive_newORF_wPseudo_v2_and_UPreference)



#################################

# experiment with the vd function, to maybe add more components to match to
# we want to find out, why uniprot has more than 700 unique peptides to itself
# in comparison to all other data bases
# find out what proteins they are from these peptides.




datasets <- list(UPreference, Genbank, wREFSEQ, db2, db2_wPseudogenes, db5)
dataset_names <- c("UPreference", "Genbank", "wREFSEQ", "db2", "db2_wPseudogenes", "db5")

# Initialize a data frame to store results
not_matched_df <- data.frame(
  Peptide = character(),
  Status = character(),
  Protein_Name = character(),
  Protein_Accessions = character()
)

# Compare each dataset to all others
for (i in seq_along(datasets)) {
  # Current dataset
  current_dataset <- datasets[[i]]
  current_name <- dataset_names[i]
  
  # Combine `Peptide Sequence` columns from other datasets
  other_peptides <- unlist(lapply(datasets[-i], function(df) df$`Peptide Sequence`))
  
  # Find peptides exclusive to the current dataset
  only_current <- setdiff(current_dataset$`Peptide Sequence`, other_peptides)
  
  # Extract rows for exclusive peptides
  unique_current <- current_dataset[current_dataset$`Peptide Sequence` %in% only_current, ]
  
  # Add exclusive peptides to the results data frame
  not_matched_df <- rbind(
    not_matched_df,
    data.frame(
      Peptide = unique_current$`Peptide Sequence`,
      Status = rep(current_name, nrow(unique_current)),
      Protein_Name = unique_current$`Protein Name`,
      Protein_Accessions = unique_current$`Protein Accessions`
    )
  )
}

# View results
print(not_matched_df)

# there are 784 filtered peptides uniquely for UniProt
# Remove rows with duplicate peptides? do not exist apparently

# find out how many proteins that are and what they are

library(dplyr)

# Filter peptides that are exclusively associated with "UPreference"
filtered_df <- not_matched_df %>%
  group_by(Peptide) %>%
  filter(all(Status == "UPreference")) %>%
  ungroup()

# View the result
print(filtered_df)


# Summarize the data frame to count occurrences of each unique protein name
protein_summary_UniProt <- filtered_df %>%
  group_by(Protein_Name) %>%
  summarize(
    Count = n(),  # Count occurrences of each protein
    Peptides = paste(unique(Peptide), collapse = "; ")  # Optional: list unique peptides associated
  ) %>%
  ungroup()

# Save the summarized data frame as a CSV file
write.csv(protein_summary_UniProt, "protein_summary_exclusive_UniProt.csv", row.names = FALSE)

# They are 545 Proteins from the 784 unique peptides identified

# now i want to check if these Proteins are also unique to uniprot, 
# or only the peptides corresponding to them

# List of other datasets to check against again the uniprot proteins
other_datasets <- list(Genbank, wREFSEQ, db2, db2_wPseudogenes, db5) 

# Combine all protein names from other datasets into a single vector
all_other_proteins <- unique(unlist(lapply(other_datasets, function(df) df$`Protein Name`)))

# Add a new column to the UniProt data frame
unique_proteins_Uniprot_or_not <- protein_summary_UniProt %>%
  mutate(
    Proteins_match = ifelse(
      Protein_Name %in% all_other_proteins,
      "matches with another dataset",
      "unique protein to uniprot"
    )
  )

# View the updated data frame
print(unique_proteins_Uniprot_or_not)

# now we filter out the Proteins that match with other datasets,
# and keep only the proteins that are unique to Uniprot

# Filter for unique proteins to UniProt
unique_proteins_to_UniProt_ <- unique_proteins_Uniprot_or_not %>%
  filter(Proteins_match == "unique protein to uniprot")

# View the result
print(unique_proteins_to_UniProt_)

# there are 528 unique proteins only in UniProt and no other database

# Save the Unique Proteins from Uniprot data frame as a CSV file
write.csv(unique_proteins_to_UniProt_, "Proteins_that_are_exclusively_in_UniProt.csv", row.names = FALSE)



############
