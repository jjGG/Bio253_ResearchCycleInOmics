# BIO 253 Research cycle in genomics

Collection of scripts that are useful for this course.

Contains scripts for the following tasks:

- Correct genome with mutations in `.vcf` file: [`~/R/operations on genome/correctMutationFromVCFToFastaSequence.R`](/R/operations%20on%20genome/correctMutationFromVCFToFastaSequence.R)
- Translate nucleotide fasta according to annotation: [`~/R/operations on genome/translateAnnotatedPseudogenesToProteins.R.`](/R/operations%20on%20genome/translateAnnotatedPseudogenesToProteins.R)
- Visualize sets of peptides: [`~/R/visualization/Pseudogene_peptides_comparison_and_visualization.R`](/R/visualization/Pseudogene_peptides_comparison_and_visualization.R)
- Create upset diagram out of sets of peptides: [`~/R/visualization/upset_diagram.R`](/R/visualization/upset_diagram.R)
- Compare fasta sequence to detect if they are different: [`~/Python/compare_sequence.py`](/Python/compare_sequence.py)
- Extract sequences from fasta file by ID e.g for a blast search [`~/Python/extract_region_from_fasta.py`](/Python/extract_region_from_fasta.py)
- Filter vcf file for samples of interest and annotate the SNPs with gene IDs [`~/Python/filtervcf_byInds_andAnnotate.py`](/Python/filtervcf_byInds_andAnnotate.py)
- Create a PCA using genotype data and color the points according to any phenotypic data column [`~/R/genotype_PCA.r`](/R/genotype_PCA.r)
- Calculate nucleotide diversity and estimate upgma tree using those estimates [`~/R/Calculate_Diversity_And_Plot_A_Tree.R`](/R/Calculate_Diversity_And_Plot_A_Tree.R)
