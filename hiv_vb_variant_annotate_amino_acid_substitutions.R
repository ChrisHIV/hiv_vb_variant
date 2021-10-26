# Author: Chris Wymant
# Acknowledgment: I wrote this while funded by ERC Advanced Grant PBDR-339251
# and a Li Ka Shing Foundation grant, both awarded to Christophe Fraser.
# This script collects codon-level annotations for the HXB2 HIV genome and
# transfers them to those positions where two other amino acid sequences (for
# example in the present study, the consensus of the VB variant sequences and
# the consensus of other Dutch subtype-B sequences in BEEHIVE) disagree with
# each other. 

# Abbreviations:
# AA = amino acid
# df = dataframe
# seq = sequence
# aln = alignment

library(tidyverse)
library(seqinr)

# WARNING: the command below deletes all objects in memory in your current R
# session. Convenient when repeatedly re-running this script.
rm(list = ls())

# INPUT ------------------------------------------------------------------------

# Set the working directory to where this code lives on your machine
# setwd()

# The per-gene amino acid alignments with the VB variant consensus, the 
# consensus of other Dutch subtype-B in BEEHIVE, and HXB2, provided as
# supplementary data with the associated paper.
files_in_aln <- Sys.glob("VB_vs_DutchB_*.fasta")

file_in_hxb2_annotated_codons <- "data_external/hxb2_annotated_SpecialSites_ByCodon.csv"

file_out_vb_differences <- "VB_vs_Dutch_B_subsitutions_annotated.csv"

# MAIN -------------------------------------------------------------------------

stopifnot(file.exists(file_in_hxb2_annotated_codons))
stopifnot(length(files_in_aln) > 0L)


# Label alns by their gene
genes <- str_match(files_in_aln, "VB_vs_DutchB_([a-z]{3}).fasta")[, 2]
names(files_in_aln) <- genes

# Read in the HXB2 annotation file. Check that each gene-codon is unique.
df_hxb2 <- read_csv(file_in_hxb2_annotated_codons, col_types = cols(
  codon_in_that_gene = col_integer()
)) %>%
  rename(codon = codon_in_that_gene)
stopifnot(!anyDuplicated(df_hxb2 %>% select(gene, codon)))

# Iterate through genes, reading in that alignment, recording codons which
# differ between the VB variant and the consensus of other Dutch subtype B seqs
# in BEEHIVE, and transferring HXB2 annotations to those positions. Merge the
# results from all genes into a single df.
df_differences <- map(genes, function(gene_) {
  
  cat("Analysing alignment for", gene_, "\n")
  
  # Read in the aln, convert seqs to upper case
  aln <- read.fasta(files_in_aln[[gene_]], seqtype = "AA") %>%
    map(toupper)

  # Check the three expected seq names are present
  name_hxb2 <- "HXB2"
  name_dutch_b <- "DutchSubtypeB"
  name_lineage <- "VB"
  stopifnot(name_hxb2 %in% names(aln))
  stopifnot(name_dutch_b %in% names(aln))
  stopifnot(name_lineage %in% names(aln))
  
  seq_hxb2    <- aln[[name_hxb2]]
  seq_dutch_b <- aln[[name_dutch_b]]
  seq_lineage <- aln[[name_lineage]]
  
  # Check all seqs have the same length
  length_aln <- length(seq_hxb2)
  stopifnot(length_aln == length(seq_dutch_b))
  stopifnot(length_aln == length(seq_lineage))
  
  # Check the HXB2 seq in the aln is what we expected
  seq_hxb2_ungapped <- seq_hxb2[seq_hxb2 != "-"]
  seq_hxb2_expected <- df_hxb2 %>%
    filter(gene == gene_) %>%
    arrange(codon) %>%
    pull(HXB2_amino_acid)
  stopifnot(identical(seq_hxb2_ungapped,
                      seq_hxb2_expected))
  
  # Get the HXB2 position of aln positions where the lineage and dutch B seqs
  # disagree
  aln_pos_differences <- integer()
  hxb2_pos_differences <- integer()
  hxb2_pos_current <- 0L
  for (pos in 1:length_aln) {
    if (seq_hxb2[[pos]] != "-") hxb2_pos_current <- hxb2_pos_current + 1L
    if (seq_dutch_b[[pos]] == seq_lineage[[pos]] ||
        seq_dutch_b[[pos]] == "X" ||
        seq_lineage[[pos]] == "X") {
      next
    }
    aln_pos_differences[[length(aln_pos_differences) + 1L]] <- pos
    hxb2_pos_differences[[length(hxb2_pos_differences) + 1L]] <- hxb2_pos_current
  }

  # Record the AAs and positions of disagreements
  df_gene <- tibble(gene = gene_,
                    gene_position_hxb2 = hxb2_pos_differences,
                    aa_hxb2 = seq_hxb2[aln_pos_differences],
                    aa_DutchB_inBEEHIVE = seq_dutch_b[aln_pos_differences],
                    aa_vb = seq_lineage[aln_pos_differences])
  stopifnot(all(! df_gene$aa_DutchB_inBEEHIVE == df_gene$aa_vb))
  
  # Merge in annotation details 
  df_gene <- left_join(df_gene,
                       df_hxb2 %>% filter(gene == gene_) %>% select(-gene),
                       by = c("gene_position_hxb2" = "codon"))
  stopifnot(all(df_gene$aa_hxb2 == df_gene$HXB2_amino_acid |
                  df_gene$aa_hxb2 == "-"))
  
  df_gene %>% select(-HXB2_amino_acid)
  
}) %>%
  bind_rows() %>%
  arrange(gene)

# Categorise differences either substitutions or indels
df_differences <- df_differences %>% 
  mutate(difference_type = if_else(aa_DutchB_inBEEHIVE != "-" & aa_vb != "-", 
                                   "substitution", "indel"))
table(df_differences$difference_type)

# Restrict the set of CTL escapes to amino acids that the VB variant has, i.e.
# not just any amino acid at a site where the variant has some mutation.
max_num_escapes_per_codon <- df_differences %>% 
  pull(CTL_escape_details_for_q_less_than_1percent) %>%
  str_count(";") %>%
  max(na.rm = TRUE) + 1
df_differences_specific_escapes <- df_differences %>% 
  select(gene, gene_position_hxb2, aa_vb, CTL_escape_details_for_q_less_than_1percent) %>%
  separate(CTL_escape_details_for_q_less_than_1percent,
           into = paste0("Escape_", 1:max_num_escapes_per_codon),
           sep = "; ", fill = "right") %>%
  pivot_longer(starts_with("Escape_"), values_to = "escape") %>%
  filter(!is.na(escape)) %>%
  select(-name) %>%
  filter(substr(escape, 1, 1) == aa_vb) 

# Separate the mutations into adaptive (the amino acid is enriched in the 
# presence of the noted HLA type, indicated by a positive log odds for escape)
# vs nonadaptive (the amino acid is depleted in the presence of the noted HLA
# type, indicated by a negative log odds for escape).
df_differences_specific_escapes <- df_differences_specific_escapes %>%
  mutate(adaptive =
           str_detect(escape, "and escape log odds [0-9]*\\.{0,1}[0-9]+ for HLA")) %>%
  group_by(gene, gene_position_hxb2, aa_vb, adaptive) %>%
  summarise(escapes = paste(escape, collapse = "; "), .groups = "drop") %>%
  mutate(escapes_adaptive    = if_else( adaptive, escapes, NA_character_),
         escapes_nonadaptive = if_else(!adaptive, escapes, NA_character_)) %>%
  select(-c("adaptive", "escapes")) %>%
  group_by(gene, gene_position_hxb2, aa_vb) %>%
  summarise(escapes_adaptive = paste(na.omit(escapes_adaptive), collapse = ""),
            escapes_nonadaptive = paste(na.omit(escapes_nonadaptive), collapse = ""),
            .groups = "drop") %>%
  mutate(escapes_adaptive = if_else(escapes_adaptive == "", NA_character_,
                                    escapes_adaptive),
         escapes_nonadaptive = if_else(escapes_nonadaptive == "", NA_character_,
                                    escapes_nonadaptive))

# Merge the VB-specific CTL mutations back into the table with all mutations
nrow_df_differences_before_merge <- nrow(df_differences)
df_differences <- left_join(df_differences, df_differences_specific_escapes,
                            by = c("gene", "gene_position_hxb2", "aa_vb"))
if (nrow(df_differences) != nrow_df_differences_before_merge) {
  stop(paste0("Problem merging VB-specific CTL escape info into the df with all info"))
}

# Write output
write_csv(df_differences %>% select(gene,
                                    gene_position_hxb2,
                                    HXB2_base_positions,
                                    aa_hxb2,
                                    aa_DutchB_inBEEHIVE,
                                    aa_vb,
                                    difference_type,
                                    DRM_Los_Alamos,
                                    DRM_IAS,
                                    Other_features,
                                    CTL_escape_details_for_q_less_than_1percent,
                                    escapes_adaptive,
                                    escapes_nonadaptive) %>%
            replace_na(list(DRM_Los_Alamos = "",
                            DRM_IAS = "",
                            Other_features = "",
                            CTL_escape_details_for_q_less_than_1percent = "",
                            escapes_adaptive = "",
                            escapes_nonadaptive = "")) %>%
            rename(CTL_escape_all = CTL_escape_details_for_q_less_than_1percent,
                   CTL_escape_vb_adaptive = escapes_adaptive,
                   CTL_escape_vb_nonadaptive = escapes_nonadaptive),
          file_out_vb_differences)

# Print some numbers of mutations of interest
num_escapes <- df_differences %>% 
  mutate(escapes_either = (!is.na(escapes_adaptive)) | (!is.na(escapes_nonadaptive))) %>%
  pull(escapes_either) %>%
  sum()
num_escapes
num_escapes_both_directions <- df_differences %>% 
  mutate(escapes_both = (!is.na(escapes_adaptive)) & (!is.na(escapes_nonadaptive))) %>%
  pull(escapes_both) %>%
  sum()
num_escapes_adaptive <- sum(!is.na(df_differences$escapes_adaptive))
num_escapes_adaptive
num_escapes_nonadaptive <- sum(!is.na(df_differences$escapes_nonadaptive))
num_escapes_nonadaptive
num_mutations_at_ctl_site <- sum(!is.na(df_differences$CTL_escape_details_for_q_less_than_1percent))
num_mutations_at_ctl_site

# What are the HLA types for which the VB variant has adaptations?
escapes_adaptive_hlas <- df_differences$escapes_adaptive %>%
  na.omit() %>%
  str_split("; ") %>%
  unlist()
escapes_adaptive_hlas <- str_match(escapes_adaptive_hlas, pattern = "for HLA (.+)")[,2]
escapes_adaptive_hlas
table_ <- table(escapes_adaptive_hlas)
tibble(hla = names(table_),
       count = table_) 
