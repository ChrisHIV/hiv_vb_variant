# Author: Chris Wymant
# Acknowledgment: I wrote this while funded by ERC Advanced Grant PBDR-339251
# and a Li Ka Shing Foundation grant, both awarded to Christophe Fraser.
# This script analyses a whole-genome nucleotide alignment with the VB variant
# consensus, the consensus of other Dutch subtype-B in BEEHIVE, and HXB2. It
# splits the genome into equally sized segments, calculates how many positions
# in each segment are SNPs (single nucleotide polymorphisms) between the VB variant and the
# consensus of other Dutch subtype-B in BEEHIVE, and compares this to a null
# model. The null distribution is constructed as follows. Conditioning on the 
# total number of nucleotide changes, the number in each segment is taken to be
# multinomially distributed, with the probability for each bin being the mean
# evolutionary rate over positions inside that segment, normalised to 1. 

library(tidyverse)
library(seqinr)

# WARNING: the command below deletes all objects in memory in your current R
# session. Convenient when repeatedly re-running this script.
rm(list = ls())

# INPUT ------------------------------------------------------------------------

# Set seed for random number generation, for reproducibility 
set.seed(12345)

# Set the working directory to where this code lives on your machine
# setwd()

# The whole-genome nucleotide alignment with the VB variant consensus, the 
# consensus of other Dutch subtype-B in BEEHIVE, and HXB2, provided as
# supplementary data with the associated paper.
file_in_aln <- "nucleotide_VB_vs_DutchB_wHXB2.fasta"

# The 'normalisation' file, measuring between-host evolutionary rate as a
# function of position in the HIV genome (HXB2 coordinates).
file_in_normalisation <- "data_external/HIV_DistanceNormalisationOverGenome.csv"

# Files we'll create
file_out_snps <- "SNPs_genome_distribution.pdf"
file_out_snps_raw <- "SNPs_genome_distribution_raw.pdf"

num_segments <- 10L

num_multinomial_draws <- 200000L

plot <- FALSE

# MAIN -------------------------------------------------------------------------

stopifnot(file.exists(file_in_aln))
stopifnot(file.exists(file_in_normalisation))

# positions (wrt HXB2) between which we have sequence for both the lineage and
# the Dutch subtype B consensus (defined by the Gall et al 2012 primers)
pos_start <- 497L
pos_end <- 9496L

segment_boundaries <- seq(pos_start, pos_end, length.out = num_segments + 1L)
segment_boundaries[[1]] <- pos_start - 1L # Since bins are defined as (min,max]

# Percentiles of the null distribution that we'll calculate. Subsequent code
# expects these specific values.
percentiles <- c(0.5, 2.5, 25, 50, 75, 97.5, 99.5) 

# Read in the aln, convert bases to upper case
aln <- read.fasta(file_in_aln, seqtype = "DNA") %>%
  map(toupper)

# Check the three expected seq names are present
name_hxb2 <- "B.FR.83.HXB2_LAI_IIIB_BRU.K03455"
name_dutch_b <- "consensus_Dutch_subtypeB_inBEEHIVE"
name_lineage <- "consensus_VB_variant"
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

# Get the HXB2 position of aln positions where the lineage and dutch B seqs
# disagree. Check the bases are what we expect.
pos_snps <- integer()
bases_ok <- c("A", "C", "G", "T")
hxb2_pos_current <- 0L
for (pos in 1:length_aln) {
  if (seq_hxb2[[pos]] != "-") hxb2_pos_current <- hxb2_pos_current + 1L
  if (seq_dutch_b[[pos]] == seq_lineage[[pos]] ||
      seq_dutch_b[[pos]] == "N" ||
      seq_lineage[[pos]] == "N" ||
      seq_dutch_b[[pos]] == "-" ||
      seq_lineage[[pos]] == "-") {
    next
  }

  if (! seq_dutch_b[[pos]] %in% bases_ok) {
    stop(paste("Unexpected base", seq_dutch_b[[pos]], "in seq", name_dutch_b, 
               "at position", pos))
  }
  if (! seq_lineage[[pos]] %in% bases_ok) {
    # Manual hack, Y means C or T so that's a SNP
    if (seq_lineage[[pos]] != "Y" && seq_dutch_b[[pos]] != "A") {
    stop(paste("Unexpected base", seq_lineage[[pos]], "in seq", name_lineage, 
               "at position", pos))
    }
  }
  pos_snps[[length(pos_snps) + 1L]] <- hxb2_pos_current
}
stopifnot(all(pos_snps >= pos_start))
stopifnot(all(pos_snps <= pos_end))
num_snps <- length(pos_snps)

# Read in the normalisation file
df_normalisation <- read_csv(file_in_normalisation, col_types = cols(
  POSITION_WRT_HXB2 = col_integer(),
  MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS = col_double()
  )) %>%
  rename(pos = POSITION_WRT_HXB2,
         distance = MEDIAN_PAIRWISE_DISTANCE_BETWEEN_STANDARD_REFS)

# Add extra rows before the start and after the end of the normalisation file,
# which are within the [pos_start, pos_end] range we consider but where
# normalisation values are missing. Then check we have all the positions needed.
df_normalisation <- df_normalisation %>%
  bind_rows(tibble(pos = c(seq(pos_start, min(df_normalisation$pos) - 1L),
                           seq(max(df_normalisation$pos) + 1L, pos_end)),
                   distance = NA)) %>%
  arrange(pos)
stopifnot(identical(df_normalisation$pos,
                    pos_start:pos_end))

# Group the normalisation file into segments; take the mean normalisation in
# each segment, and then re-normalise so the segment values sum to 1.
df_normalisation <- df_normalisation %>%
  mutate(segment = as.character(cut(pos, breaks = segment_boundaries)))
df_segments <- df_normalisation %>%
  group_by(segment) %>%
  summarise(distance = mean(distance, na.rm = TRUE),
            segment_start = min(pos),
            segment_end = max(pos),
            .groups = "drop") %>%
  arrange(segment_start) %>%
  mutate(fraction_snps_expected = distance / sum(distance))

# Make a set of multinomial draws, with probabilities given by the normalisation
# value for each segment. i.e. simulate how many SNPs expected in each segment
# by chance: our null model.
snps_over_segments_simulation <- rmultinom(
  n = num_multinomial_draws,
  size = num_snps,
  prob = df_segments$fraction_snps_expected
)

# Calculate the desired percentiles for the null model.
df_segments[paste0("num_snps_percentile_", percentiles)] <-
  apply(snps_over_segments_simulation, 1, function(snp_counts) {
  quantile(snp_counts, probs = percentiles / 100)
  }) %>%
  t()

# Merge the actual SNP counts with the null model percentiles
snp_segment_counts <- table(cut(pos_snps, segment_boundaries)) 
stopifnot(identical(names(snp_segment_counts),
                    df_segments$segment))
df_segments <- full_join(df_segments,
                         tibble(segment = names(snp_segment_counts),
                                snp_count = snp_segment_counts),
                         by = "segment")

# Plot!
if (plot) {
p <- ggplot(df_segments) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_x_continuous(limits = c(pos_start, pos_end), breaks = 1:9 * 1000, expand = c(0, 0)) +
  geom_rect(aes(xmin = segment_start,
                xmax = segment_end,
                ymin = num_snps_percentile_25,
                ymax = num_snps_percentile_75), fill = "grey", alpha = 1) +
  geom_rect(aes(xmin = segment_start,
                xmax = segment_end,
                ymin = num_snps_percentile_2.5,
                ymax = num_snps_percentile_97.5), fill = "grey", alpha = 0.7) +
  geom_rect(aes(xmin = segment_start,
                xmax = segment_end,
                ymin = num_snps_percentile_0.5,
                ymax = num_snps_percentile_99.5), fill = "grey", alpha = 0.4) +
  geom_segment(aes(x = segment_start, xend = segment_end,
                   y = snp_count, yend = snp_count), colour = "black") +
  labs(x = "genome position (HXB2)",
       y = "number of nucleotide changes")
ggsave(file_out_snps,
       p, height = 4.5, width = 5.5)

p <- ggplot() +
  theme_classic() + 
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  scale_x_continuous(limits = c(pos_start, pos_end), breaks = 1:9 * 1000, expand = c(0, 0)) +
  geom_point(data = tibble(x = pos_snps, y = 0), aes(x, y), shape = "|", size = 1.5) +
  geom_line(data = df_normalisation, aes(x = pos, y = distance)) +
  labs(x = "genome position (HXB2)",
       y = "evolutionary rate (substitutions per site)")
ggsave(file_out_snps_raw,
       p, height = 4.5, width = 5.5)
}

# After visual inspection, when 10 segments are used, segment 5 (with 32 snps) 
# deviates most from the null. Calculate its two-tailed p value, unadjusted for
# multiple testing:
2 * sum(snps_over_segments_simulation[5,] >= 32) / num_multinomial_draws
