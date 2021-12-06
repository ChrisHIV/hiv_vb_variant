# Authors: Chris Wymant and Francois Blanquart 
# Acknowledgment: Chris wrote this while funded by ERC Advanced Grant PBDR-339251
# and a Li Ka Shing Foundation grant, both awarded to Christophe Fraser.
# Abbreviations: VL = viral load (log10 copies per ml), PC = principal component,
# df = dataframe

library(ape)
library(phytools)
library(tidyverse)

# WARNING: deletes all objects in your R session's memory
rm(list = ls()) 

# Change directory appropriately:
setwd("replace this bit by the path to your directory with the data")

# Dataframe with VLs and PCs. Loadinf takes a few seconds - 50MB 
df <- read_csv("BEEHIVE_VLs.csv")

tree <- read.tree("BEEHIVE_tree.tree")

# Test association with VL for the first N PCs, with N chosen to taste
number_pcs_to_test <- 50
df.pcs <- purrr::map(1:number_pcs_to_test, function(i) {
  lm_ <- lm(df[["BEEHIVE_LVL"]] ~ df[[paste0("PC", i)]])
  c(summary(lm_)[[4]][2, ], confint(lm_)[2,])
}) %>%
  bind_rows() %>%
  mutate(pc = factor(1:number_pcs_to_test)) %>%
  mutate(pc = fct_reorder(pc, `Pr(>|t|)`)) %>%
  mutate(is_lineage_pc = if_else(pc == "5", "lineage", "other"))

View(df.pcs)

# Plot PC effect sizes on VL, with 95% confidence intervals
p <- ggplot(df.pcs %>% 
              mutate(`2.5 %` = if_else(Estimate < 0, -1 * `2.5 %`, `2.5 %`),
                     `97.5 %` = if_else(Estimate < 0, -1 * `97.5 %`, `97.5 %`)) %>%
              mutate(Estimate = abs(Estimate))) +
  geom_point(aes(x = pc, y = Estimate, col = is_lineage_pc)) +
  geom_errorbar(aes(x = pc, ymin = `2.5 %`, ymax = `97.5 %`, col = is_lineage_pc)) +
  theme_classic() +
  geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
  labs(x = "Principal component (ordered by p value)",
       y = "Size of association with VL",
       col = "PC:")
p

# Choose a PC to plot on the tree. PC5 is the lineage.
PC <- "PC5"

# Multiply the PC by minus one, purely for consistent red/blue aesthetics with
# the paper. Set colours for the tree.
df[[PC]] <- -df[[PC]]
cols <- df[[PC]]
names(cols) <- df$id_paper

# Show the PC on the tree. We see PC5 is associated with a single long branch:
# this branch defines the lineage.
plotBranchbyTrait(ladderize(tree),
                  cols, 
                  mode = "tips", 
                  edge.width = 0.4, 
                  show.tip.label = F, 
                  title = PC)
