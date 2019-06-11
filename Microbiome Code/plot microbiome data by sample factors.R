# Canine nosomes: Plot 16S data by sample factors
#
# April 23, 2019

rm(list = ls())
graphics.off()


library(RColorBrewer)
library(phyloseq)
library(metacoder)


# Script parameters -------------------------------------------------------

pathToAssembledData <- "~/Dropbox/Research/Active/Canine nosomes/Data/Assembled/"







# Load and convert data ---------------------------------------------------

setwd(pathToAssembledData)
load("16S_combined.Rdata")
tm <- parse_phyloseq(sxs)
obj <- tm








# Some cleanup ------------------------------------------------------------

# Zero then remove low counts
obj$data$otu_table <- zero_low_counts(obj, data = "otu_table", min_count = 5)
no_reads <- rowSums(obj$data$otu_table[, 2:ncol(obj$data$otu_table)]) == 0
obj <- filter_obs(obj, "otu_table", !no_reads, drop_taxa = TRUE)


# Drop lower taxonomic orders, but keep their supertaxa
obj <- taxa::filter_taxa(obj, taxon_ranks == "Family", supertaxa = TRUE)                   



# Remove "odd" taxa
obj <- taxa::filter_taxa(obj, grepl(pattern = "^[a-zA-Z]+$", taxon_names))               








# Rarefaction to proportions ----------------------------------------------

obj$data$otu_table <- calc_obs_props(obj, data = "otu_table")





# Taxon abundance ---------------------------------------------------------

obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table",
                                       cols = 2:ncol(obj$data$otu_table))






# Generate comparisons ----------------------------------------------------

is_included <- obj$data$sample_data$SampleType == "nasal"
treatment_groups <- obj$data$sample_data$Symptomatic[is_included]
selected_samples <- obj$data$sample_data$sample_id[is_included]

obj$data$diff_table <- compare_groups(obj,
                                      data = "tax_abund",
                                      cols = colnames(obj$data$tax_abund) %in% selected_samples,
                                      groups = treatment_groups)





# Work with p values on comparisons ---------------------------------------

obj$data$diff_table$wilcox_p_value_adjusted <- p.adjust(p = obj$data$diff_table$wilcox_p_value, method = 'fdr')

obj$data$diff_table$mean_diff_adjusted <- obj$data$diff_table$mean_diff
obj$data$diff_table$mean_diff_adjusted[obj$data$diff_table$wilcox_p_value_adjusted > 0.2] <- 0 

obj$data$diff_table$log2_median_ratio_adjusted <- obj$data$diff_table$log2_median_ratio
obj$data$diff_table$log2_median_ratio_adjusted[obj$data$diff_table$wilcox_p_value_adjusted > 0.2] <- 0

#obj$data$diff_table$mood <- obj$data$diff_table$log2_median_ratio * (1 - obj$data$diff_table$wilcox_p_value_adjusted)
#obj$data$diff_table$mood[obj$data$diff_table$mood == Inf] <- max(obj$data$diff_table$mood[is.finite(obj$data$diff_table$mood)])
#obj$data$diff_table$mood[obj$data$diff_table$mood == -Inf] <- min(obj$data$diff_table$mood[is.finite(obj$data$diff_table$mood)])
#obj$data$diff_table$mood[!is.finite(obj$data$diff_table$mood)] <- 0




# Plot 2-colored heat tree ------------------------------------------------

node_color_range <- c("seagreen", "grey", "orange")

setwd("~/Dropbox/Research/Active/Canine nosomes/Figures/")
heat_tree(obj,
         node_label = taxon_names,
         node_size = n_obs(obj, "diff_table"),
         node_color = mean_diff,
         node_size_interval = c(10,200),
         node_size_range = c(0.008, .08),
         node_color_range = node_color_range,
         node_color_interval = c(-1, 1) * max(abs(obj$data$diff_table$mean_diff), na.rm = T),
         node_size_axis_label = "Number of taxa",
         node_color_axis_label = "Mean difference",
         node_label_size_range = c(0.01, 0.05),
         output_file = 'symptomatic heattree.pdf')



# Plot heat tree that colors symptomatic ----------------------------------

node_color_range <- c("grey", "orange", "red", "darkred")

obj$data$diff_table$mean_diff_pos <- obj$data$diff_table$mean_diff
obj$data$diff_table$mean_diff_pos[obj$data$diff_table$mean_diff_pos < 0] <- 0

heat_tree(obj,
          node_label = taxon_names,
          node_size = n_obs(obj, "diff_table"),
          node_size_interval = c(10,200),
          node_color = obj$data$diff_table$mean_diff_pos,
          node_size_range = c(0.008, .08),
          node_color_range = node_color_range,
          node_color_interval = c(0, 1) * max(abs(obj$data$diff_table$mean_diff), na.rm = T),
          node_size_axis_label = "Number of taxa",
          node_color_axis_label = "Mean difference",
          node_label_size_range = c(0.01, 0.05),
          output_file = 'symptomatic heattree pos.pdf')





# Plot heat tree that colors asymptomatic ---------------------------------

node_color_range <- c("grey", "seagreen", "green")

obj$data$diff_table$mean_diff_neg <- -obj$data$diff_table$mean_diff
obj$data$diff_table$mean_diff_neg[obj$data$diff_table$mean_diff_neg < 0] <- 0

heat_tree(obj,
          node_label = taxon_names,
          node_size = n_obs(obj, "diff_table"),
          node_size_interval = c(10,200),
          node_color = obj$data$diff_table$mean_diff_neg,
          node_size_range = c(0.008, .08),
          node_color_range = node_color_range,
          node_color_interval = c(0, 1) * max(abs(obj$data$diff_table$mean_diff), na.rm = T),
          node_size_axis_label = "Number of taxa",
          node_color_axis_label = "Mean difference",
          node_label_size_range = c(0.01, 0.05),
          output_file = 'symptomatic heattree neg.pdf')






