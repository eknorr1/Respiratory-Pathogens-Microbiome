# Canine nosomes: Assemble and clean 16S data 
#
# April 22, 2019

rm(list = ls())
graphics.off()



# Libraries ---------------------------------------------------------------

library(phyloseq)




# Working directories ------------------------------------------------------

path16Srun1 <- "~/Dropbox/Research/Active/Canine nosomes/Data/Microbiome/run1_output_f2_r2/"
path16Srun2 <- "~/Dropbox/Research/Active/Canine nosomes/Data/Microbiome/run2_output_f2_r2/"
pathOutput <- "~/Dropbox/Research/Active/Canine nosomes/Data/Assembled/"




# Set SampleIndices that link these samples back to the rest of the data -----
#
# Also set SampleIDs = rownames for the sequence matrices, which overwrite the values they had on import
# These also include samples generated at CGRB ie blanks
# Not every 'sample' downstream has a sample index.
#
# We did the 16S data in two runs (TODO: insert details here including CGRB ref numbers)
# In the first run we submitted ran samples 1 to 80. Sample 81 was a CGRB blank of some kind.
# In the second run we submitted samples 81 to 149. Sample 150 was a PBS extraction blank. 151-154 were NTCs from library prep. 
# The blanks have SampleID numbers NA 

SampleIndex <- c(1:80, NA, 81:149, NA, NA, NA, NA) 
RunID <- c(rep(1,81),rep(2,73))

sample_id1 <- paste(c(rep('s', 80), rep('b',1)), 1:81, sep = "")
sample_id2 <- paste(c(rep('s', 69), rep('b',4)), 81:153, sep = "")
sample_id <- c(sample_id1, sample_id2)

# Import 16S data and convert to phyloseq objects -------------------------
# The rownames are where downstream methods sometimes look for sample_id

# First half of samples
setwd(path16Srun1)

seqmat1 <- readRDS('seqtab_nochim.rds')
taxmat1 <- readRDS('tax.rds')

rownames(seqmat1) <- sample_id1

seqtab1 <- otu_table(seqmat1, taxa_are_rows = FALSE)
taxtab1 <- tax_table(taxmat1)

ps1 <- phyloseq(seqtab1, taxtab1)


# Second half of samples
setwd(path16Srun2)

seqmat2 <- readRDS('seqtab_nochim.rds')
taxmat2 <- readRDS('tax.rds')

rownames(seqmat2) <- sample_id2

seqtab2 <- otu_table(seqmat2, taxa_are_rows = FALSE)
taxtab2 <- tax_table(taxmat2)

ps2 <- phyloseq(seqtab2, taxtab2)





# Merge data to a single phyloseq object ----------------------------------

ps <- merge_phyloseq(ps1, ps2)





# Add in sample data, which for now is just sample and run IDs ------------

sampledata <- data.frame(SampleID = sample_id, SampleIndex = SampleIndex, RunID = RunID)
row.names(sampledata) <- sample_id

sample_data(ps) <- sampledata





# Attempt to remove potentially spurious ASVs -----------------------------

# 1. Remove those found in fewer than k samples 

k <- 3
isPresent <- colSums(otu_table(ps)>0)
ps <- prune_taxa(isPresent > k, ps)


# 2. Remove those found fewer than k times across all samples

k <- 20
ps <- prune_taxa(taxa_sums(ps) > k, ps) 








# Save --------------------------------------------------------------------

sxs <- ps

setwd(pathOutput)
save(sxs, file = "16S.Rdata")



