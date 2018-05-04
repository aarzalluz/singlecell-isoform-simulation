# LOAD REQUIRED LIBRARIES
library(readr)
library(dplyr)
library(ggplot2)

# LOAD REQUIRED CUSTOM FUNCTIONS
source("lr_simulation_functions.R")

# PREPARE DATA
# load and clean expression values to convert into number of reads per transcript
load("data/transcript_expression.rda") # neural stem cell and oligodendrocyte expression values for two replicates (data frame)
expr <- data.frame(transcript = quantification$PacBio, NSC = quantification$POST_NSC1, OLD = quantification$POST_OLD1) # select one replicate per sample for the simulation

# generate reduced-expression vectors for a range of cell numbers
cells <- c(2,6,10,16,20)
reduced_expression.NSC <- lapply(cells, reduce_expression, expression = expr$NSC) # using these vectors, calls to the simulate_long_reads()
reduced_expression.OLD <- lapply(cells, reduce_expression, expression = expr$OLD) # function were made, and reads stored in fasta files
reads <- as.integer(round(1000000 / cells)) # aprox. no. of reads per cell in simulation

# SUMMARIZE SIMULATION RESULTS
# simulated expression matrices
load("data/gene.isoform_table.rda") # load isoform-to-gene table (data frame)
expr_matrices.NSC <- list()
expr_matrices.OLD <- list()
# use list of reduced expression to make isoform expression matrices for each simulation
for (i in seq(length(cells))){
  # 1. filter out undetected isoforms
  # 2. order by gene and by expression
  expr_matrices.OLD[[i]] <- data.frame(gene.isoform_table, expression = reduced_expression.OLD[[i]]) %>% filter(expression!=0)
  expr_matrices.NSC[[i]] <- data.frame(gene.isoform_table, expression = reduced_expression.NSC[[i]]) %>% filter(expression!=0)
}

# SIMULATION RESULTS ANALYSIS
# see which genes show more than one isoform in each simulation, and count them
counts.NSC <- list()
counts.OLD <- list()
for (i in seq(length(expr_matrices.NSC))){
  counts.NSC[[i]] <- rle(expr_matrices.NSC[[i]]$associatedGene) %>% unclass %>% as.data.frame %>% filter(lengths!=1)
  counts.OLD[[i]] <- rle(expr_matrices.OLD[[i]]$associatedGene) %>% unclass %>% as.data.frame %>% filter(lengths!=1)
}
# count number of MIG per simulation and make data frame
mig.summary <- data.frame(cell.no = rep(cells, 2) %>% as.factor, 
                          MIG = c(lapply(counts.NSC, nrow) %>% unlist, lapply(counts.OLD, nrow) %>% unlist),
                          sample = c(rep("NSC", 5), rep("OLD", 5)))
mig.summary$sample <- factor(mig.summary$sample, labels = c("Neural stem cells", "Oligodendrocytes"))

iso_switch_genes <- list() # make empty list

# recursively call new function and store list of genes
for(i in seq(length(expr_matrices.NSC))){
  iso_switch_genes[[i]] <- compare_mig(expr_matrices.NSC, expr_matrices.OLD, counts.NSC, counts.OLD, i)
}
# count no. of podium changes per simulation and make data frame for plotting
iso_switch_counts <- lapply(iso_switch_genes, length) %>% unlist
pch.summary <- data.frame(cell.no = cells %>% as.factor, isoform_switches = iso_switch_counts)

# RESULT PLOTS
# plot number of MIG detected per simulation
mig.bulk <- 2423 # total no. of genes with >1 isoform in Tardaguila et al. 2018 data

ggplot(mig.summary, aes(x = cell.no, y = MIG, fill = sample)) +
  geom_bar(stat = "identity",  position = position_dodge(), colour = "black", width = 0.75, size = 2) +
    geom_hline(yintercept = mig.bulk, linetype = "dashed", size = 1, colour = "gray20") + 
  labs(x = "\n Number of single cells per Sequel run (total = 1 million reads) \n", 
       y = "Multi-isoform genes detected \n", fill = "Cell type: ") + ylim(0, 2500) +
  theme_minimal() +
  theme(text = element_text(family = "AvantGarde"), axis.text = element_text(size = 56), axis.title = element_text(size = 56), 
        legend.text = element_text(size = 40), legend.title = element_text(size = 40),
        legend.position = "bottom") +
  scale_x_discrete(breaks = c("2", "6", "10", "16", "20")) +
  annotate("text", x = 2, y = mig.bulk, vjust = -0.5, hjust = -0.2, label = "Original data: 2423 MIG", 
           size = 14, family = "AvantGarde", color = "gray20")

# plot number of isoform switches per simulation
pch.bulk <- 337   # number of isoform switches in Tardaguila et al. 2018 data

ggplot(pch.summary, aes(x = cell.no, y = isoform_switches)) + 
  geom_bar(stat = "identity", position = position_dodge(), colour = "black", width = 0.5, size = 2, fill = "#D2701F") +
  geom_hline(yintercept = pch.bulk, linetype = "dashed", size = 1, colour = "gray20") + 
  labs(x = "\n Number of single cells per Sequel run (total = 1 million reads) \n", 
       y = "Isoform switches detected \n") +
  theme_minimal() +
  theme(text = element_text(family = "AvantGarde"),  axis.text = element_text(size = 56), axis.title = element_text(size = 56)) +
  scale_x_discrete(breaks = c("2", "6", "10", "16", "20")) +
  annotate("text", x = 2, y = pch.bulk, vjust = 1.5, hjust = -0.2, label = "Original data: 337 isoform switches", 
           size = 14, family = "AvantGarde", color = "gray20")
