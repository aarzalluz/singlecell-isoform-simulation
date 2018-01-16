# LOAD REAQUIRED LIBRARIES
library(readr)
library(polyester)
library(dplyr)
library(hash)
library(ggplot2) 
library(scales)
library(ggradar)

# LOAD REQUIRED CUSTOM FUNCTIONS
source("sr_simulation_functions.R")

# LOAD TRANSCRIPTOME (not fasta formatted, data frame with headers and seqs)
load("data/sr_transcriptome.rda")

# MAKE TRANSRIPT SUBSTRINGS
for(i in c(100, 200, 300, 500, 1000)){
  subset_tr(tr_set = transcriptome, fragment_length = i, 
            filename = paste0("fragmented.", as.character(i), "nt.pacbio",
                   selected_end = "3end"))
}

# SIMULATE SHORT READS FROM FRAGMENTED TRANSCRIPTS
# load simulation parameters selected for every trimmed transcript length
load("data/sr_sim.parameters.rda")
# create fold_changes vector for simulate_experiment (polyester), a single replicate and a single group (just one sample)
fold_changes <- rep(1, length(transcriptome)/2)
# generate no. of reads per transcript values
load("data/transcript_expression.rda") # neural stem cell and oligodendrocyte expression values for two replicates (data frame)
NSC <- round(quantification$POST_NSC1) # select one replicate and convert expr values to whole numbers
OLD <- round(quantification$POST_OLD1)
# run simulate_experiment for all fragmented transcriptomes (generate simulated UMI reads)
# NOTE: here we present the code for NSC analysis, which can be rerun using oligodendrocyte data (OLD expression vector)
for(i in seq(nrow(sim_parameters))){
  message(paste("Simulating reads for", as.character(sim_parameters$length[i]), "nt fragments...")) # print progress
  simulate_experiment(fasta = paste0("fragmented.", as.character(sim_parameters$length[i]), "nt.pacbio.fa"), 
                      outdir = paste0("NSC_", as.character(sim_parameters$length[i]), "nt"), 
                      num_reps = 1, reads_per_transcript = NSC, readlen = sim_parameters$read_length[i],
                      fold_changes = fold_changes, paired = sim_parameters$paired[i])
}
# run simulate_experiment for unfragmented transcriptome (generate simulated full-length reads)
simulate_experiment(fasta = "transcriptome/transcriptome.fasta", 
                    outdir = "NSC_full_length", num_reps = 1, reads_per_transcript = NSC, readlen = 250,
                    fold_changes = fold_changes, paired = TRUE)

#-----------------------------------------------------------
# align and quantify isoform expression using RSEM + STAR
#-----------------------------------------------------------

# PROCESS RSEM RESULTS
# load and process REFERENCE data (i.e. quantification of reads from full-length transcripts)
# filename as output by RSEM (file not provided)
NSC_h <- process_sim_results("NSC_full.isoforms.results")

# load and process SIMULATED data (i.e. quantification of reads from fragmented transcripts)
# filenames output by RSEM are written as SAMPLE_NUMBERnt.isoforms.results
ids.n <- paste0("NSC_", c(100, 200, 300, 500, 1000), "nt")
sim.results.NSC <- list()
for (i in seq(length(ids.n))){
  sim.results.NSC[[i]] <- process_sim_results(paste0(ids.n[i], ".isoforms.results"))  # each position in the list stores a hash with the results of a simulation
}

# MAKE REFERENCE HASHES
  # NO. OF ISOFORMS PER MIG IN FULL-LENGTH QUANTIFICATION:
  # preserve NSC_h object, and create new reference hash to be modified by the make_mig_ref function
  # mig.NSC contains the no. of isoforms detected per MIG when quantifying with full-length simulated reads
mig.NSC <- process_sim_results("NSC_full.isoforms.results") %>% make_mig_ref 
  # NO. OF ISOFORMS PER MIG IN THE ANNOTATION:
  # load the entire set of gene-to-isoform relationships annotated in the transcriptome
load("data/gene.isoform_table.rda")
annot_count <- rle(gene.isoform_table$associatedGene) # count no. of isoforms per gene in the annotation
mig.annot <- make_mig_ref(hash(annot_count$values, annot_count$lengths)) # make hash of MIG in annotation

#------------------------------------------------------------------------------------------------------------------------------
# NOTE: The bit of code above is used to load and process the results of the RSEM run on the simulated and full-length reads. 
# Alternatively, the already processed data can be loaded:
load("data/sr_sim.isoform.results.NSC.rda") # load sim.results.NSC list of hashes
load("data/sr_full-length.isoform.results.NSC.rda") # load NSC_h hash
load("data/sr_full-length.reference.hash.rda") # load mig.NSC, full-length reference hash to compare UMI simulations to
load("data/sr_annot.reference.hash.rda") # load mig.annot, annotation reference hash to compare full-length simulation to
#------------------------------------------------------------------------------------------------------------------------------

# ANALYSE NO. OF ISOFORMS DETECTED PER MULTI-ISOFORM GENE
  # UMI simulations (reads from fragmented transcripts) vs full-length simulation (reads from across entire transcripts)
resolution_pcnt.NSC <- lapply(sim.results.NSC, calc_MIG_percent, mig.NSC) # resolution percentage for each MIG
  # binning of genes according to percentage of resolution intervals
tmp25n <- lapply(lapply(resolution_pcnt.NSC, as.data.frame), sub_percent, 0, 25) %>% unlist()
tmp50n <- lapply(lapply(resolution_pcnt.NSC, as.data.frame), sub_percent, 25, 50) %>% unlist()
tmp75n <- lapply(lapply(resolution_pcnt.NSC, as.data.frame), sub_percent, 50, 75) %>% unlist()
tmp100n <- lapply(lapply(resolution_pcnt.NSC, as.data.frame), sub_percent, 75, 100) %>% unlist()

  # full-length simulation vs annotated number of isoforms per MIG
resolution_pcnt.NSCannot <- calc_MIG_percent(NSC_h, mig.annot) %>% as.data.frame()
tmp25a <- sub_percent(resolution_pcnt.NSCannot, 0, 25)
tmp50a <- sub_percent(resolution_pcnt.NSCannot, 25, 50)
tmp75a <- sub_percent(resolution_pcnt.NSCannot, 50, 75)
tmp100a <- sub_percent(resolution_pcnt.NSCannot, 75, 100)
  
# SUMMARIZE IN DATA FRAME
pcnt_intervals.NSC <- data.frame(rep(c(100, 200, 300, 500, 1000, 7000), 4), 
                                 c(c(tmp25n, tmp25a), c(tmp50n, tmp50a), c(tmp75n, tmp75a), c(tmp100n, tmp100a)), 
                                 c(rep(4,6), rep(3,6), rep(2,6), rep(1,6)),
                                 c(rep(1, 5), 2))
colnames(pcnt_intervals.NSC) <- c("nucleotides", "percentage", "interval", "library")
pcnt_intervals.NSC$interval <- factor(pcnt_intervals.NSC$interval, 
                                      labels = rev(c("0-25%", "25-50%", "50-75%", "75-100%")))
pcnt_intervals.NSC$library <- factor(pcnt_intervals.NSC$library, labels = c("UMI-based", "Smart\n-seq"))
pcnt_intervals.NSC$percentage <- as.integer(pcnt_intervals.NSC$percentage)
pcnt_intervals.NSC$nucleotides <- factor(pcnt_intervals.NSC$nucleotides, 
                                         labels = c("100", "200", "300", "500", "1000", "Full-length"))

# PLOT RESULTS
# plot percentage of resolution intervals for UMI vs full-length simulations and full-length simulation vs annotation
ggplot() + ggtitle("Neural stem cells") + 
  geom_bar(data = pcnt_intervals.NSC, aes(fill = interval, x = nucleotides, y = percentage), width = 0.75,
           position = "fill", stat = "identity", colour = "black") +
  labs(x = "Nucleotides sequenced from 3' end", y = "% of MIG in sample") +
  scale_y_continuous(labels = percent_format()) + 
  theme(text = element_text(family = "AvantGarde"),
        axis.text = element_text(size = 14), axis.title = element_text(size = 20), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 16),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size= 14),
        strip.background = element_rect(colour = "black")) +
  scale_fill_manual(name = "Isoform resolution per MIG", values = c("#0C6A2B", "#8D1D1D", "#D2701F", "#044771")) +
  facet_grid(. ~ library, scales = "free", space = "free")

# qualitative spider graph of UMI-based, Smart-based and long-read's performance in the different requirements for
# isoform expression studies
spider_data <- data.frame(technology = c("UMI-based", "Smart-based", "Long reads"), depth = c(0.5, 1, 0), quantification = c(1, 0.5, 0),
                          errors = c(0, 0, 1), cells = c(1, 0.5, 0), isoforms = c(0, 0.5, 1))

png("spider_graph.png", width = 1500, height = 1000)  # save plot with right dimensions
ggradar(spider_data, 
        axis.labels = c("Sequencing depth", "Expression \nquantification", "Sequencing errors", 
                        "No. of cells sequenced", "No. of \nisoforms \ndetected"), 
        axis.label.size = 10, legend.text.size = 22, grid.label.size = 0,
        label.gridline = FALSE)
dev.off()

# NOTE: a comparative plot of the simulations from the 3' and 5' ends can be produced by loading a data frame where results for 
# the 5' end are also included (see figure 4 in our manuscript). These results can be generated by running this script on fragments 
# generated from the 5' end of transcripts using NSC expression values, and manually combining both data frames.
load("data/sr_pcnt.interval.results.both.ends.rda")

ggplot() + ggtitle("Neural stem cells") + 
  geom_bar(data = pcnt_intervals.NSC.3and5end, aes(fill = interval, x = nucleotides, y = percentage), width = 0.75, 
           position = "fill", stat = "identity", colour = "black") +
  labs(x = "No. of nucleotides sequenced", y = "% of MIG in sample") +
  scale_y_continuous(labels = percent_format()) +
  theme(text = element_text(family = "AvantGarde"), axis.text = element_text(size = 28), 
        axis.title.x = element_text(size = 30), axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 28), legend.title = element_text(size = 28),
        #legend.position = "bottom",
        plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),
        strip.text.x = element_text(size= 30), strip.text.y = element_text(size= 30),
        strip.background = element_rect(colour = "black")) +
  scale_fill_manual(name = "Isoform resolution \nper MIG", values = c("#0C6A2B", "#8D1D1D", "#D2701F", "#044771")) +
  facet_grid(end ~ library, scales = "free", space = "free")


