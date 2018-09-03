# Long read simulations

This vignette contains the necessary code to reproduce the long read simulations in our manuscript "Single-cell RNAseq for the
study of isoforms -how is that possible?", which can be found [here](https://doi.org/10.1186/s13059-018-1496-z).

### Prepare environment

Load the required libraries:

```
library(readr)
library(dplyr)
library(ggplot2)
```

Load the custom functions in the `lr_simulation_functions.R` file:

`source("lr_simulation_functions.R")`

Load isoform-to-gene table and Tardaguila et al. (2018) transcript expression values:

```
# load isoform-to-gene table (data.frame)
load("data/gene.isoform_table.rda") 

# load neural stem cell and oligodendrocyte expression values for two replicates (data.frame)
load("data/transcript_expression.rda") 

# select one replicate per sample for the simulation
expr <- data.frame(transcript = quantification$PacBio, NSC = quantification$POST_NSC1, OLD = quantification$POST_OLD1) 
```

### Simulate multiplexed long-read sequencing results

In order to explore the number of cells vs sequencing depth trade-off impacting long read sequencing in single-cell RNAseq, we simulate a situation where a growing number of cells is sequenced in a Sequel flow cell with a total throughput of 1 million long reads. This is done by downsampling the previously loaded expression values, which are scaled in TPM units, using the `reduce_expression()` custom function:

```
# create a range of cell numbers for the simulation experiment
cells <- c(2,6,10,16,20)

# using these vectors, perform expression value downsampling for both samples
reduced_expression.NSC <- lapply(cells, reduce_expression, expression = expr$NSC) 
reduced_expression.OLD <- lapply(cells, reduce_expression, expression = expr$OLD) 

# aproximate no. of reads per cell in each simulation
reads <- as.integer(round(1000000 / cells)) 
```

In order to assess the detection potential of the technology as a growing number of cells is sequenced, we tidy and filter the simulation results as follows:

``` 
expr_matrices.NSC <- list()
expr_matrices.OLD <- list()

# use list of reduced expression to create expression matrices (i.e. number of times each transript is detected, as per TPM scaling) for each simulation
for (i in seq(length(cells))){
  expr_matrices.OLD[[i]] <- data.frame(gene.isoform_table, expression = reduced_expression.OLD[[i]]) %>% filter(expression!=0)
  expr_matrices.NSC[[i]] <- data.frame(gene.isoform_table, expression = reduced_expression.NSC[[i]]) %>% filter(expression!=0)
}
```

### Analyse and plot simulation results

First, we count the number of genes for which more than one transcript isoform is detected (using the `rle()` function), based on the transcript expression results simulated. Higher multiplexing will gradually lose the less abundant transcripts, and therefore the number of genes detected as multi-isoform. This is done for both samples.

```
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
```

Then, we plot the results and compare the number of genes for which more than isoform is detected to the 2423 multi-isoform genes detected by Tardaguila et al. (2018) (see Figure 5.c in our manuscript).

```
# total no. of genes with >1 isoform in Tardaguila et al. 2018 data
mig.bulk <- 2423 

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
```

The second simulation analysis is dedicated to the amount of isoform switches observed between NSC and oligodendrocytes, and how many are observed in higher multiplexing situations. To evaluate this, we used the `compare_mig()` custom function, which performs a gene-by-gene comparison of the most expressed isoform per gene in the two samples, and returns the genes for which a switch is produced. 

```
# make empty list
iso_switch_genes <- list() 

# call compare_mig() for each simulation, and store isoform switching genes in each case
for(i in seq(length(expr_matrices.NSC))){
  iso_switch_genes[[i]] <- compare_mig(expr_matrices.NSC, expr_matrices.OLD, counts.NSC, counts.OLD, i)
}

# count no. of podium changes per simulation
iso_switch_counts <- lapply(iso_switch_genes, length) %>% unlist

# make data frame for plotting
pch.summary <- data.frame(cell.no = cells %>% as.factor, isoform_switches = iso_switch_counts)
```

In this case, we plot the number of isoform switches detected in each multiplexing scenario, and compare it to the 337 podium changes present in the original Tardaguila et al. (2018) data.

```
# number of isoform switches in Tardaguila et al. 2018 data
pch.bulk <- 337  

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
```


