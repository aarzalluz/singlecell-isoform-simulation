# Single-cell RNAseq for the study of isoforms: simulation code

This repository contains the custom code used to produce the simulation 
that is shown in our manuscript "Single-cell RNAseq for the study of isoforms:
how is that possible?", by Ana Conesa and Ángeles Arzalluz-Luque.

### .R files

- The prefixes `sr` and `lr` are short for short and long reads, respectively,
and designate which simulation the code is associated with.
- `simulation_function` files contain custom functions that implement some of the
main steps of the simulation.
- `simulation_pipeline` files are scripts that contain the necessary code to run
the simulation from beginning to end. It also contains plenty of information on the logic behind
the code in the form of comments.

### Data

The respository includes the necessary data to run the pipeline from the beginning,
in the form of .RData files. These files are also designated whith appropriate
prefixes where needed.

- `sr_trancsriptome.rda`: contains the transcript sequences from in Tardaguila et al. (2017),
upon which we based our simulation, and their fasta headers.
- `gene.isoform_table.rda`: contains the isoform name and the gene they belong to, also obtained
from the Tardaguila et al. (2017) transcriptome.
- `transcript_expression.rda`: contains the isoform expression table (bulk RNAseq data) produced
by Tardaguila et al. (2017), for the two samples (neural stem cells and oligodendrocytes) and 
their two replicates.
- `sr_sim.isoform.results.NSC.rda` contains the output of RSEM+STAR, and can be loaded into R to
avoid running this part of the pipeline (see next section).
- The rest of the files containing the `sr` prefix include intermediate step data. What they 
correspond to and how they can be loaded is specified in the `sr_simulation_pipeline.R` script.

### Transcriptome files

- `transcriptome.fasta` is the transcriptome fasta file from Tardaguila et al. (2017), necessary
to simulate reads from full-length transcripts using the polyester package.
- `annotation.gtf` is the annotation generated by Tardaguila et al. (2017) for this transcriptome.
Although it is not necessary to run the R code, it is required to run RSEM + STAR.

### Instructions for running the simulations

Download this repository, and change your working directory in R to the corresponding folder. Then, source the 
.R files containing the custom simulation functions and load the data. For instance, to start running the short-read simulation, execute the following in the R terminal

```
source("sr_simulation_functions.R")
load("data/sr_transcriptome.rda")
```

as specified in the frist lines of the `sr_simulation_pipeline.R` script. Then, run the rest of the code chunk-by-chunk. 
Running the entire script at once is strongly advised against, as the code contains options to run every step 
of the pipeline as well as alternatives to load the already processed data (.rda files) to skip the most computationally 
costly steps. Read the script carefully and follow the instructions in the comments, paying special attention to
occassions where it will be necessary to load .rda files to proceed to the next step.

### Running RSEM + STAR for the short read simulation

The `sr_simulation_pipeline.R` script is divided into two parts: the generation of the simulated short
reads, and the analysis of the expression matrix resulting from quantification of these reads. Although
the processed RSEM results can be loaded from .rda files as described in the `sr_simulation_pipeline.R` script, 
here we specify the commands used to replicate the short read mapping and quantification we performed:

```
# prepare the reference genome (must have the chromosome sequences of the mouse genome as fasta files)
rsem-prepare-reference data/mouseM10/chromosomes rsem_ref/mouseM10 --gtf transcriptome/annotation.gtf

# reads simulated across full-length transcripts are paired-end to ensure sufficient coverage
rsem-calculate-expression --paired-end NSC_full_length/sample_01_1.fasta simulation/NSC_full/sample_01_2.fasta rsem_ref/mouseM10 rsem_results/NSC_full --no-qualities

# example quantification of short reads generated from 100 nt fragments on the 3' end
rsem-calculate-expression NSC_100nt/sample_01.fasta rsem_ref/mouseM10 rsem_results/NSC_100nt --no-qualities

# move to rsem_results folder to begin the analysis detailed in the second part of the script
cd rsem_results
``` 

### References

Tardaguila, Manuel, Lorena de la Fuente, Cristina Marti, Cécile Pereira, Hector Risco, 
Maravillas Mellado, Marissa Macchietto, Kenneth Verheggen, and Mariola Edelmann. 2017. 
“SQANTI : Extensive Characterization of Long Read Transcript Sequences for Quality Control 
in Full-Length Transcriptome Identification and Quantification.” bioRxiv, August, 1–28. 
https://doi.org/10.1101/118083.
