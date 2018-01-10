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

- `sr_transriptome.rda`: contains the transcript sequences from in Tardaguila et al. (2017),
upon which we based our simulation, and their fasta headers.
- `transcriptome.fasta` is the transcriptome fasta file from Tardaguila et al. (2017), necessary
to simulate reads using the polyester package.
- `gene.isoform_table.rda`: contains the isoform name and the gene they belong to, also obtained
from the Tardaguila et al. (2017) transcriptome.
- `transcript_expression.rda`: contains the isoform expression table (bulk RNAseq data) produced
by Tardaguila et al. (2017), for the two samples (neural stem cells and oligodendrocytes) and 
their two replicates.
- `sr_sim.isoform.results.NSC.rda` contains the output of RSEM+STAR, and can be loaded into R to
avoid running this part of the pipeline.
- The rest of the files containing the `sr` prefix include intermediate step data. What they 
correspond to is specified in the `sr_simulation_pipeline.R` script.

### References

Tardaguila, Manuel, Lorena de la Fuente, Cristina Marti, Cécile Pereira, Hector Risco, 
Maravillas Mellado, Marissa Macchietto, Kenneth Verheggen, and Mariola Edelmann. 2017. 
“SQANTI : Extensive Characterization of Long Read Transcript Sequences for Quality Control 
in Full-Length Transcriptome Identification and Quantification.” bioRxiv, August, 1–28. 
https://doi.org/10.1101/118083.
