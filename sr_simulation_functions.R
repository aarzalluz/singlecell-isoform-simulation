# This file contains the set of custom functions necessary to run the code in the sr_simulation_pipeline script.

#' Reformat fasta file into character vector
#' 
#' Creates a character vector containing headers in the odd indices, and their nucleotide sequence as one string in the even indices.
#'
#' @param file_path A string indicating the path to the target fasta file.
#' @param output_path A string indicating the path where output, an unformatted sequence file, will be written.
#'
#' @return No object is returned, but a file is created in the designated output path.
undo_fasta <- function(file_path, output_path){
  # dependencies
  library(readr)
  library(stringr)
  
  # LOAD TRANSCRIPTS AS CHARACTER VECTORS
  fasta <- read_lines(file_path)
  
  # UNDO FASTA FORMAT: ONE STRING PER TRANSCRIPT SEQUENCE
  headers_logical <- str_detect(fasta, ">") # positions with TRUE will contain the fasta headers
  headers <- fasta[headers_logical == TRUE] # returns only the fasta headers
  fasta_lines <- which(headers_logical) # returns indices containing fasta headers
  
  # create empty character vector to store all the sequences (one per header)
  fullseqs <- character(length(fasta_lines))
  
  # merge fasta lines into the same string and store in character vector
  for(i in seq_len(length(fasta_lines))){
    # make condition for last element
    if(i == length(fasta_lines)){
      idx1 <- fasta_lines[i] + 1
      idx2 <- length(fasta) # indices beyond the last element in fasta_lines are not accessible by sum
      fullseqs[i] <- paste(fasta[idx1:idx2], collapse = '')
      # unless it is the last iteration, do this:
    } else {
      idx1 <- fasta_lines[i] + 1  # first index containing a sequence fragment
      idx2 <- fasta_lines[i+1] -1 # last index containing a sequence fragment
      fullseqs[i] <- paste(fasta[idx1:idx2], collapse = '') # make a single seq string
    }
  }
  
  # result: fullseqs contains one transcript per position, headers contains the corresponding headers in order
  
  # MERGE HEADERS AND SEQS INTO OBJECT
  unformatted_tr <- character(length(fullseqs)+length(headers)) # create empty vector to fill with IDs and subset transcripts
  idx_id <- seq(1, length(unformatted_tr)-1, 2) # get uneven indices
  idx_seqs <- seq(2, length(unformatted_tr), 2) # get even indices
  unformatted_tr[idx_id] <- headers # add IDs to uneven indices
  unformatted_tr[idx_seqs] <- fullseqs # add full transcripts to even indices
  
  # WRITE UNFORMATTED TRANSCRIPTS TO FILE
  write_lines(unformatted_tr, output_path)
}

#' Create transcript fragments
#' 
#' Nucleotide sequences from a set of transcripts are subset, starting at the 3' or 5' end and ending at the indicated fragment size.
#'
#' @param tr_set A character vector containing fasta headers in the odd indices, and transcript sequences in the even. If not available,
#' a regular fasta file can be reformatted using the \code{"undo_fasta} function.
#' @param fragment_length An integer indicating the selected length of the transcript fragments.
#' @param filename A string. Path and file name where the fragmented transcripts will be output. Do not include the ".fasta" or ".fa"
#' extension, as it is provided by the function to avoid extension errors.
#' @param select_end A string indicating the transcript end in which the fragments will be generated.
#'
#' @return No object is returned. Instead, a file containing the fragmented transcripts and their headers is output to the specified 
#' location.
subset_tr <- function(tr_set, fragment_length, filename, select_end = c("3end", "5end")){
  # dependencies
  require(stringr)
  
  # SEPARATE HEADERS AND SEQUENCES
  idx_id <- seq(1, length(tr_set)-1, 2) # get uneven indices (headers)
  idx_seqs <- seq(2, length(tr_set), 2) # get even indices (transcripts)
  seqs <- tr_set[idx_seqs] # get transcripts
  headers <- tr_set[idx_id] # get headers
  
  # MAKE SUBSTRINGS AND MERGE SEQS AND HEADERS IN FASTA FILE
  # subset from selected end
  select_end <- match.arg(select_end)
  
  if(select_end == "3end"){
    subseqs <- str_sub(seqs, start = -fragment_length, end = -1) # select the last x nucleotides from the transcripts in fullseqs
  } 
  else if (select_end == "5end"){
    subseqs <- str_sub(seqs, start = 1, end = fragment_length) # select the first x nucleotides from the transcripts in fullseqs
  }
  
  result <- character(length(tr_set)) # create empty vector to fill with IDs and subset transcripts
  result[idx_id] <- headers # add IDs to uneven indices
  result[idx_seqs] <- subseqs # add subsequences to even indices
  
  write_lines(result, paste0(filename, ".fa"))
}

#' Load and process isoform expression results
#'
#' The isoform expression file output by RSEM is loaded and the number of isoforms per gene counted.
#'
#' @param results_file A string indicating the path and filename of the tab-delimited file output by RSEM containing isoform expression.
#'
#' @return A hash object (class defined in the hash library) containing gene IDs as keys and number of isoforms expressed by the
#' gene as values.
process_sim_results <- function(results_file){
  # dependencies
  require(readr)
  require(hash)
  # load isoform results and filter out those with zero expression value
  result <- read_tsv(results_file) %>% filter(TPM!=0)
  # count isoforms per gene
  count <- rle(result$gene_id)
  # gene_id sorted data frame
  processed_result <- hash(count$values, count$lengths)
  
  return(processed_result)
}

#' Edit a hash object to preserve multi-isoform genes only
#'
#' Removes mono-isoform genes from a hash object. The resulting hash is the reference (i.e. true/maximum number of isoforms per
#' multi-isoform gene that can be detected in the biological context under study) of the simulation analysis.
#' 
#' @param editable_hash Hash object containing the isoform expression results obtained from full-length simulated data quantification
#' with RSEM.
#'
#' @return A hash object (class defined in the hash library) containing the name of each multi-isoform gene, 
#' and the number of isoforms that are expressed by it.
make_mig_ref <- function(editable_hash){
  # dependencies
  require(hash)
  
  for (key in keys(editable_hash)){
    # genes with only one isoform are removed
    if(editable_hash[[key]] == 1){
      del(key, editable_hash)
    }
  }
  return(editable_hash)
}

#' Calculate percentage of resolution per gene
#'
#' Calculates the percentage of isoforms that are detected per multi-isoform gene (MIG) compared to a second situation. This 
#' "reference" is assumed to be a "true" number of isoforms per gene.
#'
#' @param iso_hash 
#' @param mig_hash 
#'
#' @return A numeric vector containing a percentage of resolution (case no. of isoforms divided by reference no. of isoforms) 
#' per MIG.
calc_MIG_percent <- function(iso_hash, mig_hash){
  # dependencies
  require(hash)
  # retrieve keys common to both hashes
  common_keys <- intersect(keys(iso_hash), keys(mig_hash))
  # make reference hash containing only common MIG
  common_h <- mig_hash[common_keys]
  # use modified reference hash to keep only MIG in simulation results
  iso_hash_mig <- iso_hash[common_keys]
  # calculate resolution percentage for all MIG
  res_percentages <- values(iso_hash_mig, USE.NAMES = F)*100 / values(common_h, USE.NAMES = F) 
  return(res_percentages)
}

#' Explore the distribution of resolution percentages within intervals
#' 
#' Calculates the percentage of genes that fall within an interval of resolution (for example, percentage of genes resolved between
#' 0 and 25%).
#'
#' @param percentages Data frame containing the percentage of resolution for each gene in the rows, and in the first column.
#' @param start Integer between 0 and 100 indicating the beginning of the resolution percentage interval.
#' @param end Integer between 0 and 100 (and > start) indicating the end of the resolution percentage interval.
#'
#' @return A single value (double) indicating the percentage of genes that fall within the designated interval of resolution.
sub_percent <- function(percentages, start, end){
  hits <- intersect(which(percentages[,1] > start), which(percentages[,1] <= end))
  sub_pcnt <- length(hits)*100 / nrow(percentages)
  return(sub_pcnt)
}
