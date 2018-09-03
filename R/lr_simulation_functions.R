#' Simulate sequencing of several cells per SMRT cell
#'
#' Reduces the expression values provided (1 cell or 1 bulk experiment) to simulate a scenario where several single cells are sequenced
#' in a SMRT cell, and the subsequent loss of depth per cell.
#'
#' @param expression A vector of expression values for a sample.
#' @param reduction An integer providing the number of cells per SMRT cell for which expression is simulated.
#'
#' @return A vector containing expression values for one cell in the simulated sequencing scenario.
reduce_expression <- function(expression, reduction){
  # divide TPM values by selected rate
  red_expression <- expression/reduction
  # set all expression values that are below 1 after division to zero
  red_expression[which(red_expression < 1)] <- 0
  return(red_expression)
}

#' Analysis of podium changes between samples
#' 
# This file contains the set of custom functions necessary to run the code in the sr_simulation_pipeline script.

#' Performs a gene-by-gene comparison of the most expressed isoform in two samples to detect podium changes (i.e. switches in 
#' the most highly expressed isoform).
#'
#' @param expr_NSC, expr_OLD Lists containing the isoform expression matrix for every simulation (one list per sample).
#' @param mig_OLD, mig_NSC Lists containing the number of isoforms detected per MIG in every simulation (one list per sample).
#' @param sim Integer indicating the list slot (i.e. simulation) for which the comparison will be performed.
#'
#' @return A vector containing, as character, the names of the genes that switch their most expressed isoform between samples.
compare_mig <- function(expr_NSC, expr_OLD, mig_OLD, mig_NSC, sim){
  # dependencies
  require(dplyr)
  
  # table for NSC
  mig_expr.NSC <- expr_NSC[[sim]]
  mig_expr.NSC <- mig_expr.NSC[mig_expr.NSC$associatedGene %in% mig_NSC[[sim]]$values, ]  # keep MIG only
  top_iso.NSC <- mig_expr.NSC %>% group_by(associatedGene) %>% filter(expression==max(expression))  # keep only the most expressed isoform
  # table for OLD
  mig_expr.OLD <- expr_OLD[[sim]]
  mig_expr.OLD <- mig_expr.OLD[mig_expr.OLD$associatedGene %in% mig_OLD[[sim]]$values, ] # keep MIG only
  top_iso.OLD <- mig_expr.OLD %>% group_by(associatedGene) %>% filter(expression==max(expression))  # keep only the most expressed isoform
  
  # MOST EXPRESSED ISOFORM COMPARISON
  common_mig <- intersect(top_iso.NSC$associatedGene, top_iso.OLD$associatedGene) # make list of mig common to both samples
  top_iso.NSC <- top_iso.NSC %>% filter(associatedGene %in% common_mig) # select only the genes that are common between samples
  top_iso.OLD <- top_iso.OLD %>% filter(associatedGene %in% common_mig)
  # NOTE: there might be genes that co-express two isoforms in OLD or NSC, but have a single most expressed isoform in the other
  dup <- top_iso.OLD[which(duplicated(top_iso.OLD$associatedGene)), ]$associatedGene # select the names of the duplicated genes
  top_iso.OLD <- top_iso.OLD %>% filter(!(associatedGene %in% dup)) # remove the genes that follow co-expression patterns
  top_iso.NSC <- top_iso.NSC %>% filter(!(associatedGene %in% dup)) # in both tables
  
  # compare the most expressed isoform in both samples, and select indices of genes with a podium change
  changes <- which(top_iso.NSC$isoform != top_iso.OLD$isoform)
  # select genes that suffer podium changes
  podium_ch <- top_iso.NSC[changes, ]$associatedGene
  # return the genes that suffer podium changes
  return(podium_ch)
}