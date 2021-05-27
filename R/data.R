# data.R - bam - 8/13/18

#' RNA-Seq from the paper referenced below.
#'
#'
#' @docType data
#' @keywords datasets
#' @name mdd.RNAseq
#' @usage data(mdd.RNAseq)
#' @references
#' Trang T. Le, et al. “Identification and replication of RNA-Seq gene network modules associated with depression severity,” Translational Psychiatry. 2018.
#' @format
#' Data frame with 157 mdd and hc subjects, 5912 genes.
"mdd.RNAseq"

#' RNA-Seq from the paper referenced below.
#'
#'
#' @docType data
#' @keywords datasets
#' @name mdd.RNAseq.small
#' @usage data(mdd.RNAseq.small)
#' @format
#' Data frame with 157 mdd and hc subjects with 
#' 500 randomly selected genes from the original 5912 genes.
"mdd.RNAseq.small"

#' case-control interaction effect data
#'
#' @docType data
#' @keywords datasets
#' @name case.control.3sets
#' @usage data(case.control.3sets)
#' @format
#' Data frame with 300 samples (100 samples in train/holdout/test) and 
#' 100 features with moderate (interaction) effect size.
"case.control.3sets"

#' quantitative trait main effect data
#'
#' @docType data
#' @keywords datasets
#' @name qtrait.3sets
#' @usage data(qtrait.3sets)
#' @format
#' Data frame with 300 samples (100 samples in train/holdout/test) and 
#' 100 features with moderate (main) effect size.
"qtrait.3sets"

#' Test data from Rinbix
#'
#' @docType data
#' @keywords datasets
#' @name testdata10
#' @usage data(testdata10)
"testdata10"
