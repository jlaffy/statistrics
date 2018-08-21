#!/usr/bin/env Rscript

fun1 <- function(mat) {
  geme <- statistrics::logTPM(sort(rowMeans(statistrics::TPM(mat))))
  DF <- data.frame(gene = names(geme),
                   val = geme,
                   lev = gl(31, 276))
  saveRDS(DF, file = 'DF.gene.bins.rds')
  DF
}


fun2 <- function(DF) {
  Groups <- lapply(split(DF$gene, DF$lev), as.character)
  saveRDS(Groups, 'bins.rds')
  Groups
}


fun3 <- function(mat, Groups) {
  DF <- apply(mat, 2, function(col) {
                sapply(Groups, function(Group) {
                         rep(mean(col[Group]), 276) })
             })
  genes <- unlist(Groups, use.names = F)
  DF.ord <- DF[order(genes), ]
  rownames(DF.ord) <- sort(genes)
  saveRDS(DF.ord, 'malignant.logtpm.binned_controls.rds')
  DF.ord
}


bin_indices <- function(v) {
  DF <- readRDS('/Volumes/tirosh/sbjulie/GBM_2/binning/DF.gene.bins.rds')
  sapply(v, function(vi) {
		   DF %>%
			 dplyr::filter(gene %in% vi) %>%
			 dplyr::select(lev) %>%
			 unlist(use.names = F) %>%
			 as.integer },
		   simplify = FALSE,
		   USE.NAMES = TRUE)
}


.bin_sample <- function(x, n = 100, replace = FALSE) {
  unlist(sapply(x, sample, n, replace = replace,
				simplify = FALSE,
				USE.NAMES = TRUE),
		 use.names = F)

}

bin_sample <- function(indices, n = 100, replace = FALSE) {
  bins <- readRDS('/Volumes/tirosh/sbjulie/GBM_2/binning/bins.rds')

  if (any(unlist(lengths(bins) > n))) {
	replace = TRUE
  }

#   sapply(indices, function(i) {
# 		   .bin_sample(x = bins[i], n = n, replace = replace)},
# 		 simplify = FALSE,
# 		 USE.NAMES = TRUE)

  sapply(indices, function(i) {
		   .bin_sample(x = bins[i], n = n, replace = replace)},
		 USE.NAMES = TRUE)
}


control <- function(v, n = 100, replace = FALSE) {
  indices <- bin_indices(v = v)
  bin_sample(indices = indices, n = n, replace = replace)
}
