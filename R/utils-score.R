#' Top Genes across Programs
#'
#' @param Programs list of character vectors. Each vector is a program
#' @param cutoff a numeric value indicating how many top genes to return.
#'
#' @return table of top genes and their counts across programs
#' @export
#'
top_program_genes <- function(Programs, cutoff=50) {

  sort(table(unlist(Programs)), decreasing=T)[1:cutoff]

}


#' Programs
#'
#' @param List list of named numeric vectors. Numbers are p-values and names are genes.
#' @param cutoff a numeric value indicating the cutoff for program size

#' @return list of character vectors. Each vector is a program
#' @export
#'
program <- function(List, cutoff=50) {

  if (is.null(cutoff)) {
	Programs <- lapply(List, function(genes) names(genes))
  }

  Programs <- lapply(List, function(genes) {
					   Genes <- names(genes)[1:cutoff]
					   Genes[!is.na(Genes)] })

  names(Programs) <- paste("P", 1:length(Programs), sep="")

  cat(paste("\nTop", cutoff, "program genes:\n"))
  print(top_program_genes(Programs, cutoff=50))

  Programs

}


#' Program Scores
#'
#' @param mat a matrix of gene expression data (cells by genes)
#' @param programs either a single or a list of character vectors (genes). Each vector is a program.
#' @param center if TRUE, the resulting score matrix is centered.
#' @param center.rowWise if TRUE, centering will be performed row-wise. Else, centering by column.
#'
#' @return list (if many=F) or matrix (if many=T) of program scores
#' @export
#'
score <- function(mat, programs, center = FALSE, center.rowWise=TRUE, na.rm = FALSE) {

  result <- sapply(programs, function(program) {colMeans(mat[program, ], na.rm = na.rm)},
#				   USE.NAMES=T,
				   simplify=F)

  if (length(result) > 1) {
	result <- do.call(cbind, result)
  }

  if (is.null(colnames(result))) {
	colnames(result) <- paste("P", 1:ncol(result), sep="")
  }

  if (is.null(rownames(result))) {
	rownames(result) <- colnames(mat)
  }

  if (isTRUE(center)) {
	result <- center(result, rowWise = center.rowWise)
  }

  result

}
