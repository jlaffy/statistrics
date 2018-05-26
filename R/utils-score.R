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

  cat(paste("\nTop", cutoff, "program genes:\n"))
  print(top_program_genes(Programs, cutoff=50))

  Programs

}


#' Program Scores
#'
#' @param mat a matrix of gene expression data (cells by genes)
#' @param programs either a single or a list of character vectors (genes). Each vector is a program.
#' @param many boolean. If TRUE, <programs> arg is a list of vectors, each of which scores will be calculated for.
#' @param center if TRUE, the resulting score matrix is centered.
#'
#' @return list (if many=F) or matrix (if many=T) of program scores
#' @export
#'
score <- function(mat, programs, many=TRUE, center=TRUE) {

  .score <- function(mat, program) {
    sapply(1:ncol(mat), function(cell) mean(mat[program,cell], na.rm=TRUE))
  }

  if (isTRUE(many)) {
    result <- lapply(programs, function(program) .score(mat=mat, program=program))
  }


  else if (!isTRUE(many) | class(programs) == 'character') {
    result <- .score(mat=mat, program=programs)
  }

  if (isTRUE(center)) result <- center(result)

  result

}
