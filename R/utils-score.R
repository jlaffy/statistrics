top_program_genes <- function(Programs, cutoff=50) {

  sort(table(unlist(Programs)), decreasing=T)[1:cutoff]

}


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
