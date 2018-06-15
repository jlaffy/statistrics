Args <- function(FUN, except = c('mat', 'sep', 'collapse', 'df', 'cachePath'), ...) {

  arguments <- formals(FUN)
  dots <- list(...)

  if (!is.null(except)) {
    arguments <- arguments[!names(arguments) %in% except]
  }

  if (is.null(names(dots)) & is.null(names(arguments))) {
    return(dots)
  }

  if (length(dots) > 0) {

    if (is.null(names(dots)) & !is.null(names(arguments))) {
      names(dots) <- names(arguments)
      arguments <- dots
    }

    changing <- dots[names(dots) %in% names(arguments)]
    arguments[names(changing)] <- changing
  }

  arguments
}


#' Generate one string from a list of strings
#'
#' @param args named list of arguments
#' @param argNames boolean indicating whether argument names should be included in the final output.
#' @param sep a character to separate arguments and their names. Defaults to ":".
#' @param collapse a character to separate arguments(+names) from one another. Defaults to "__".
#'
#' @return string of concatenated arguments
#' @export
#'
string <- function(args, argNames=FALSE, sep=":", collapse="_") {

  if (isTRUE(argNames)) {
    args <- lapply(1:length(args), function(i) paste(names(args)[i], args[[i]], sep=sep))
  }

  paste(args, collapse=collapse)
}


#' Generate filename to match Cache files
#'
#' @param pipeName the name of main function
#' @param fnName the name of function that is a step in the main function
#' @param args a list of arguments to the main function
#' @param cachePath a character string providing path to the Cache directory.
#' @param sep passed to \code{string}. a character to separate arguments and their names. Defaults to ":".
#' @param collapse passed to \code{string}. a character to separate arguments(+names) from one another. Defaults to "__".
#' @param argNames passed to \code{string}. Boolean indicating whether argument names should be included in the final output.
#'
#' @return character string that is the generated filename. It includes the main function and step function names, and the arguments.
#' @export
#'
makeFilename <- function(args, pipeName=NULL, fnName=NULL, cachePath=".", sep=":", collapse="__", argNames=FALSE) {

  if (stringr::str_sub(cachePath, -1) != "/") {
    cachePath <- paste(cachePath, "/", sep="")
  }

  components <- c(pipeName, fnName, args)
  String <- string(args=components, argNames=argNames, sep=sep, collapse=collapse)

  print(paste(cachePath, String, ".rds", sep=""))
  paste(cachePath, String, ".rds", sep="")
}


#' Combine data from files
#'
#' @param path a character vector of full path names; the default corresponds to the working directory, \code{getwd()}.
#' @param pattern an optional regular expression. Only file names which match the regular expression will be returned.
#'
#' @return Combined files' input - each file's contents contained in one list.
#' @export
#'
combine <- function(path, pattern) {
  f <- list.files(path = path, pattern = pattern, full.names = T)
  out <- sapply(f, function(ff) readRDS(ff), USE.NAMES = F, simplify = FALSE)
  unlist(out, recursive = FALSE)
}
