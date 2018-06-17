cluster_reorder <- function(param1, param2, param3, decreasing=TRUE) {
  order(lengths(param1), lengths(param2), lengths(param3), decreasing=decreasing)
}


list_to_df <- function(List) {
  dfs <- lapply(1:length(List), function(i) {
    name <- names(List[[i]])
    val <- List[[i]]
    id <- rep(names(List[i]), length(List[[i]]))
    data.frame(name, val, id) })
  do.call(rbind, dfs)
}


df_to_list <- function(df, col=NULL) {
  if (is.null(col)) col <- df$id

  sapply(split(df, df$id), function(s) {
    v <- s$val
    names(v) <- s$name
    v })
}


input_check <- function(obj, Names=NULL) {

  if (!all(sapply(obj[1:3], class) == "list")) {
    stop("Objects within <obj> must be lists.
         One or more objects within <obj> are not of class list.")
  }

  if (length(obj) < 3) {
    stop(paste("<obj> must comprise three lists:
               1. k (Clusters' cell members)
               2. sig.1 (Clusters' significant gene members)
               3. sig.2 (Clusters' most significant gene members).
               Only" , length(obj), "lists were found."))
  }


  if (is.null(Names) | length(Names) < 3) {
    Names <- c("k", "sig.1", "sig.2")
  }

  if (!all(Names %in% names(obj))) {
    names(obj)[1:3] <- Names[1:3]
  }

  obj

}


#' Apply function to every pair between x and a list
#'
#' @param a value or vector
#' @param vect list of values or vectors of the same type as 'a'
#' @param FUN function to be applied on every pair. eg `sum(a, vect[[1]])` if a and vect are numeric
#'
#' @return result of function call for each pair
#' @export
#'
pairs_apply <- function(a, vect, FUN=overlap) {
  sapply(vect, function(b) FUN(a=a, b=b))
}


#' Apply function to every pair in a list
#'
#' @param vect a list. All elements must be of the same type.
#' @param FUN function to be applied on every pair. eg sum(a, vect[[1]]) if a and vect are numeric
#'
#' @return result of function call for each pair
#' @export
#'
pairs_apply_all <- function(vect, FUN=overlap) {

  overlaps <- NULL

  for (i in 1:length(vect)) {
    x <- vect[[i]]
    overlap.scores <- list(pairs_apply(a=x,
                                       vect=vect,
                                       FUN=FUN))
    overlaps <- c(overlaps, overlap.scores)
  }

  overlaps
}
