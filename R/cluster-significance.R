# mat <- readRDS("xeno.mat.logtpm.rds")
# k <- hcluster::hcutree(hcluster::hcluster(mat))
# hcsigObj <- hcsig(k=k, mat=mat, fc.value=3, p.value=10^(-4), reorder=TRUE)
# hcsigCut <- hcsig_cut(hcsigObj, n.sig.1=50, n.sig.2=10)


cluster_reorder <- function(param1, param2, decreasing=TRUE) {
  order(lengths(param1), lengths(param2), decreasing=decreasing)
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



most_significant <- function(List, by='small') {
  df <- list_to_df(List)

  if (by == 'small') {
    # sort by id and abs(value)
    ordered.df <- df[order(df$name, abs(df$val)), ]
  }
  else if (by == 'large') {
    # sort by id and reverse of abs(value)
    ordered.df <- df[order(df$name, - abs(df$val)), ]
  }

  # take the first row within each id
  return.df <- ordered.df[!duplicated(ordered.df$name), ]
  return.df.2 <- return.df[order(return.df$id), ]

  df_to_list(return.df.2, col=return.df.2$id)
}


hcsig <- function(k, mat, fc.value=3, p.value=10^(-4), reorder=TRUE, fc.sort=T, pval.sort=F) {
  sig.1 <- lapply(k, function(kk) DEgenes(k=kk, mat=mat, fc.value=fc.value, p.value=p.value, fc.sort=fc.sort, pval.sort=pval.sort))
  sig.2 <- most_significant(sig.1)
  sig.1.zeros <- sig.1[!names(sig.1) %in% names(sig.2)]
  sig.2 <- c(sig.2, sig.1.zeros)

    if (isTRUE(reorder)) {
    ord <- cluster_reorder(sig.1, sig.2)
    k <- k[ord]
    sig.1 <- sig.1[ord]
    sig.2 <- sig.2[ord]
  }

  List <- list(k=k, sig.1=sig.1, sig.2=sig.2)
  names(List) <- c("k", "sig.1", "sig.2")
  List

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


hcsig_cut <- function(obj, n.sig.1=50, n.sig.2=10, by='either'){

  obj <- input_check(obj)

  sig.1.bool <- lengths(obj$sig.1) >= n.sig.1
  sig.2.bool <- lengths(obj$sig.2) >= n.sig.2

  if (by=='either') {
    cutter <- sig.1.bool | sig.2.bool
  }
  else if (by=='both') {
    cutter <- sig.1.bool & sig.2.bool
  }

  list(k=obj$k[cutter], sig.1=obj$sig.1[cutter], sig.2=obj$sig.2[cutter])
}


hcsim_cut <- function(obj, jac.out=NULL, jac.cut=0.75) {

  obj <- input_check(obj)
  k.cut <- filter_by_jaccard(x=obj$k, out=jac.out, jac.cut=jac.cut)
  kept <- obj$k %in% k.cut

  list(k=k.cut,  sig.1=obj$sig.1[kept], sig.2=obj$sig.2[kept])
}


# hcut <- function(obj, by=c("significance", "similarity"), n.sig.1=50, n.sig.2=10, jac.cut=0.75) {
#   if (length(by) > 1) by=by[[1]]
# 
#   if (by == "significance") {
#     hcsig_cut(obj, n.sig.1=n.sig.1, n.sig.2=n.sig.2)
#   }
# 
#   if (by == "similarity") {
#     hcsim_cut(obj, jac.cut=0.75)
#   }
# }
# 
