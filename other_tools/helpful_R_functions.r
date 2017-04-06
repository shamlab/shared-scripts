#' Count pairwise observations
#' @param dat a data.frame object
#' @return a matrix with the number of pairwise observations 
count.pairwise <- function(dat) {
	x <- !is.na(dat)
	return(crossprod(x))
}
