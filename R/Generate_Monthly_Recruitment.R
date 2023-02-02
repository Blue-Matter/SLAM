#' Title
#'
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
Generate_Monthly_Recruitment <- function(mu=4, sigma=2) {
  n <- length(mu)
  rec <- list()
  for (i in 1:n) {
    tt <- dnorm(0:11, 5.5, sigma[i])
    tt <- tt/sum(tt)
    tt <- rep(tt,5)
    rec_pattern <- rep(0, 12)
    rec_pattern[floor(mu[i])] <- tt[6]
    j <- floor(mu[i])+1
    k <- 1

    while(j<=12) {
      rec_pattern[j] <- tt[6+k]
      j <- j+1
      k <- k+1
    }
    for (j in 1:(floor(mu[i])-1)) {
      rec_pattern[j] <- tt[6+k]
      k <- k+1
    }
    rec[[i]] <- rec_pattern
  }
  rec <- do.call('rbind', rec)
  apply(rec, 2, sum)/sum(apply(rec, 2, sum))
}
