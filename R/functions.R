assignGroups <- function(df, var="mean") {
  group <- 1
  df$groups <- NA
  for (i in 1:nrow(df)) {
    if (i == 1) {
      if (!is.na(df[[var]][i])) {
        df$groups[i] <- group + 1
      }
    } else {
      if (!is.na(df[[var]][i])) {
        df$groups[i] <- group
      } else if (is.na(df[[var]][i]) & !is.na(df[[var]][i-1])) {
        group <- group+1
      }
    }
  }
  df$x <- 1:nrow(df)
  df
}

dropLastNAs <- function(df, var='mean') {
  ind <- df[[var]] %>% is.na() %>% rev()
  ind <- min(which(!ind))
  head(df, -ind+1)
}

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

h2CR <- function(h) {
  (4*h)/(1-h)
}

CR2h <- function(CR) {
  CR/(CR+4)
}


