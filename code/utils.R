

## IQR function
inter_quantile <- function(x, probs = c(0.25, 0.5, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}


## standardizing function (Sacha qgraph)
scale2 <- function(x) {
  if (all(is.na(x))) return(NA)
  if (sd(x,na.rm=TRUE)!=0){
    return((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
  } else {
    return(rep(0, length(x)))
  }
}

## define the node names used throughout
nodenames <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")
