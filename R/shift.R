shift <- function(vec,n) {
  if (n == 0){
    vec
  }else if (n > 0){
    c(rep(NA, n), head(vec, -n))
  }else if (n < 0){
    if( length(tail(vec,n)) == 0 ){
      return(rep(NA,length(vec)))
    }
    c(tail(vec, n), rep(NA, -n))
  }
}