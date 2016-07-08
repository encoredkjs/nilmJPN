grepl.pattern <- function(united_state, targetPattern){
  
  n <- length(targetPattern)
  targetPattern.boolean <- sapply( 0:(n-1), function(i) shift( united_state == targetPattern[i+1],-i))
  return( which( rowSums( targetPattern.boolean ) == n ) )
}
