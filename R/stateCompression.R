stateCompression <- function( state ){
  
  stateChangePos <- which( head(state,-1) != tail(state,-1) )
  
  state.str <- c(0,stateChangePos) + 1
  state.end <- c(stateChangePos,length(state))
  state.val <- state[ state.str ]
  
  return( data.frame( str = state.str, end = state.end, val = state.val) )
}
