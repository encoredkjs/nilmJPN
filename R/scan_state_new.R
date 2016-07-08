scan_state_new <- function(speed_log){
  
  idx.nneg.speed <- which( speed_log >= 0 ) # index of non-negative speed
  idx.neg.speed  <- which( speed_log <  0 ) # index of     negative speed
  
  acceleration   <- c(NA, diff( speed_log ))
  idx.pos.accel  <- which( acceleration >   0 ) # index of     positive acceleration 
  idx.npos.accel <- which( acceleration <=  0 ) # index of non-positive acceleration
  
  idx.state1 <- intersect( idx.nneg.speed, idx.pos.accel  ) 
  idx.state2 <- intersect( idx.nneg.speed, idx.npos.accel ) 
  idx.state3 <- intersect( idx.neg.speed,  idx.npos.accel ) 
  idx.state4 <- intersect( idx.neg.speed,  idx.pos.accel  ) 
  
  state <- rep(0, length(speed_log))
  state[idx.state1] <- 1
  state[idx.state2] <- 2
  state[idx.state3] <- 3
  state[idx.state4] <- 4
  
  return(state)
}
