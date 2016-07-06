
# ENVELOPEDETECTOR This function detects upper and lower envelope of active 
# power usage sequence. 
#  
# ** Output variables
#  data$LowerEnvelope: Lower envelope
#  data$UpperEnvelope: Upper envelope
#  diffApLowerEnvelope: difference between dataAp and dataApLowerEnvelope
#
# TODO: Make decisions based on time interval, not number of valid samples
# to the left and rite.
envelopeDetector <- function( data
                              , SIDE_LENGTH_TO_CONSIDER_LOWER = 15 # Number of samples to consider for lower enevelope decision
                              , SIDE_LENGTH_TO_CONSIDER_UPPER = 15 # Number of samples to consider for upper enevelope decision
                              , use.active.power = TRUE
                              , DEBUG_MODE = FALSE  # true to see figures
                              ){ 

  if( use.active.power ) power <- data$active_power else power <- data$reactive_power 
  envelope <- envelopeDetector_sub( power, SIDE_LENGTH_TO_CONSIDER_LOWER, SIDE_LENGTH_TO_CONSIDER_UPPER )
  data     <- cbind( data, data.frame(envelope) ) 
  
  # Plot if in DEBUG_MODE
  if( DEBUG_MODE ){
    
    fig.data <- data.frame( 'timestamp' = data$timestamp, 'original' = power )
    fig.data <- cbind( fig.data, envelope )
    fig.data.melt <- reshape::melt( fig.data, id.vars='timestamp' )
      
    require(ggplot2)
    fig <- ggplot( fig.data.melt, aes(x=timestamp,y=value,group=variable,colour=variable) ) + 
      geom_line() + facet_grid(variable~., scales='free') + theme(legend.position='NONE') 
    print(fig)
  }
  return(data)
}

envelopeDetector_sub <- function( power
                                  , SIDE_LENGTH_TO_CONSIDER_LOWER = 15 # Number of samples to consider for lower enevelope decision
                                  , SIDE_LENGTH_TO_CONSIDER_UPPER = 15 # Number of samples to consider for upper enevelope decision
                                  ){ 
  
  # Use 2-iteration technique to address 'cannot touch both sides' problem
  NUM_ITERATION <- 2
  if( length(power) <= SIDE_LENGTH_TO_CONSIDER_LOWER + SIDE_LENGTH_TO_CONSIDER_UPPER ) return(power)
  
  LowerEnvelope <- power
  UpperEnvelope <- power
  
  ind.cent <- 1:length(power)
  ind.left <- sapply(SIDE_LENGTH_TO_CONSIDER_LOWER:1, function(i) c(tail(ind.cent,-i),rep(NA,i)))
  ind.rite <- sapply(1:SIDE_LENGTH_TO_CONSIDER_UPPER, function(i) c(rep(NA,i),head(ind.cent,-i)))
  ind.left <- cbind( ind.left, ind.cent )  
  ind.rite <- cbind( ind.cent, ind.rite )
  
  # The iterations
  for( t in 1:NUM_ITERATION ){
    LowerEnvelope <- pmax( apply( matrix( LowerEnvelope[ind.left], nrow=nrow(ind.left) ), 1, min ),
                           apply( matrix( LowerEnvelope[ind.rite], nrow=nrow(ind.rite) ), 1, min ) )
    LowerEnvelope[is.na(LowerEnvelope)] <- power[is.na(LowerEnvelope)]
    
    UpperEnvelope <- pmin( apply( matrix( UpperEnvelope[ind.left], nrow=nrow(ind.left) ), 1, max ),
                           apply( matrix( UpperEnvelope[ind.rite], nrow=nrow(ind.rite) ), 1, max ) )
    UpperEnvelope[is.na(UpperEnvelope)] <- power[is.na(UpperEnvelope)]
  }
  
  diff.LowerEnvelope <- power - LowerEnvelope
  diff.UpperEnvelope <- power - UpperEnvelope
  
  # Take care of head and tail parts
  LowerEnvelope[is.na(rowSums(ind.left))] <- min(power[is.na(rowSums(ind.left))])
  LowerEnvelope[is.na(rowSums(ind.rite))] <- min(power[is.na(rowSums(ind.rite))])
  UpperEnvelope[is.na(rowSums(ind.left))] <- max(power[is.na(rowSums(ind.left))])
  UpperEnvelope[is.na(rowSums(ind.rite))] <- max(power[is.na(rowSums(ind.rite))])
  
  Envelope <- matrix( c(LowerEnvelope, UpperEnvelope, diff.LowerEnvelope, diff.UpperEnvelope), ncol=4)
  colnames(Envelope) <- c('LowerEnvelope', 'UpperEnvelope', 'diff.LowerEnvelope', 'diff.UpperEnvelope')
  return(Envelope)
}


envelopeDetector.old <- function( data
                              , SIDE_LENGTH_TO_CONSIDER_LOWER = 15 # Number of samples to consider for lower enevelope decision
                              , SIDE_LENGTH_TO_CONSIDER_UPPER = 15 # Number of samples to consider for upper enevelope decision
                              , DEBUG_MODE = FALSE  # true to see figures
                              , use.active.power = TRUE
){ 
  
  require(plyr)
  # Use 2-iteration technique to address 'cannot touch both sides' problem
  NUM_ITERATION <- 2
  
  # Envelope sequences to return
  if( use.active.power ){
    data$LowerEnvelope <- data$active_power
    data$UpperEnvelope <- data$active_power
  }else{
    data$LowerEnvelope <- data$reactive_power
    data$UpperEnvelope <- data$reactive_power
  }
  
  lengthData <- nrow(data)
  
  # The iterations
  for( t in 1:NUM_ITERATION ){
    
    index <- data.frame( cent = 1:lengthData )
    index$left <- index$cent - SIDE_LENGTH_TO_CONSIDER_LOWER
    index$rite <- index$cent + SIDE_LENGTH_TO_CONSIDER_LOWER
    interior <- which((index$left >= 1) & (index$rite <= lengthData))
    if( length(interior) > 0 ){
      index <- index[ interior, ]
      
      #data$LowerEnvelope[interior] <- 
      #  daply( index, 1, function(row){ max( min(data$LowerEnvelope[row$left:row$cent]),
      #                                       min(data$LowerEnvelope[row$cent:row$rite]) ) })
      #data$UpperEnvelope[interior] <- 
      #  daply( index, 1, function(row){ min( max(data$UpperEnvelope[row$left:row$cent]),
      #                                       max(data$UpperEnvelope[row$cent:row$rite]) ) })
      
      data$LowerEnvelope[interior] <- 
        sapply( 1:nrow(index), function(x) max( min(data$LowerEnvelope[index$left[x]:index$cent[x]]),
                                                min(data$LowerEnvelope[index$cent[x]:index$rite[x]]) ))
      data$UpperEnvelope[interior] <- 
        sapply( 1:nrow(index), function(x) min( max(data$UpperEnvelope[index$left[x]:index$cent[x]]),
                                                max(data$UpperEnvelope[index$cent[x]:index$rite[x]]) ))    
    }
  }
  
  # Calculate diff
  if( use.active.power )
    data$diff.LowerEnvelope <- data$active_power   - data$LowerEnvelope
  else
    data$diff.LowerEnvelope <- data$reactive_power - data$LowerEnvelope
  
  # Take care of head and tail parts
  data$LowerEnvelope[1:min(nrow(data),SIDE_LENGTH_TO_CONSIDER_LOWER)] <- 
    min(data$LowerEnvelope[1:min(nrow(data),SIDE_LENGTH_TO_CONSIDER_LOWER)])
  data$LowerEnvelope[max(1,lengthData-SIDE_LENGTH_TO_CONSIDER_LOWER):lengthData] <-
    min(data$LowerEnvelope[max(1,lengthData-SIDE_LENGTH_TO_CONSIDER_LOWER):lengthData])
  
  data$UpperEnvelope[1:min(nrow(data),SIDE_LENGTH_TO_CONSIDER_UPPER)] <- 
    max(data$UpperEnvelope[1:min(nrow(data),SIDE_LENGTH_TO_CONSIDER_UPPER)])
  data$UpperEnvelope[max(1,lengthData-SIDE_LENGTH_TO_CONSIDER_UPPER):lengthData] <-
    max(data$UpperEnvelope[max(1,lengthData-SIDE_LENGTH_TO_CONSIDER_UPPER):lengthData])
  
  # Plot if in DEBUG_MODE
  if( DEBUG_MODE ){
    
    require(ggplot2)
    if( use.active.power ){
      ggplot( data, aes(x=timestamp,y=active_power) ) + geom_line() + 
        ggtitle('Original')
    }else{
      ggplot( data, aes(x=timestamp,y=reactive_power) ) + geom_line() + 
        ggtitle('Original')
    }
    
    ggplot( data ) +
      geom_line(aes( x = timestamp, y = active_power,   colour = 'active power')) + 
      geom_line(aes( x = timestamp, y = LowerEnvelope,  colour = 'UpperEnvelope')) + 
      geom_line(aes( x = timestamp, y = UpperEnvelope,  colour = 'LowerEnvelope')) + 
      geom_line(aes( x = timestamp, y = reactive_power, colour = 'reactive power')) + 
      ggtitle('Envelopes')
    
    ggplot( data, aes(x=timestamp,y=diff.LowerEnvelope) ) + geom_line() +
      ggtitle('Impulsive Energy')
    
    # Find base power
    threshold = 5 # 5-percentile power level
    powerBaseValue  <- quantile( data$LowerEnvelope, probs = threshold/100 )
    powerBox        <- sum( data$LowerEnvelope - powerBaseValue )
    powerImpulse    <- sum( data$diff.LowerEnvelope )
    
    sprintf("Power Base : %.2f        Impulse : %.2f       Box : %.2f", 
            powerBaseValue*lengthData, powerImpulse, powerBox )
  }
  return(data)
}
