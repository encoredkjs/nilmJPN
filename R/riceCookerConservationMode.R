rice.cooker.conservation.mode <- function(data){
  
  box.reverse <- find.box.shape.reverse( data$active_power )
  box.shape   <- pmax( data$active_power - box.reverse, 0 )
  box.series  <- series.to.box.lists( data$timestamp, box.shape, 0 )
  
  n <- nrow(box.series[[1]])
  rice.cooker.logical <- which( (box.series[[1]]$duration < 15) & (box.series[[1]]$max.p < 800) )
  
  rice.cooker.logical <- rice.cooker.logical[min(which(diff(rice.cooker.logical) == 1)):length(rice.cooker.logical)]
  box.series[[1]] <- box.series[[1]][min(rice.cooker.logical):n,]
  
  if( max(rice.cooker.logical) < n )
  rice.cooker.logical <- rice.cooker.logical[ (rice.cooker.logical-min(rice.cooker.logical)) < 
                                                min(which( (box.series[[1]]$str[-1] - box.series[[1]]$end[-nrow(box.series[[1]])]) > 30 )) ] 
                                
  rice.cooker <- rbind.fill( box.series[[2]][ rice.cooker.logical ] )
  rice.cooker <- rbind( rice.cooker, 
                        data.frame( timestamp = data$timestamp[ !(data$timestamp %in% rice.cooker$timestamp) ],
                                    p = 0, q = 0))
  rice.cooker <- rice.cooker[ order(rice.cooker$timestamp), ]
  return(rice.cooker)
}

rice.cooker.conservation.mode.2nd.type <- function(data){
  
  data.old <- data
  tmp <- envelopeDetector(data)
  rice.cooker <- data.frame( timestamp = tmp$timestamp, p = tmp$diff.LowerEnvelope, q = 0 )
  box.series  <- series.to.box.lists( rice.cooker$timestamp, rice.cooker$p, 0 )
  rice.cooker <- rbind.fill( box.series[[2]][ which(box.series[[1]]$max.p > 1) ])
  rice.cooker <- rbind( rice.cooker, 
                        data.frame( timestamp = data$timestamp[ !(data$timestamp %in% rice.cooker$timestamp) ],
                                    p = 0, q = 0))
  rice.cooker <- rice.cooker[ order(rice.cooker$timestamp), ]
  data$active_power <- data$active_power - rice.cooker$p
  
  conservation.interval <- rice.cooker$p >= 10
  interval.str <- which( diff(conservation.interval) == diff(c(F,T)) ) + 1
  interval.end <- which( diff(conservation.interval) == diff(c(T,F)) ) 
  
  # 다음 신호가 10초 이내에 생성  
  conservation.interval[ unlist(sapply( interval.end, function(i) 
    if(any(interval.str < (i+10))) (i+1):max(interval.str[interval.str<(i+10)]) )) ] <- T
  interval.end <- which( diff(conservation.interval) == diff(c(T,F)) ) 
  
  conservation.interval[ unlist(lapply( interval.end, function(i) 
    (i+1):min((i+6),length(conservation.interval)) )) ] <- T
  
  interval.str <- which( diff(conservation.interval) == diff(c(F,T)) ) + 1
  interval.end <- which( diff(conservation.interval) == diff(c(T,F)) ) 
  if( length(interval.str) > 0 ){
    if( head(interval.end,1) < head(interval.str,1) ) interval.end <- interval.end[-1]    
  }
  if( length(interval.end) > 0 ){
    if( tail(interval.str,1) > tail(interval.end,1) ) interval.str <- interval.str[-length(interval.str)]  
  }
  if( length(interval.str) == 0 ) return(rice.cooker)
  
  if( any(interval.str[-1] - interval.end[-length(interval.end)] > 1000) ){
    cooker.end.pt <- min(interval.str[-1][ which(interval.str[-1] - interval.end[-length(interval.end)] > 1000) ])
    conservation.interval[ cooker.end.pt:length(conservation.interval) ] <- 0  
  }
  conservation.interval <- conservation.interval * 100
  
  box.reverse <- find.box.shape.reverse2(active.power=data$active_power, reactive.power=conservation.interval )[[1]]
  box.shape   <- pmax( data$active_power - box.reverse, 0 )

  rice.cooker$p <- rice.cooker$p + box.shape
  if( any(interval.str[-1] - interval.end[-length(interval.end)] > 1000) )
    rice.cooker$p[cooker.end.pt:nrow(rice.cooker)] <- 0
  
#   box.series  <- series.to.box.lists( data$timestamp, rice.cooker$p, 0 )
#   n <- nrow(box.series[[1]])
#   under.box.median <- median(box.series[[1]]$max.p)
#   if( any(under.box.median*1.2 < box.series[[1]]$max.p)){
#     for( i in which(under.box.median*1.2 < box.series[[1]]$max.p))
#       box.series[[2]][[i]]$p <- pmax( box.series[[2]][[i]]$p - max(box.series[[2]][[i]]$p) + under.box.median, 0 )
#   }else if( any(under.box.median*.8 > box.series[[1]]$max.p)){
#     for( i in which(under.box.median*.8 < box.series[[1]]$max.p))
#       box.series[[2]][[i]]$p <- pmax( box.series[[2]][[i]]$p - max(box.series[[2]][[i]]$p) + under.box.median, 0 )
#   }
#   rice.cooker <- rbind.fill( box.series[[2]] )
#   rice.cooker <- rbind( rice.cooker, 
#                         data.frame( timestamp = data$timestamp[ !(data$timestamp %in% rice.cooker$timestamp) ],
#                                     p = 0, q = 0))
#   rice.cooker <- rice.cooker[ order(rice.cooker$timestamp), ]
#   
  return(rice.cooker)
}

rice.cooker.conservation.mode.3rd.type <- function(data){
  
  data.old <- data
  tmp <- envelopeDetector(data)
  rice.cooker <- data.frame( timestamp = tmp$timestamp, p = data$active_power - tmp$LowerEnvelope, q = 0 )
  data$active_power <- data$active_power - rice.cooker$p
  
  return(rice.cooker)
}



rice.cooker.cooking.mode <- function(data, mode1.threshold = 1000, mode2.threshold = 500, 
                                     mode1.logical = T, mode2.logical = T, mode3.logical = T ){
  
  data.original <- data
  
  # first mode
  box.reverse <- find.box.shape.reverse( data$active_power, outer.threshold = mode1.threshold, inner.threshold = 50 )
  box.shape   <- pmax( data$active_power - box.reverse, 0 )
  if(!all(box.shape==0)){
    box.series  <- series.to.box.lists( data$timestamp, box.shape, 0 )
    rice.cooker.mode1.logical <- box.series[[1]]$str[2:(nrow( box.series[[1]]))] - box.series[[1]]$end[1:(nrow( box.series[[1]])-1)]
    if( any(rice.cooker.mode1.logical > 500) )
    rice.cooker.mode1.logical <- 1:(min(which(rice.cooker.mode1.logical > 500)))
    rice.cooker.mode1 <- rbind.fill( box.series[[2]][1:length(rice.cooker.mode1.logical)] )
  }else{
    rice.cooker.mode1 <- data.frame(timestamp=min(data$timestamp), p = 0, q = 0)
  }
  
  # third mode
  if( mode3.logical ){
    data <- data[ which(data$timestamp > max(rice.cooker.mode1$timestamp)),]
    data <- data[ 1:min(nrow(data),3600), ]
    box.reverse <- find.box.shape.reverse( data$active_power, outer.threshold = mode2.threshold )
    box.reverse <- envelopeDetector( data.frame( active_power = box.reverse, timestamp=data$timestamp))$LowerEnvelope
    box.shape   <- pmax( data$active_power - box.reverse, 0 )
    if(!all(box.shape==0)){
      box.series  <- series.to.box.lists( data$timestamp, box.shape, 0 )
      tmp <- rbind.fill( box.series[[2]][which(box.series[[1]]$max.p<1)] )
      box.shape[unlist(apply( box.series[[1]][box.series[[1]]$max.p<1,], 1, function(x){x[1]:x[2]}))] <- 0
      box.series  <- series.to.box.lists( data$timestamp, box.shape, 0 )
      rice.cooker.mode3 <- rbind.fill( box.series[[2]][1:(which.min(diff(box.series[[1]]$max.p))+1)] )
    }else{
      rice.cooker.mode3 <- data.frame(timestamp=min(data$timestamp), p = 0, q = 0)
    }
    
    # second mode
    if( mode2.logical ){
      buffer <- 10
      data <- data[ which(data$timestamp <= max(rice.cooker.mode3$timestamp)),]
      data$active_power[ which(data$timestamp %in% rice.cooker.mode3$timestamp) ] <- 
        data$active_power[ which(data$timestamp %in% rice.cooker.mode3$timestamp) ] - rice.cooker.mode3$p
      box.reverse <- find.box.shape.reverse( data$active_power )
      box.shape   <- pmax( data$active_power - box.reverse, 0 )
      if(!all(box.shape==0)){
        box.series  <- series.to.box.lists( data$timestamp, box.shape, 0 )
        rice.cooker.mode2 <- rbind.fill( box.series[[2]] )
      }else{
        rice.cooker.mode2 <- data.frame(timestamp=min(data$timestamp), p = 0, q = 0)
      }
    }else{
      rice.cooker.mode2 <- data.frame(timestamp=min(data$timestamp), p = 0, q = 0)
    }
  }else{
    rice.cooker.mode2 <- data.frame(timestamp=min(data$timestamp), p = 0, q = 0)
    rice.cooker.mode3 <- data.frame(timestamp=min(data$timestamp), p = 0, q = 0)
  }
  
  rice.cooker <- rbind( rice.cooker.mode1, rice.cooker.mode2, rice.cooker.mode3 )
  rice.cooker <- ddply( rice.cooker, .(timestamp), summarize, p=sum(p), q=sum(q))
  rice.cooker <- rbind( rice.cooker, 
                        data.frame( timestamp = data.original$timestamp[!(data.original$timestamp %in% rice.cooker$timestamp)],
                                    p = 0, q = 0))
  rice.cooker <- rice.cooker[ order(rice.cooker$timestamp), ]
  return(rice.cooker)
}

rice.cooker <- function(data, mode1.threshold = 1000, mode2.threshold = 500){
  
  first.step <- rice.cooker.cooking.mode(data, mode1.threshold, mode2.threshold)
  data <- data[ data$timestamp >= max(first.step$timestamp[ first.step$p > 0]), ]
  
  secon.step <- rice.cooker.conservation.mode(data)
  
  total.step <- rbind( first.step, secon.step )
  total.step <- ddply( total.step, .(timestamp), summarize, p=sum(p), q=sum(q))
}

rice.cooker.2nd.type <- function(data, mode1.threshold = 1000, mode2.threshold = 500,
                                 mode1.logical = T, mode2.logical = T, mode3.logical = T){
  
  first.step <- rice.cooker.cooking.mode(data, mode1.threshold, mode2.threshold,
                                         mode1.logical = mode1.logical, 
                                         mode2.logical = mode2.logical, 
                                         mode3.logical = mode3.logical)
  data <- data[ data$timestamp >= max(first.step$timestamp[ first.step$p > 0]), ]
  
  secon.step <- rice.cooker.conservation.mode.2nd.type(data)
  
  total.step <- rbind( first.step, secon.step )
  total.step <- ddply( total.step, .(timestamp), summarize, p=sum(p), q=sum(q))
}

rice.cooker.3rd.type <- function(data, mode1.threshold = 1000, mode2.threshold = 500,
                                 mode1.logical = T, mode2.logical = T, mode3.logical = T){
  
  first.step <- rice.cooker.cooking.mode(data, mode1.threshold, mode2.threshold,
                                         mode1.logical = mode1.logical, 
                                         mode2.logical = mode2.logical, 
                                         mode3.logical = mode3.logical)
  data <- data[ data$timestamp >= max(first.step$timestamp[ first.step$p > 0]), ]
  
  secon.step <- rice.cooker.conservation.mode.3rd.type(data)
  
  total.step <- rbind( first.step, secon.step )
  total.step <- ddply( total.step, .(timestamp), summarize, p=sum(p), q=sum(q))
}