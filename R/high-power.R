high.power.detect <- function(data, min.threshold, max.iter = 3, debug.mode = F ){
  generate.HeavyLoad.meta(data=data,min.threshold=min.threshold,max.iter=max.iter, debug.mode=debug.mode)  
}

generate.HeavyLoad.meta <- function( data, min.threshold, max.iter = 3, debug.mode = F ){
  
  preprocessed.data <- data.frame( timestamp      = data$timestamp, 
                                   active_power   = remove.consecutive.jump.pt2(data$active_power, level=2),
                                   reactive_power = remove.consecutive.jump.pt2(data$reactive_power, level=2))
  
  if( !missing(min.threshold) ) max.iter <- 1
  
  meta <- list()
  para <- NULL
  high.power.info.total <- data.frame()
  for( iter in 1:max.iter ){
    
    if( missing(min.threshold) | max.iter > 1 ){
      
      least.max <- 300
      if( diff(range(preprocessed.data$active_power)) < least.max ) return(meta)
      data.lists <- split.data( preprocessed.data )
      r.edge.height <- laply( data.lists, 
                              function(x){ edge.height <- sort(pmax(least.max,diff(x$active_power)))
                                           return(mean(edge.height[which.max(diff(edge.height)) + c(0,1)])) })
      f.edge.height <- laply( data.lists, 
                              function(x){ edge.height <- sort(pmin(-least.max,diff(x$active_power)))
                                           return(mean(edge.height[which.max(diff(edge.height)) + c(0,1)])) })
      
      data.lists[which(abs(r.edge.height + f.edge.height) > least.max)] <- NULL
      edge.height <- unlist(lapply( data.lists, 
                                    function(x) sort(pmax(least.max,abs(diff(x$active_power)))) ), use.names=F)
      if( length(meta) > 0 & is.null(para) ) edge.height <- edge.height[ edge.height < min.threshold ]
      min.threshold <- mean(edge.height[which.max(diff(edge.height)) + c(0,1)])
      print(min.threshold)
      if( min.threshold < least.max ) break
      
    } 
    
    box.reverse <- find.box.shape.reverse( power = preprocessed.data$active_power, 
                                           recursive = T, 
                                           inner.threshold=min.threshold/10, 
                                           outer.threshold=min.threshold)
    box.result <- preprocessed.data$active_power - box.reverse
    under.peak <- data$active_power - preprocessed.data$active_power
    under.peak[ box.result <= 0 ] <- 0
    
    box.lists <-  series.to.box.lists( data$timestamp, box.result, 0 )[[1]]
    box.temp  <- box.result + under.peak
    
    under.peak.correction <- 
      data.frame( value = box.temp[which( box.temp < 0 )],
                  idx   = sapply( which( box.temp < 0 ), 
                                  function(x) which((box.lists$str <= x) & (x <= box.lists$end)) ))
    under.peak.correction <- ddply( under.peak.correction, .(idx), summarize, value = min(value))
    
    box.result <- box.temp
    if( nrow(under.peak.correction) != 0 ){
      for( i in 1:nrow(under.peak.correction) ){
        str.ind <- box.lists$str[under.peak.correction$idx[i]]
        end.ind <- box.lists$end[under.peak.correction$idx[i]]
        box.result[ str.ind:end.ind ] <- 
          box.result[ str.ind:end.ind ] + abs(under.peak.correction$value[i])
      }  
    }
    box.lists <- box.lists[ box.lists$str > 15, ]
    box.lists <- box.lists[ box.lists$end < length(data$timestamp) - 15, ]
    
    ap.rising <- sapply( box.lists$str, function(t) 
      find.stationary.state( box.result[t+c(0:15)], forward=T)[[2]] - 
        find.stationary.state( box.result[t - c(15:0)], forward=F)[[2]] )
    rp.rising <- sapply( box.lists$str, function(t) 
      find.stationary.state( data$reactive_power[t+c(0:15)], forward=T)[[2]] - 
        find.stationary.state( data$reactive_power[t - c(15:0)], forward=F)[[2]] )
    ap.falling <- sapply( box.lists$end, function(t) 
      find.stationary.state( box.result[t+c(0:15)], forward=T)[[2]] - 
        find.stationary.state( box.result[t - c(15:0)], forward=F)[[2]] )
    rp.falling <- sapply( box.lists$end, function(t) 
      find.stationary.state( data$reactive_power[t+c(0:15)], forward=T)[[2]] - 
        find.stationary.state( data$reactive_power[t - c(15:0)], forward=F)[[2]] )
    
    high.power.info <- data.frame( str = data$timestamp[ box.lists$str ],
                                   end = data$timestamp[ box.lists$end ],
                                   ap.rising, rp.rising, ap.falling, rp.falling )
    high.power.info <- subset( high.power.info, 
                               (ap.rising > min.threshold) | (ap.falling < -min.threshold) )
    
    if( nrow(high.power.info) != 0 ){
      
      if( nrow(high.power.info.total) == 0 ){
        high.power.info.total <- high.power.info
      }else{
        high.power.info.total <- rbind(high.power.info.total, high.power.info)
      }
      
      box.reverse[ data$timestamp < high.power.info$str[1] ] <- 
        data$active_power[ data$timestamp < high.power.info$str[1] ]
      if( nrow(high.power.info) > 1){
        for( i in 1:(nrow(high.power.info)-1) ){
          box.reverse[ data$timestamp %within% interval(high.power.info$end[i],high.power.info$str[i+1]) ] <- 
            data$active_power[ data$timestamp %within% interval(high.power.info$end[i],high.power.info$str[i+1]) ]
        }
      }
      box.reverse[ data$timestamp > high.power.info$end[nrow(high.power.info)] ] <- 
        data$active_power[ data$timestamp > high.power.info$end[nrow(high.power.info)] ]
      data$active_power <- box.reverse
      
      para <- search.heavy.load.meta( high.power.info, debug.mode = debug.mode, bin.threshold=5)
      print( para )
      #if(debug.mode) plot( high.power.info.total$ap.rising, high.power.info.total$rp.rising )
      if( !is.null(para) ){
        para.json <- lapply( para, function(x) 
          from.para.to.JSON2( data, r.edge=x[[1]], f.edge=x[[2]], box.no = x[[3]] ))
        meta <- append( meta, para.json )
        names(meta) <- 1:length(meta)
        
        high.power <- data.frame( timestamp = data$timestamp, active_power = box.result )
        high.power$logical <- high.power$active_power < min.threshold
        remove.data <- ddply( high.power, .(day=day(timestamp)), summarize, logical = all(logical))
        high.power <- subset( high.power, day(timestamp) %in% remove.data$day[!remove.data$logical] )
        
        if( debug.mode ){
          tmp <- ddply( high.power, .(day=as.Date(timestamp, tz=tz(timestamp[1]))), summarize, 
                        str.t = min(timestamp[!logical]), 
                        end.t = max(timestamp[!logical]))
          tmp$str.t <- tmp$str.t - 10
          tmp$end.t <- tmp$end.t + 10
          high.power <- rbind.fill( apply(tmp,1,function(x) subset(high.power,timestamp %within% interval(x['str.t'],x['end.t']))))
          high.power$min30 <- floor_date(high.power$timestamp,'hour') + 
            minutes((minute(high.power$timestamp) %/% 30) * 30)
          remove.data <- ddply( high.power, .(min30), summarize, active_power = sum(active_power) )
          high.power <- subset( high.power, !(min30 %in% remove.data$min30[ remove.data$active_power == 0 ]) )
          
          if( nrow(high.power) > 0 ) {
            high.power$facet <- cumsum(c(TRUE,diff(high.power$timestamp) > 60))
            print(ggplot( high.power, aes(x=timestamp, y=active_power)) + 
                    geom_line() + facet_wrap(~facet, scales='free'))
          }
        }
        #data$facet <- floor_date( data$timestamp, 'day')
        #ggplot( data, aes(x=timestamp, y=active_power)) + geom_line() + 
        #  facet_wrap(~facet, scales='free')
        
        preprocessed.data$active_power <- remove.consecutive.jump.pt2(box.reverse)
        
      }    
    }else para <- NULL
  }
  return( meta )
}


# split.data1 <- function( data, probs = .05, debug.mode = F, level = 2, apply.envelopeDetector = F ){
#   
#   data.original <- data
#   if( apply.envelopeDetector ){
#     tmp <- envelopeDetector( data, use.active.power = T ); data$active_power <- tmp$LowerEnvelope
#   }
#   
#   n <- nrow(data)
#   baseline              <- data.frame(over.baseline  = (data$active_power > quantile(data$active_power, probs) + 10))
#   baseline$state.change <- c(ifelse( baseline$over.baseline[1], 1, 0),diff(baseline$over.baseline))
#   baseline$block.ind    <- cumsum(abs(baseline$state.change))
#   baseline$global.ind   <- 1:n
#   
#   block.summary <- ddply( baseline, .(block.ind,over.baseline), summarize, 
#                           str = min(global.ind), end = max(global.ind))
#   
#   merge.interval <- which( !block.summary$over.baseline & (block.summary$end - block.summary$str) < 20 )
#   if( head(merge.interval,1) ==                   1 ) block.summary$str[2] <- 1
#   if( tail(merge.interval,1) == nrow(block.summary) ) block.summary$end[nrow(block.summary)-1] <- n
#   merge.interval <- merge.interval[ (1 < merge.interval) & (merge.interval < nrow(block.summary)) ]
#   block.summary$end[ merge.interval[ c(0,which(diff(merge.interval) > 4)) + 1] - 1 ] <- 
#     block.summary$end[ merge.interval[ c(which(diff(merge.interval) > 4),length(merge.interval))] + 1 ]
#   block.summary <- block.summary [ - unlist(mapply( function(s,e) s:e, 
#                                                     merge.interval[ c(0,which(diff(merge.interval) > 4)) + 1],
#                                                     merge.interval[ c(which(diff(merge.interval) > 4),length(merge.interval))] + 1 )),]
#   block.subset <- subset( block.summary, over.baseline )
#   
#   data.lists   <- apply( block.subset, 1, function(x){ ind <- (x['str']-10):(x['end']+10)
#                                                        ind <- ind[(1 <= ind)&(ind <= n)]
#                                                        data <- data[ind,]} )
#   
#   need.subdivision <- sapply( data.lists, function(x) difftime(max(x$timestamp),min(x$timestamp),units='hours')) > 12
#   if( length(which(need.subdivision)) > 0  & level > 0 ){
#     recursive.d <- lapply( data.lists[ need.subdivision ], 
#                            function(df) split.data1(df,probs=probs,debug.mode=debug.mode,level=(level-1) ))
#     data.lists  <- append( data.lists[ !need.subdivision ], unlist(recursive.d,recursive=F) )
#   }
#   
#   ap.range   <- sapply( data.lists, function(x) diff(range(x$active_power)) )
#   
#   data.lists[ (ap.range < 10) | (ap.range > 500) ] <- NULL 
#   data.lists[ sapply( data.lists, nrow) < 60 ] <- NULL  
#   
#   if( debug.mode & level == 2 ){
#     period <- ldply( data.lists, function(x) data.frame('str'=min(x$timestamp), 'end'=max(x$timestamp)) )
#     gaps   <- which(tail( period$str,-1 ) - head( period$end, -1 ) > 1)
#     period <- data.frame( str = c(head(period$str,1), period$str[gaps+1]),
#                           end = c(period$end[gaps],tail(period$end,1) ))
#     data.subset <- adply( period, 1, function(x) subset(data,timestamp %within% interval(x$str,x$end) ))
#     data.subset <- merge( data.frame(timestamp=data$timestamp),data.subset, all.x=T)
#     
#     par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
#     plot( data$timestamp, data.original$active_power, type='l', xaxt='n')
#     abline( v = period$str, col='red'  )
#     abline( v = period$end, col='blue' )
#     plot( data.subset$timestamp, data.subset$active_power, type='l', xlab='Timestamp', ylab='Active power')
#     title(paste('Step0 : 사용량이 낮은 구간 찾기\n', min(data$timestamp), '--', max(data$timestamp)), 
#           outer=TRUE)
#     mtext("Timestamp", 1, 2, outer=TRUE)
#     mtext("After / Before", 2, 3, outer=TRUE, las=0)
#     par(mfrow=c(1,1))
#   }
#   return(data.lists)
# }

split.data1 <- function( data, probs = .05, debug.mode = F, level = 2, apply.envelopeDetector = F ){
  
  data.original <- data
  if( apply.envelopeDetector ){
    tmp <- envelopeDetector( data, use.active.power = T ); data$active_power <- tmp$LowerEnvelope
  }
  
  n <- nrow(data)
  baseline              <- data.frame(over.baseline  = (data$active_power > quantile(data$active_power, probs) + 10))
  baseline$state.change <- c(ifelse( baseline$over.baseline[1], 1, 0),diff(baseline$over.baseline))
  baseline$block.ind    <- cumsum(abs(baseline$state.change))
  baseline$global.ind   <- 1:n
  
  block.summary <- ddply( baseline, .(block.ind,over.baseline), summarize, 
                          str = min(global.ind), end = max(global.ind))
  block.summary$divide <- mapply( function(s,e) which.min(data$active_power[s:e]), 
                                  block.summary$str, block.summary$end ) + block.summary$str - 1 
  
  move.str.pt <- which( !block.summary$over.baseline )
  move.str.pt <- move.str.pt[ move.str.pt + 1 <= nrow(block.summary) ]
  block.summary$str[ move.str.pt + 1 ] <- block.summary$divide[ move.str.pt ] + 1
  
  move.end.pt <- which( !block.summary$over.baseline )
  move.end.pt <- move.end.pt[ move.end.pt - 1 >= 1 ]
  block.summary$end[ move.end.pt - 1 ] <- block.summary$divide[ move.end.pt ]
  
  block.subset <- subset( block.summary, over.baseline )
  data.lists   <- apply( block.subset, 1, function(x) data[x['str']:x['end'],] )
  
  need.subdivision <- sapply( data.lists, function(x) difftime(max(x$timestamp),min(x$timestamp),units='hours')) > 12
  if( length(which(need.subdivision)) > 0  & level > 0 ){
    recursive.d <- lapply( data.lists[ need.subdivision ], 
                           function(df) split.data1(df,probs=probs,debug.mode=debug.mode,level=(level-1) ))
    data.lists  <- append( data.lists[ !need.subdivision ], unlist(recursive.d,recursive=F) )
  }
  
  ap.range   <- sapply( data.lists, function(x) diff(range(x$active_power)) )
  
  data.lists[ (ap.range < 10) | (ap.range > 500) ] <- NULL 
  data.lists[ sapply( data.lists, nrow) < 60 ] <- NULL  
  
  if( debug.mode & level == 2 ){
    period <- ldply( data.lists, function(x) data.frame('str'=min(x$timestamp), 'end'=max(x$timestamp)) )
    gaps   <- which(tail( period$str,-1 ) - head( period$end, -1 ) > 1)
    period <- data.frame( str = c(head(period$str,1), period$str[gaps+1]),
                          end = c(period$end[gaps],tail(period$end,1) ))
    data.subset <- adply( period, 1, function(x) subset(data,timestamp %within% interval(x$str,x$end) ))
    data.subset <- merge( data.frame(timestamp=data$timestamp),data.subset, all.x=T)
    
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot( data$timestamp, data.original$active_power, type='l', xaxt='n')
    abline( v = period$str, col='red'  )
    abline( v = period$end, col='blue' )
    plot( data.subset$timestamp, data.subset$active_power, type='l', xlab='Timestamp', ylab='Active power')
    title(paste('Step0 : 사용량이 낮은 구간 찾기\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 2, outer=TRUE)
    mtext("After / Before", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  return(data.lists)
}


split.data <- function( data, least.max.threshold = 500, interval.length = 5 ){
  
  baseline <- data.frame( over.baseline  = (data$active_power > quantile( data$active_power, .05)), 
                          over.least.max = (data$active_power > least.max.threshold) )
  
  baseline$state.change <- F; baseline$state.change[1] <- T
  baseline$state.change[ which(diff(baseline$over.baseline)!=0) + 1 ] <- T
  baseline$ind <- cumsum(baseline$state.change)
  
  baseline.summary <- ddply( baseline, .(ind), summarize, log = all(over.baseline) & any(over.least.max))
  
  subinterval     <- data.frame(str = which(baseline$state.change))
  subinterval$end <- c( subinterval$str[-1] - 1, nrow(data)) 
  subinterval     <- subinterval[ baseline.summary$ind[ baseline.summary$log ], ]
  
  data.lists <- apply( subinterval, 1, function(x) data[x['str']:x['end'],] )
  return(data.lists)
  
}

find.box.shape.reverse5 <- function( data, json.para, threshold = 10, level = 2, recursive = T, 
                                     debug.mode = F ){
  
  ap <- data$active_power
  rp <- data$reactive_power
  
  if( (length(ap) < 5) || (level==0) ) return(ap)
  
  r.edge.parameters <- json.para$rising_edge
  f.edge.parameters <- json.para$falling_edge
  
  ap.clean    <- remove.consecutive.jump.pt2( ap, threshold, level = 2 )
  rp.clean    <- remove.consecutive.jump.pt2( rp, threshold, level = 2 )
  ap.diff     <- diff( ap.clean )
  rp.diff     <- diff( rp.clean )
  box.reverse <- ap.clean
  
  state.change <- function( power, i ){
    rite.ind <- i + c(0:15); rite.ind <- rite.ind[ rite.ind <= length(power) ]
    left.ind <- i - c(0:15); left.ind <- left.ind[ left.ind >=             1 ]
    stationary.rite <- find.stationary.state(power[ rite.ind ], forward=T)
    stationary.left <- find.stationary.state(power[ left.ind ], forward=F)
    edge.length     <- stationary.rite[[2]] - stationary.left[[2]]
    stationary.pt   <-  i + c( -stationary.left[[1]], stationary.rite[[1]] )
    return(c(edge.length,stationary.pt))
  }
  
  rp.is.negligible <- abs(r.edge.parameters['rp_height']) < 30
  if( !rp.is.negligible ){
    r.lower.bound <- qnorm( .01, r.edge.parameters['rp_height'], r.edge.parameters['rp_sigma'] )
    r.upper.bound <- qnorm( .99, r.edge.parameters['rp_height'], r.edge.parameters['rp_sigma'] )
    f.lower.bound <- qnorm( .01, f.edge.parameters['rp_height'], f.edge.parameters['rp_sigma'] )
    f.upper.bound <- qnorm( .99, f.edge.parameters['rp_height'], f.edge.parameters['rp_sigma'] )
    
    r.edge.candidate <- NULL
    f.edge.candidate <- NULL
    if( r.edge.parameters['rp_height'] > 0 ){
      r.edge.candidate      <- which( rp.diff > r.lower.bound / 3 )
      f.edge.candidate      <- which( rp.diff < pmin(f.upper.bound / 3, 0) )
      
      if( is.null(r.edge.candidate) | is.null(f.edge.candidate) ){
        return( list( data$active_power, r.edge, f.edge) )
      }

      tmp <- sapply( r.edge.candidate, function(x) state.change(rp.clean, x) )
      tmp[1, which((diff(tmp[2,]) == 0) & (diff(tmp[3,]) == 0)) ] <- Inf
      r.edge.candidate      <- r.edge.candidate[ which( tmp[1,] > r.lower.bound ) ]
      
      tmp <- sapply( f.edge.candidate, function(x) state.change(rp.clean, x) )
      tmp[1, which((diff(tmp[2,]) == 0) & (diff(tmp[3,]) == 0)) ] <- Inf
      f.edge.candidate      <- f.edge.candidate[ which( tmp[1,] < f.upper.bound ) ]
    }else{
      r.edge.candidate      <- which( rp.diff < r.upper.bound / 3 )
      f.edge.candidate      <- which( rp.diff > f.lower.bound / 3 )
      
      if( (length(r.edge.candidate)==0) | (length(f.edge.candidate) == 0) ){
        return( data$active_power )#return( list( data$active_power, r.edge, f.edge) )
      }

      tmp <- sapply( r.edge.candidate, function(x) state.change(rp.clean, x) )
      tmp[1, which((diff(tmp[2,]) == 0) & (diff(tmp[3,]) == 0)) ] <- Inf
      r.edge.candidate      <- r.edge.candidate[ which( tmp[1,] < r.upper.bound ) ]
      
      tmp <- sapply( f.edge.candidate, function(x) state.change(rp.clean, x) )
      tmp[1, which((diff(tmp[2,]) == 0) & (diff(tmp[3,]) == 0)) ] <- Inf
      f.edge.candidate      <- f.edge.candidate[ which( tmp[1,] > f.lower.bound ) ]
    }  
  }
  
  r.lower.bound <- qnorm( .01, r.edge.parameters['ap_height'], r.edge.parameters['ap_sigma'] )
  r.upper.bound <- qnorm( .99, r.edge.parameters['ap_height'], r.edge.parameters['ap_sigma'] )
  f.lower.bound <- qnorm( .01, f.edge.parameters['ap_height'], f.edge.parameters['ap_sigma'] )
  f.upper.bound <- qnorm( .99, f.edge.parameters['ap_height'], f.edge.parameters['ap_sigma'] )
  
  r.edge      <- which( ap.diff > r.lower.bound / 3 )
  f.edge      <- which( ap.diff < f.upper.bound / 3 )
  
  if( length(r.edge) == 0 | length(f.edge) == 0 ) return( list( data$active_power, integer(), integer()) )
  
  tmp     <- sapply( r.edge, function(x) state.change(ap.clean, x) )
  tmp[1, which((diff(tmp[2,]) == 0) & (diff(tmp[3,]) == 0)) ] <- -Inf
  r.edge  <- r.edge[ which( tmp[1,] > r.lower.bound ) ]
  
  tmp     <- sapply( f.edge, function(x) state.change(ap.clean, x) )
  tmp[1, which((diff(tmp[2,]) == 0) & (diff(tmp[3,]) == 0)) ] <- Inf
  f.edge  <- f.edge[ which( tmp[1,] < f.upper.bound ) ]
  
  if( !rp.is.negligible ){
    # reactive power가 변하는 지점 기준으로 +-5 이내로 active power가 변한다고 가정
    r.edge <- r.edge[ r.edge %in% outer(r.edge.candidate, c(-5:5),'+') ]
    f.edge <- f.edge[ f.edge %in% outer(f.edge.candidate, c(-5:5),'+') ]
  }
  r.edge <- r.edge[ r.edge < max(f.edge) ]
  f.edge <- f.edge[ f.edge > min(r.edge) ]
  
  if( debug.mode ){  
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(data$active.power, type='l', ann=FALSE,   xaxt="n") 
    abline( v=r.edge, col='red' ); abline( v=f.edge, col='blue')
    plot(data$reactive.power, type='l', ann=FALSE,   xaxt="n") 
    abline( v=r.edge, col='red' ); abline( v=f.edge, col='blue')
    title(paste('Heavy load step1 : Edge detection 결과\n', 
                min(data$timestamp), '--', max(data$timestamp)), outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  
  #   r.edge.exist <- rep( T, length(r.edge) )
  #   for( i in 1:(length(r.edge)-1) ){
  #     if( !any((r.edge[i] < f.edge) & (f.edge < r.edge[i+1])) ) r.edge.exist[i] <- F
  #   }
  #   r.edge.exist[ r.edge > max(f.edge) ] <- F; r.edge <- r.edge[ r.edge.exist ]
  #   
  #   f.edge.exist <- rep( T, length(f.edge) )
  #   for( i in 1:(length(f.edge)-1) ){
  #     if( !any((f.edge[i] < r.edge) & (r.edge < f.edge[i+1])) ) f.edge.exist[i] <- F
  #   }
  #   f.edge.exist[ f.edge < min(r.edge) ] <- F; f.edge <- f.edge[ f.edge.exist ]
  
  jump.position <- sort(c( r.edge, f.edge ))
  if( length(jump.position) == 0 ) return( list( data$active_power, integer(), integer()) )
  #edge          <- data.frame( re = r.edge, fe = f.edge )
  #if( nrow(edge) == 0 ) return(ap)
  
  increase <- integer()
  decrease <- integer()
  interval <- integer()
  jump.str <- integer()
  jump.end <- integer()
  for( i in 1:length(jump.position) ){
    if( sign(ap.diff[jump.position[i]]) == 1){
      if( i == length(jump.position) ){
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.left <- jump.position[i]   - stationary.left
        box.reverse[ (stationary.left+1):length(box.reverse) ] <- box.reverse[stationary.left]
      }else if( sign(ap.diff[jump.position[i]]) != sign(ap.diff[jump.position[i+1]]) ){
        increase <- c( increase,ap.diff[jump.position[i]] )
        decrease <- c( decrease,ap.diff[jump.position[i+1]] )
        interval <- c( interval,jump.position[i+1] - jump.position[i] )
        jump.str <- c( jump.str, jump.position[i]+1)
        jump.end <- c( jump.end, jump.position[i+1])
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        rite.ind <- jump.position[i+1]:(jump.position[i+1]+14);  rite.ind <- pmin(rite.ind,length(ap))
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F,threshold=threshold) 
        left.pt <- stationary.left[[2]]
        stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T,threshold=threshold) 
        rite.pt <- stationary.rite[[2]]
        stationary.left <- jump.position[i]   - stationary.left[[1]]
        stationary.rite <- jump.position[i+1] + stationary.rite[[1]]
        
        subinterval.ind <- c((stationary.left+1):(stationary.rite-1))
        subinterval <- ap.clean[subinterval.ind]
        
        test <- find.stationary.state( subinterval,forward=T, threshold=threshold) 
        remove.left.pt <- test[[1]]
        if( remove.left.pt > 0 ) subinterval[1:remove.left.pt] <- test[[2]]
        test <- find.stationary.state( subinterval,forward=F, threshold=threshold) 
        remove.rite.pt <- test[[1]]
        if( remove.rite.pt > 0 ) subinterval[(length(subinterval)+1) - (1:remove.rite.pt)] <- test[[2]]
        
        box.reverse[ subinterval.ind ] <- min(c(left.pt,rite.pt))
        box.reverse[ box.reverse > ap ] <- ap[ box.reverse > ap ]
        #         if( remove.left.pt > 0 ){
        #           subinterval.ind <- subinterval.ind[-(1:remove.left.pt)]
        #           subinterval     <- subinterval[-(1:remove.left.pt)]
        #         } 
        #         if( remove.rite.pt > 0 ){
        #           subinterval.ind <- subinterval.ind[-(length(subinterval.ind+1)-(1:remove.rite.pt))]
        #           subinterval     <- subinterval[-(length(subinterval+1)-(1:remove.rite.pt))]
        #         } 
        
        if( recursive ){
          if(level==1){
            box.reverse[ which(box.reverse > ap) ] <- ap[ which(box.reverse > ap) ]
            return( box.reverse )
          } 
          if(length(subinterval)>5)
            if( max(abs(diff(subinterval))) > threshold ){
              while( length(which(abs(diff(subinterval)) > threshold)) > 0 ){
                subinterval.old <- subinterval
                subbox <- ( subinterval - find.box.shape.reverse( subinterval, level = (level-1) ))
                if(subbox[1]>0){
                  #  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox - subbox[1], 0 )
                  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox , 0 )
                }else if(subbox[length(subbox)]>0){
                  #  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox - subbox[length(subbox)], 0 )
                  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox , 0 )
                }else{
                  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox, 0 )
                }
                subinterval <- subinterval - subbox
                if( all(subinterval-subinterval.old == 0) ) break
              }
            } 
        }
      } 
    }else if( i == 1 & jump.position[1] < length(ap)) 
      box.reverse[ 1:(jump.position[1])] <- ap.clean[jump.position[1]+1]
  }
  box.reverse[ box.reverse > ap ] <- ap[ box.reverse > ap ]
  
  box.shape <- pmax(ap.clean - box.reverse,0)
  box.lists <- series.to.box.lists( data$timestamp, box.shape, 0 )[[1]]
  box.lists <- box.lists[ which(box.lists$end - box.lists$str > 5), ]
  
  if( nrow(box.lists) == 0 ) return( ap )
  box.reverse[ 1:(box.lists$str[1]-1) ] <- ap[ 1:(box.lists$str[1]-1) ]
  if( nrow(box.lists) > 1){
    for( i in 1:(nrow(box.lists)-1)){
      box.reverse[ (box.lists$end[i]+1):(box.lists$str[i+1]-1) ] <- 
        ap[ (box.lists$end[i]+1):(box.lists$str[i+1]-1) ]
    }
  }
  box.reverse[ (box.lists$end[nrow(box.lists)]+1):length(box.reverse) ] <- 
    ap[ (box.lists$end[nrow(box.lists)]+1):length(box.reverse) ]
  
  rising.edge  <- ddply( box.lists, .(str), function(x){find.jump.value(data,x$str)} )
  box.lists    <- box.lists[ r.lower.bound > rising.edge$ap | rising.edge$ap > r.upper.bound, ]
  
  falling.edge <- ddply( box.lists, .(end), function(x){find.jump.value(data,x$end)} )
  box.lists    <- box.lists[ f.lower.bound > falling.edge$ap | falling.edge$ap > f.upper.bound, ]
  
  
  #box.lists <- subset(box.lists, (r.lower.bound > mode.p) | (mode.p > r.upper.bound))
  if( nrow(box.lists) > 0 ){
    for( i in 1:nrow(box.lists) ){
      box.reverse[ box.lists$str[i]:box.lists$end[i] ] <- ap[ box.lists$str[i]:box.lists$end[i] ]
    }
  }
  return( list( box.reverse, r.edge, f.edge) )
}
