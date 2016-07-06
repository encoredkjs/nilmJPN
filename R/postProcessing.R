period.index     <- function(eventTime, total.period){
  eventTime.diff <- diff(eventTime)
  if( missing(total.period) ) total.period   <- median(eventTime.diff)
  return(list( total.period = total.period, idx = round( eventTime.diff / total.period )))
}

undetected.event <- function(eventTime, str.idx, end.idx, total.period){
  
  if( missing(str.idx) || str.idx > head(eventTime,1) ) str.idx <- head(eventTime,1)
  if( missing(end.idx) || end.idx < tail(eventTime,1) ) end.idx <- tail(eventTime,1)
  
  event.no   <- period.index( eventTime, total.period )[['idx']]
  period.t   <- period.index( eventTime, total.period )[['total.period']]
  
  undetected <- list()
  if( any(event.no != 1) )
    undetected <- lapply( which(event.no > 1), function(i) seq( eventTime[i], eventTime[i+1], length.out=event.no[i]+1))
  
  if( head(event.no,1) > str.idx ) undetected <- append( undetected, seq(head(eventTime,1),str.idx,by=-period.t))
  if( tail(event.no,1) < end.idx ) undetected <- append( undetected, seq(tail(eventTime,1),end.idx,by= period.t))
  
  undetected <- round((unlist(undetected)))
  undetected <- setdiff( undetected, eventTime )
  
  return(undetected)
}

total.event <- function(eventTime, total.period){ 
  
  if( length(eventTime) == 1 ) return( data.frame( eventTime = eventTime, is.detected = T) )
  
  result     <- data.frame( eventTime = eventTime, is.detected = T )
  undetected <- undetected.event(eventTime,total.period = total.period)
  if( length(undetected) == 0 ) return(result)
  
  result     <- rbind( result, data.frame( eventTime = undetected, is.detected = F ) )
  result     <- result[ order(result$eventTime), ]
  return(result)
}

post.processing <- function(meta, data, debug.mode = F, calcWTbound = T){
  
  if(!('ap.box' %in% names(data))){
    data$ap.box <- data$p
    data$rp.box <- data$q
    data$ap.box.reverse <- 0
    data$rp.box.reverse <- 0
  }
  tmp         <- series.to.box.lists( data$timestamp, data$ap.box, data$rp.box )
  box.summary <- tmp[[1]]
  box.summary <- box.summary[sapply( box.summary, function(x) !is.null(x) )]
  box.summary <- box.summary[ box.summary$max.p >= 35,]
  
  if( all(data$ap.box==0) || nrow(box.summary) <= 1 ){
    return( data )
  } 
    
  total.period   <- period.index( box.summary$str )$total.period
  
  if( calcWTbound ){
    wt.lower.bound <- qnorm( .01, meta$cycle['working_time'], meta$cycle['working_time_sigma'] )
    wt.upper.bound <- qnorm( .99, meta$cycle['working_time'], meta$cycle['working_time_sigma'] )
    
    if( wt.lower.bound < 0){
      wt.lower.bound <- boxplot.stats(box.summary$duration)$conf[1]
    }
    #  box.summary    <- subset(box.summary, ((wt.lower.bound <= duration) & 
    #                                           (duration <= wt.upper.bound)) |
    #                             (abs(c(NA, diff(str)) / total.period - 1) < .2))
    box.summary    <- subset(box.summary, ((wt.lower.bound <= duration) & 
                                             (duration <= wt.upper.bound))) 
    box.summary    <- box.summary[ order(box.summary$str),]
  }
  
  n <- nrow(box.summary)
  
  box.series  <- tmp[[2]][box.summary$no]
  clean.box   <- rbind.fill(box.series)


  post.processing.data <- data.frame( 'timestamp' = data$timestamp
                                      , 'ap.box.reverse' = data$ap.box.reverse
                                      , 'rp.box.reverse' = data$rp.box.reverse )
  post.processing.data <- merge( post.processing.data, clean.box, all.x = T ) 
  post.processing.data[ is.na(post.processing.data) ] <- 0
  
  box.series  <- ldply( box.series, function(x){x$idx <- 1:nrow(x); return(x)})
  
  p.value <- unique(box.series$p)
  q.value <- unique(box.series$q)
  
  if( length(p.value) == 1 ){
    box.series  <- data.frame( idx = 1:pmin(round(meta$cycle['working_time']), nrow(box.series)),
                               p.median = p.value )
  }else{
    box.series  <- ddply( box.series, .(idx), summarize, p.median = median(p[ !(p %in% boxplot.stats(p)$out)]))
    box.series  <- box.series[1:pmin(round(meta$cycle['working_time']), nrow(box.series)),]
  }
  total.period   <- period.index( box.summary$str )$total.period
  total.period   <- meta$cycle['tot_time']
  ind <- box.summary$end[1] - length(box.series$p.median):0
  while( any(ind >= 1) ){
    ind <- ind - total.period
    ind <- ind[ ind >= 1]
    if( length(p.value) != 1 ) post.processing.data$p[ ind ] <- rev(rev(box.series$p.median)[1:length(ind)])
    if( length(p.value) == 1 ) post.processing.data$p[ ind ] <- p.value
  }
    
  ind <- box.summary$str[n] + 0:length(box.series$p.median)
  while( any(ind <= nrow(data)) ){
    ind <- ind + total.period
    ind <- ind[ ind <= nrow(data)]
    if( length(p.value) != 1 ) post.processing.data$p[ ind ] <- box.series$p.median[1:length(ind)]
    if( length(p.value) == 1 ) post.processing.data$p[ ind ] <- p.value
  }
  
  add.box <- undetected.event(box.summary$str)
  if( length(add.box) > 0 ){
    for( i in 1:length(add.box) ){
      str.new <- add.box[i]
      end.new <- pmin( str.new + nrow(box.series) - 1, nrow(data) )        
      if( length(p.value) != 1 ) post.processing.data$p[str.new:end.new] <- box.series$p.median[ 1:(end.new-str.new+1)]
      if( length(p.value) == 1 ) post.processing.data$p[str.new:end.new] <- p.value
    }
  }

  post.processing.data$p[ is.na(post.processing.data$p) ] <- 0

  if( debug.mode ){
        
    compare <- merge( post.processing.data, data )
    par(mfrow=c(4,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$ap.box, ann=FALSE, xaxt="n", type='l')
    plot(x=post.processing.data$timestamp, 
         y=post.processing.data$p, ann=FALSE, xaxt="n", type='l')
    plot( compare$timestamp, pmax( compare$p - compare$ap.box, 0), ann=FALSE, xaxt="n", type='l')
    plot( compare$timestamp, pmin( compare$p - compare$ap.box, 0), ann=FALSE, xaxt="n", type='l')
    title(paste('Step3 : Post processing 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Removed / Added / After / Before", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
    
  }
  
  names(post.processing.data)[names(post.processing.data)=='p'] <- 'ap.box'
  names(post.processing.data)[names(post.processing.data)=='q'] <- 'rp.box'
  return(post.processing.data)
}

