generate.PatternScan.meta_1Hz <- function (data
                                           , genResolution.n = 1000
                                           , genEffSize.n = 12
                                           , staPeriodicity.p = 0.2
                                           , endEffSlot.p = 0.1
                                           , endConsistency.n = 1.0
                                           , min_sigmag.active.n    = 20
                                           , min_sigmag.reactive.n  = 3
                                           , s.periodsig.lossrate.p = 0.8 # more than this loss would not be considered in finding the start groups. 
                                           , clustering.method = 3
                                           , debug.mode = F
                                           , main_type = "reactive"
                                           , sub_type = "active"
){
  
  result1 <- generate.PatternScan.meta_summit_1Hz_new( data = data
                                                       , genResolution.n=genResolution.n
                                                       , genEffSize.n=genEffSize.n
                                                       , staPeriodicity.p=staPeriodicity.p
                                                       , endEffSlot.p=endEffSlot.p
                                                       , endConsistency.n=endConsistency.n
                                                       , min_sigmag.active.n=min_sigmag.active.n
                                                       , min_sigmag.reactive.n=min_sigmag.reactive.n
                                                       , min_lossrate.p=s.periodsig.lossrate.p
                                                       , sign_flag.b=1
                                                       , clustering.method=clustering.method
                                                       , debug.mode=debug.mode
                                                       , main_type=main_type
                                                       , sub_type=sub_type) 
  
  if( main_type != sub_type ){
    rev_data <- data
    rev_data$reactive_power <- (-rev_data$reactive_power)
    result2 <- generate.PatternScan.meta_summit_1Hz_new( data = rev_data
                                                         , genResolution.n=genResolution.n
                                                         , genEffSize.n=genEffSize.n
                                                         , staPeriodicity.p=staPeriodicity.p
                                                         , endEffSlot.p=endEffSlot.p
                                                         , endConsistency.n=endConsistency.n
                                                         , min_sigmag.active.n=min_sigmag.active.n
                                                         , min_sigmag.reactive.n=min_sigmag.reactive.n
                                                         , min_lossrate.p=s.periodsig.lossrate.p
                                                         , sign_flag.b=-1
                                                         , clustering.method=clustering.method
                                                         , debug.mode=debug.mode
                                                         , main_type=main_type
                                                         , sub_type=sub_type) 
    
    return(append(result1, result2))
  }else{
    return(result1)
  }
}




meta2PatternScan_1Hz <- function (data, meta.json, extension.p, flexibility.p, consistency.n, postprocessing = T, debug.mode, filename = "meta2pattern" ) {
  
  if(meta.json[['shape_type']] != 'pattern_scan') stop("Meta is not correct : meta2PatternScan_1Hz")
  
  if( debug.mode ) png( paste0(filename,"-%d.png"), width=3000)
  
  result <- meta2PatternScan_summit_1Hz_new(data, meta.json, extension.p, flexibility.p, consistency.n, postprocessing=postprocessing, debug.mode = debug.mode)
  
  if( debug.mode ) dev.off()
  
  return(result)  
}

meta2PatternScan_summit_1Hz_new <- function (data, meta, extension.p, flexibility.p, consistency.n, postprocessing, debug.mode, minsupportRatio = .5) {
  
  if( meta$supportRatio < minsupportRatio ) postprocessing <- F
  
  # parameter setting statement
  min_sigmag.active.n   <- 20
  min_sigmag.reactive.n <- 3
  
  # pattern matching
  InfoExt.p <- extension.p
  m.flex.p  <- flexibility.p
  r.edge    <- meta$rising_edge
  f.edge    <- meta$falling_edge
  
  s.rp_min <- r.edge['rp_min'] *(1-InfoExt.p*InfoExt.p*sign(r.edge['rp_min']))
  s.rp_max <- r.edge['rp_max'] *(1+InfoExt.p*InfoExt.p*sign(r.edge['rp_min']))
  s.rp_med <- r.edge['rp_med']
  
  s.ap_min <- r.edge['ap_min'] *(1-InfoExt.p*InfoExt.p*sign(r.edge['ap_min']))
  s.ap_max <- r.edge['ap_max'] *(1+InfoExt.p*InfoExt.p*sign(r.edge['ap_max']))
  s.ap_med <- r.edge['ap_med']
  
  s.t_min <- r.edge['min.t'] *(1-InfoExt.p)
  s.t_med <- r.edge['med.t']
  s.t_min.med.rate <- r.edge['min.med.rate']
  
  e.rp_min <- f.edge['rp_min'] *(1-InfoExt.p*sign(f.edge['rp_min']))
  e.rp_max <- f.edge['rp_max'] *(1+InfoExt.p*sign(f.edge['rp_max']))
  
  e.ap_min <- f.edge['ap_min'] *(1-InfoExt.p*sign(f.edge['ap_min']))
  e.ap_max <- f.edge['ap_max'] *(1+InfoExt.p*sign(f.edge['ap_max']))
  e.ap_med <- f.edge['ap_med']
  
  e.eff_t_min <- f.edge['EffTimeOn.min'] *(1-InfoExt.p)
  e.eff_t_med <- f.edge['EffTimeOn.med']
  
  e.eff_ap_med <- f.edge['EffAP_Drop.med']
  e.eff_ap_sd  <- f.edge['EffAP_Drop.med']
  
  e.eff_rp_med <- f.edge['EffRP_Drop.med']
  e.eff_rp_sd  <- f.edge['EffRP_Drop.sd']
  
  summit_flag <- meta$rising_edge['summit_flag']
  main_type <- switch( summit_flag, '1' = 'reactive', '-1' = '-reactive', '2' = 'active'   )
  sub_type  <- switch( summit_flag, '1' = 'active',   '-1' = 'active',    '2' = 'reactive' )
  
  # search start points
  s.pattern <- DetectPattern_1Hz_new( data = data
                                      , position = "start"
                                      , main_type = main_type
                                      , sub_type  = sub_type
                                      , debug.mode = F )
                                    
  s.pattern <- subset( s.pattern, (delta >= s.rp_min) & (delta <= s.rp_max) 
                       & (sub.delta >= s.ap_min) & (sub.delta <= s.ap_max))
  s.pattern <- s.pattern[order(s.pattern$start.timestamp),] 
  
  # examine s.pattern
  tmp_flag <- 1
  while (tmp_flag == 1) {
    
    if (min(diff(as.numeric(s.pattern$start.timestamp))) < s.t_min) {
      
      t.idx <- which.min(diff(as.numeric(s.pattern$start.timestamp)))
      
      likelihood1 <- sqrt( ( (s.rp_med - s.pattern$delta[t.idx] )/s.rp_med )^2 +
                             ( (s.ap_med - s.pattern$sub.delta[t.idx] )/s.ap_med )^2 )   
      
      likelihood2 <- sqrt( ( (s.rp_med - s.pattern$delta[t.idx+1] )/s.rp_med )^2 +
                             ( (s.ap_med - s.pattern$sub.delta[t.idx+1] )/s.ap_med )^2 )
      
      if (likelihood1 > likelihood2) {
        s.pattern <- s.pattern[-t.idx,]
      } else if (likelihood1 < likelihood2) {
        s.pattern <- s.pattern[-(t.idx+1),]
      } else {
        print("A start point has been removed randomly")
        s.pattern <- s.pattern[-t.idx,]
      }
    } else {tmp_flag <- 0}
  }
  
  # search end points
  e.pattern <- DetectPattern_1Hz_new( data = data
                                      , position = "end"
                                      , main_type = "reactive"
                                      , sub_type = "active"
                                      , c_factor = consistency.n
                                      , debug.mode = F )
  
  e.pattern <- subset( e.pattern, -delta >= min_sigmag.reactive.n & -sub.delta >= min_sigmag.active.n )
  
  e.pattern <- subset(e.pattern, (h2 >= e.rp_min) & (h2 <= e.rp_max) & 
                        (sub.delta >= e.ap_min) & (sub.delta <= e.ap_max) )
  e.pattern <- e.pattern[ order(e.pattern$start.timestamp),]
  
  if(debug.mode){
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
    if( nrow(s.pattern) > 0 ) abline( v=s.pattern$start.timestamp, col='red'  ) 
    if( nrow(e.pattern) > 0 ) abline( v=e.pattern$start.timestamp, col='blue' ) 
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  
  # edge selection
  t.duration <- list()
  chosen.ends <- list()
  if (nrow(s.pattern) != 0 && nrow(e.pattern) != 0) {
    
    e.pattern          <- subset( e.pattern, start.timestamp >= min(s.pattern$end.timestamp) )
    slot_number        <- findInterval( e.pattern$start.timestamp, vec=s.pattern$end.timestamp )
    e.pattern$duration <- as.numeric(e.pattern$start.timestamp) - as.numeric(s.pattern$end.timestamp[ slot_number ])
    
    flexibility <- e.eff_t_med * m.flex.p
    e.pattern <- subset( e.pattern, (e.eff_t_min < duration) & (duration < (e.eff_t_med + flexibility)))
    
    slot_number      <- findInterval( e.pattern$start.timestamp, vec=s.pattern$end.timestamp )
    slot_number.rle  <- rle(slot_number)
    slot_Property    <- slot_number.rle$lengths
    slot_ones        <- length( which( slot_Property == 1 ) ) 
    
    if( any(slot_Property != 0) ){
      matched.ends <- lapply( 1:length(s.pattern$end.timestamp), function(x){
        e.pattern.sub <- e.pattern[ slot_number == x, ]
        if( nrow(e.pattern.sub) == 0 ) return(NULL)
        calc_likelihood <- sqrt( ( (e.eff_rp_med - as.numeric(e.pattern.sub$h2) )/e.eff_rp_med )^2 +
                                   ( (e.eff_ap_med  - as.numeric(e.pattern.sub$sub.delta))/e.eff_ap_med )^2 )
        e.pattern.sub[ which.min(calc_likelihood), ]
      })
      tmp_cnt <- length(which(sapply( matched.ends, is.null )))
    }
    resultant.info <- data.frame( start.timestamp = s.pattern$start.timestamp[which(!sapply( matched.ends, is.null ))], 
                                  end.timestamp = ldply(matched.ends)$end.timestamp )
    if( any(sapply( matched.ends, is.null )) ){
      resultant.info <- rbind( resultant.info, 
                               data.frame( start.timestamp = s.pattern$start.timestamp[which(sapply( matched.ends, is.null ))], 
                                           end.timestamp = s.pattern$start.timestamp[which(sapply( matched.ends, is.null ))] + e.eff_t_med))
    }
    resultant.info <- resultant.info[order(resultant.info$start.timestamp),]
    cat("detection rate for appliance: ", (nrow(resultant.info) - tmp_cnt)/nrow(s.pattern) , "\n")
    
    if(debug.mode){
      par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
      plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
      if( nrow(resultant.info) > 0 ) abline( v=resultant.info$start.timestamp, col='red'  ) 
      if( nrow(resultant.info) > 0 ) abline( v=resultant.info$end.timestamp, col='blue' ) 
      plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
      title(paste('Step2 : Edge Matching 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
      mtext("Timestamp", 1, 1, outer=TRUE)
      mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
      par(mfrow=c(1,1))
    }
    
    # dumb restoration
    if ( postprocessing ) {
      
      tmp_flag <- 1
      while (tmp_flag == 1) {
        
        # restoration from start point
        t.value       <- max(diff(as.numeric(resultant.info$start.timestamp)))
        t.value_start <- diff(c( as.numeric(min(data$timestamp)), as.numeric(min(resultant.info$start.timestamp)) ))
        t.value_end   <- diff(c( as.numeric(max(resultant.info$start.timestamp)), as.numeric(max(data$timestamp)) ))
        
        if (t.value >= 1.5*s.t_med) {
          
          t.idx <- which.max(diff(as.numeric(resultant.info$start.timestamp)))
          t.slot <- t.value/round(t.value/s.t_med)
          for (f_idx in 1:(round(t.value/s.t_med)-1)) {
            
            tmp_start.t <- as.numeric(resultant.info$start.timestamp[t.idx]) + t.slot*f_idx
            str.t <- as.POSIXct(tmp_start.t, origin = "1970-01-01", tz = "Asia/Seoul")
            end.t <- str.t + e.eff_t_med
            resultant.info <- rbind( resultant.info, 
                                     data.frame( start.timestamp = str.t, end.timestamp = end.t ) )
            
          }
          resultant.info <- resultant.info[order(as.numeric(resultant.info$start.timestamp)),] 
        } else {
          tmp_flag <- 0
        }
        
        
        if (t.value_start >= 1.5*s.t_med){
          t.slot <- t.value_start/round(t.value_start/s.t_med)
          for (f_idx in 1:(round(t.value_start/s.t_med)-1)) {
            tmp_start.t <- as.numeric(min(resultant.info$start.timestamp)) - t.slot*f_idx
            
            str.t <- as.POSIXct(tmp_start.t, origin = "1970-01-01", tz = "Asia/Seoul")
            end.t <- str.t + e.eff_t_med
            resultant.info <- rbind( resultant.info, 
                                     data.frame( start.timestamp = str.t, end.timestamp = end.t ) )
          }
        }
        
        if (t.value_end >= 1.5*s.t_med){
          t.slot <- t.value_end/round(t.value_end/s.t_med)
          for (f_idx in 1:(round(t.value_end/s.t_med)-1)) {
            tmp_start.t <- as.numeric(max(resultant.info$start.timestamp)) + t.slot*f_idx
            
            str.t <- as.POSIXct(tmp_start.t, origin = "1970-01-01", tz = "Asia/Seoul")
            end.t <- str.t + e.eff_t_med
            resultant.info <- rbind( resultant.info, 
                                     data.frame( start.timestamp = str.t, end.timestamp = end.t ) )
          }
        }
        
        resultant.info <- resultant.info[order(as.numeric(resultant.info$start.timestamp)),] 
        
      }
      
    }
    
    # print(resultant.info)
    
  } else {
    
    print("Error: no start or end point!")
    return(NULL)
  }
  
  if (nrow(resultant.info) != 0) {
    
    # log extraction
    p.pad <-( s.ap_med + abs(e.eff_ap_med) ) / 2
    q.pad <-( s.rp_med + abs(e.eff_rp_med) ) / 2
    
    answer.log <- data.frame( timestamp = data$timestamp, p= 0, q= 0)
    if(postprocessing){
      ts.diff <- round(as.numeric( diff(data$timestamp), units = 'secs' ))
      missingData <- lapply( which(ts.diff > 1), 
                             function(i) data.frame( timestamp = seq( data$timestamp[i], 
                                                                      data$timestamp[i+1],
                                                                      length.out=ts.diff[i]-1), 
                                                     p = 0, q = 0) )
      answer.log <- unique(rbind( answer.log, rbind.fill(missingData) ))
      answer.log <- answer.log[ order(answer.log$timestamp), ]
    }
    
    signalOn <- unlist( mapply( function(str,end){which(answer.log$timestamp %within% interval(str, end))}, 
                                resultant.info$start.timestamp, resultant.info$end.timestamp ) )
    answer.log$p[signalOn] <- p.pad
    answer.log$q[signalOn] <- q.pad
    
    
    if(debug.mode){
      par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
      plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
      plot(x=answer.log$timestamp, y=answer.log$p, ann=FALSE, xaxt="n", type='l')
      title(paste('Step3 : Post processing 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
      mtext("Timestamp", 1, 1, outer=TRUE)
      mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
      par(mfrow=c(1,1))
    }
    return(answer.log)
  } else {
    stop("Error: no end point!")
  }
}



edgeCluster <- function(data, xColName, yColName, resolution, d_factor, gp_size, 
                        vec, clusterMethod = 1, debug.mode, eps=.15, 
                        xBuffer_mag = 0, yBuffer_mag = 0){
  
  findGaps <- function(x,n){
    x <- sort(x)
    x.diff <- data.frame( val = diff(x), idx = 1:(length(x)-1) )
    wall.idx <- x.diff$idx[ order( x.diff$val, decreasing=T ) ][1:n] # choose first N gaps
    result <- sapply( wall.idx, function(i) mean(x[i+c(0,1)]))
    return( sort(result) )
  }
  
  divideByGaps <- function(x,n){
    breaks   <- findGaps( x, n )
    division <- cut( x, c(min(x), breaks, max(x)) ,labels= FALSE, include.lowest=T ) -1
    return(list(division,breaks))
  }
  
  if( missing(vec) ){
    gp.summary <- function(data){
      result <- cell.summary( cell=data, xColName=xColName, yColName=yColName )
      result$d.rate <- (result$min.med.rate2 >= d_factor) & result$leftCell.exist & result$bottomCell.exist
      return(result)
    }
  }else{
    data       <- subset( data, start.timestamp >= min(vec) )
    vec        <- sort(vec)
    gp_size    <- length(vec) * d_factor
    gp.summary <- function(data){
      result <- cell.summary2( cell=data, vec=vec, xColName=xColName, yColName=yColName )
      result$d.rate <- (result$effective.rate >= d_factor) & result$leftCell.exist & result$bottomCell.exist
      return(result)
    } 
  }
  
  if( clusterMethod == 1 ){
    logSummary <- list()
    logDetail  <- list()
    if(debug.mode) plot( data[,xColName], data[,yColName] )
    
    r_idx <- 1
    for( iter in 1:resolution ){
      
      if ( nrow(data) +1 <= r_idx ){
        print("Resolution has been reduced..")
        break
      }
      
      xDivision <- divideByGaps( data[,xColName], r_idx )
      yDivision <- divideByGaps( data[,yColName], r_idx )
      data      <- cbind( data, data.frame( g1 = xDivision[[1]], 
                                            g2 = yDivision[[1]]) )
      data$leftCell.exist   <- (data$g1 != 0)
      data$bottomCell.exist <- (data$g2 != 0)
      
      if(debug.mode){
        abline( v = xDivision[[2]], col='red'  )
        abline( h = yDivision[[2]], col='blue' )
      } 
      
      data <- ddply( data, .(g1, g2), transform, len = length(g1))
      small.group <- unique(data[ data$len < gp_size, c('g1','g2') ])
      r_idx <- max( r_idx - nrow( small.group ), 1 )
      
      ### exclude the info from isolated elements
      data <- subset( data, len >= gp_size ) 
      
      ### effective group (conservative search by default)
      if (nrow(data) > 0){
        group.info <- ddply( data, .(g1, g2), gp.summary )
        eff_group  <- group.info[ group.info$d.rate, ]
        
        if( nrow(eff_group) > 0 ){
          for( i in 1:nrow(eff_group) ){
            logDetail[[length(logDetail) +1]] <- subset( data, (g1 == eff_group$g1[i]) & (g2 == eff_group$g2[i]))
            data <- subset( data, !((g1 == eff_group$g1[i]) & (g2 == eff_group$g2[i])))
            if(debug.mode){
              tmp <- logDetail[[length(logDetail)]]
              points( tmp[,xColName], tmp[,yColName], col='green')
            }
          }
          logSummary <- rbind( logSummary, eff_group )
          r_idx <- max( r_idx - nrow(eff_group) + 1, 1 ) # index increment
        }else{
          r_idx <- r_idx + 1 # index increment
        }
      }
      data$g1  <- NULL
      data$g2  <- NULL
      data$len <- NULL
    }
    
    # return information
    if(NROW(logSummary)>0) logSummary$candidate.idx <- 1:nrow(logSummary)
    logDetail[[length(logDetail) +1]] <- logSummary
    return(logDetail)
    
  }else if( clusterMethod == 2 ){
    data      <- mutate( data, cluster = 1, leftCell.exist = F, bottomCell.exist = F )
    data$idx  <- seq_along(data[,xColName])
    data.list <- list(data)
    
    dbscanResult <- fpc::dbscanCBI( as.matrix(data[,c(xColName,yColName)])
                                    , eps = eps
                                    , MinPts=gp_size
                                    , scale = T
                                    , showplot = debug.mode)
    
    is.Periodic <- sapply( dbscanResult$clusterlist, 
                           function(log){ cell.info <- gp.summary(data[log, ]) 
                                          return(cell.info$d.rate)})
    if(all(!is.Periodic)) return(list())
    
    logDetail <- ldply( dbscanResult$clusterlist[is.Periodic], 
                        function(log) gp.summary(data[log,])) 
    data.subset <- ldply( dbscanResult$clusterlist[is.Periodic], 
                          function(log) data[log,]) 
    
    if(debug.mode & length(which(is.Periodic)) > 0){
      plot( data[,xColName], data[,yColName] )
      points( data.subset[,xColName], data.subset[,yColName], col='green')
    } 
    return( list(data.subset, logDetail))
    
  }else if( clusterMethod == 3 ){
    
    divideCell <- function( cell ){
      if( nrow(cell) <= 2 ) return(list(cell))
      if(!all(cell$converged)){
        cell.pca <- stats::prcomp( cell[,c(xColName, yColName)], center = TRUE, scale. = TRUE)
        d.level  <- cell.pca$rotation[2,1]/cell.pca$rotation[1,1]
        
        if( abs(d.level) > tan(x=pi/6) ){
          xDivision <- divideByGaps(cell[,xColName], 1)
          cell$xCluster <- xDivision[[1]]
          #cell$leftCell.exist <- length(sign(cell[,xColName])) == 1
          if( xBuffer_mag == 0 ){
            cell$leftCell.exist <- (cell$xCluster != 0)
          }else{
            for( i in unique(cell$xCluster) ){
              cell$leftCell.exist[ cell$xCluster == i ] <- 
                all( round(abs(cell[ cell$xCluster == i, xColName])) > xBuffer_mag )
            }
          }
          #cell$leftCell.exist <- (cell$xCluster != 0)
          if(debug.mode) segments( y0 = min(cell[,yColName]), y1 = max(cell[,yColName]),
                                   x0 = xDivision[[2]], x1 = xDivision[[2]], col='red'  )
        }
        if( abs(d.level) < tan(x=pi/3) ){
          yDivision <- divideByGaps(cell[,yColName], 1)
          cell$yCluster <- yDivision[[1]]
          #cell$bottomCell.exist <- length(sign(cell[,yColName])) == 1
          if( yBuffer_mag == 0 ){
            cell$bottomCell.exist <- (cell$yCluster != 0)
          }else{
            for( i in unique(cell$yCluster) ){
              cell$bottomCell.exist[ cell$yCluster == i ] <- 
                all( round(abs(cell[ cell$yCluster == i, yColName ])) > yBuffer_mag )
            }
          }
          #cell$bottomCell.exist <- (cell$yCluster != 0)
          if(debug.mode) segments( x0 = min(cell[,xColName]), x1 = max(cell[,xColName]),
                                   y0 = yDivision[[2]], y1 = yDivision[[2]], col='blue'  )
        }
      }
      return(dlply( cell, .(xCluster, yCluster)))
    }
    
    checkConvergence <- function(ith){
      cell           <- data.list[[ith]]
      cell$cluster   <- ith
      cell$converged <- gp.summary(cell)$d.rate
      return(cell)
    }
    
    data <- mutate( data, cluster = 1, xCluster = 1, yCluster = 1, leftCell.exist = F, bottomCell.exist = F )
    data$converged <- F
    data.list <- list(data)
    if(debug.mode) plot( data[,xColName], data[,yColName], xlab = xColName, ylab = yColName )
    
    for( iter in 1:resolution ){
      data.list <- lapply( data.list, divideCell )
      data.list <- unlist( data.list, recursive = F )
      data.list[ sapply( data.list, nrow ) <= gp_size ] <- NULL
      data.list <- lapply( seq_along(data.list), checkConvergence )
      converge <- sapply( data.list, function(cell) unique(cell$converged[!is.na(cell$converge)]))
      if(all(converge)) break
    }
    
    if( all(!converge) ) return(NULL)
    if( debug.mode ){
      tmp <- rbind.fill(data.list[ converge ])
      points( tmp[,xColName], tmp[,yColName], col=c(factor(tmp$cluster))+1 )
    } 
    logDetail  <- data.list[ converge ]
    logSummary <- ldply( data.list[ converge ], gp.summary )
    if(nrow(logSummary)>0){
      logSummary$candidate.idx <- 1:nrow(logSummary)
    } 
    logDetail[[length(logDetail)+1]] <- logSummary
    return( logDetail )
  }
}

generate.PatternScan.meta_summit_1Hz_new <- function ( data
                                                       , genResolution.n
                                                       , genEffSize.n
                                                       , staPeriodicity.p
                                                       , endEffSlot.p
                                                       , endConsistency.n
                                                       , min_sigmag.active.n
                                                       , min_sigmag.reactive.n
                                                       , min_lossrate.p
                                                       , sign_flag.b
                                                       , clustering.method = 1
                                                       , debug.mode = F
                                                       , main_type 
                                                       , sub_type  ) {
  
  # parameters
  s.max_resolution.n <- genResolution.n
  e.max_resolution.n <- genResolution.n
  s.min_gsize.n      <- genEffSize.n
  g.eff_sample.n     <- genEffSize.n 
  s.p_factor.p       <- staPeriodicity.p
  e.eff_rate.p       <- endEffSlot.p      # select small -> big
  e.eff_rate_ext.p   <- e.eff_rate.p / 2
  
  # search starting points
  if( main_type == 'reactive' & sub_type == 'active'){
    use_posRpAp <- T; use_negRpAp <- F; use_Ap <- F
  }else if( main_type == '-reactive' & sub_type == 'active'){
    use_posRpAp <- F; use_negRpAp <- T; use_Ap <- F
  }else if( main_type == 'active' ){
    use_posRpAp <- F; use_negRpAp <- F; use_Ap <- T
  }
  s.pattern <- DetectPattern_1Hz_extended( data
                                           , min_mag_ap = min_sigmag.active.n
                                           , min_mag_rp = min_sigmag.reactive.n
                                           , position = "start"
                                           , main_type = main_type 
                                           , sub_type  = sub_type
                                           , debug.mode = F
                                           , use_posRpAp = use_posRpAp
                                           , use_negRpAp = use_negRpAp
                                           , use_Ap = use_Ap )
  
  
  if( nrow(s.pattern) == 0 ){
    print("Start point detection failed")
    return(list())
  }
  if( debug.mode ){
    plot( data$active_power, type='l')
    abline( v = s.pattern$start.idx, col='red')
    title(paste('Step1 : rising edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
  }
  
  s.info <- list()
  s.info <- edgeCluster( s.pattern , 
                         , xColName = ifelse( use_Ap, 'ap.delta', 'rp.delta' )
                         , yColName = ifelse( use_Ap, 'ap.delta', 'ap.delta' )
                         , resolution = s.max_resolution.n
                         , d_factor = s.p_factor.p
                         , gp_size = s.min_gsize.n
                         , clusterMethod = clustering.method
                         , debug.mode = debug.mode
                         , xBuffer_mag=min_sigmag.reactive.n
                         , yBuffer_mag=min_sigmag.active.n )
  
  if (length(s.info) == 0 || NROW(s.list <- getLast(s.info)) == 0){
    print("Start point detection failed")
    return(list())
  }
  print("Detected groups for starting points")
  print(s.list)
  
  e.pattern <- DetectPattern_1Hz_extended( data
                                           , min_mag_ap = min_sigmag.active.n
                                           , min_mag_rp = min_sigmag.reactive.n
                                           , position = "end"
                                           , main_type = main_type 
                                           , sub_type  = sub_type
                                           , c_factor = endConsistency.n
                                           , debug.mode = F
                                           , use_posRpAp = use_posRpAp
                                           , use_negRpAp = use_negRpAp
                                           , use_Ap = use_Ap)
  
                                      
  if( nrow(e.pattern) == 0 ){
    print("End point detection failed")
    return(list())
  }
    
  if( debug.mode ){
    plot( data$active_power, type='l')
    abline( v = e.pattern$start.idx, col='blue')
    title(paste('Step1 : falling edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
  }

  chosen.e.idx <- c()
  chosen.e.data.list <- list()
  chosen.e.data.summary <- list()
  
  # detect corresponding end points for each start group
  for(s_idx in 1:nrow(s.list)){
    
    if (s.list$lost.sig.rate[s_idx] <= (min_lossrate.p * 100) ) {
      
      e.info <- edgeCluster( data = e.pattern
                             , xColName = ifelse( use_Ap, 'ap.delta', 'rp.delta' )
                             , yColName = ifelse( use_Ap, 'ap.delta', 'ap.delta' )
                             , resolution = e.max_resolution.n
                             , d_factor = e.eff_rate.p
                             , vec = s.info[[s_idx]]$end.timestamp
                             , clusterMethod = clustering.method
                             , debug.mode = debug.mode )
      
      need.more.search <- F
      if (length(e.info) != 0) {
        e.list <- subset(getLast(e.info), metric.info < 1)
        if (nrow(e.list) == 0) need.more.search <- T
      }else need.more.search <- T
      
      if( need.more.search ){
        print("Reducing 'effctive slot rate' to find end points again")
        tmp_e.eff_rate.p <- e.eff_rate.p - e.eff_rate_ext.p
        e.info <- edgeCluster( data = e.pattern
                               , xColName = ifelse( use_Ap, 'ap.delta', 'rp.delta' )
                               , yColName = ifelse( use_Ap, 'ap.delta', 'ap.delta' )
                               , resolution = e.max_resolution.n
                               , d_factor = tmp_e.eff_rate.p
                               , vec = s.info[[s_idx]]$end.timestamp
                               , clusterMethod = clustering.method
                               , debug.mode = debug.mode )
      }
      
      if( NROW(e.info) > 0 ){
        e.list <- subset(getLast(e.info), metric.info < 1)
      }else{
        e.list <- NULL
      }
      
      if (!is.null(e.list) && nrow(e.list) != 0){
        cat("End point candidate list for appliance ", s_idx, "\n")
        print(e.list)
        
        chooseInterval.1 <- findInterval( s.list$sum[s_idx], c(0,50,100) )
        if( chooseInterval.1 == 1 ){ # 0 <= nrow(s.info[[s_idx]] < 50
          chooseInterval.2 <- findInterval( e.list$metric.info, c(0, .9))
        }else if( chooseInterval.1 == 2 ){ #  50 <= nrow(s.info[[s_idx]] < 100
          chooseInterval.2 <- findInterval( e.list$metric.info, c(0, .6, .9))        
        }else if( chooseInterval.1 == 3 ){ # 100 <= nrow(s.info[[s_idx]]
          chooseInterval.2 <- findInterval( e.list$metric.info, c(0, .2, .6, .9))
        }
        
        e.list$metricInterval <- chooseInterval.2
        
        chosen.e.data <- data.frame()
        chosen.e.list <- data.frame()
        for( i in sort(unique(chooseInterval.2)) ){
          ex.e.list <- e.list[ e.list$metricInterval == i, ]
          if( nrow(ex.e.list) > 0 ){
            tmp_compute <- which.min(sqrt( ( (ex.e.list$eff.PowerDrop.med - s.list$med.sub.d[s_idx])/s.list$med.sub.d[s_idx] )^2 +
                                             ((ex.e.list$eff.rpDrop.med - s.list$med.d[s_idx])/s.list$med.d[s_idx] )^2     ))
            
            chosen.e.data <- e.info[[ex.e.list$candidate.idx[tmp_compute]]]
            chosen.e.list <- ex.e.list[tmp_compute, ]
            break
          }
        }
        if ( nrow(chosen.e.data) == 0 ){
          print("End point detection failed")
        }else{
          if (nrow(chosen.e.data) != 0){
            chosen.e.idx <- c(chosen.e.idx, s_idx)
            chosen.e.data.list[[length(chosen.e.data.list)+1]]       <- chosen.e.data
            chosen.e.data.summary[[length(chosen.e.data.summary)+1]] <- chosen.e.list
            #e.pattern <- e.pattern[ !(e.pattern$start.idx %in% chosen.e.data$start.idx), ]
            if( nrow(e.pattern) == 0 ) break
            #r.edge <- s.list[s_idx]
            #f.edge <- chosen.e.list
            
          }
        }       
      }else{
        print("End point detection failed") 
      }
    }
  }
  
  if (length(chosen.e.idx) != 0) {
    
    str.t <- as.character(min(data$timestamp))
    end.t <- as.character(max(data$timestamp))
    
    # extract meta data
    json.result <- list()
    for (m_idx in 1:length(chosen.e.idx)) {
      
      
      s.data.log <- s.info[[chosen.e.idx[m_idx]]]
      e.data.log <- chosen.e.data.list[[m_idx]]
      e.data.summary <- chosen.e.data.summary[[m_idx]]
      
      supportRatio <- as.numeric( diff( c( min(s.data.log$start.timestamp),
                                           max(e.data.log$end.timestamp))), units='secs')
      supportRatio <- supportRatio / as.numeric(diff(range(data$timestamp)),units='secs')
      
      # start data info
      #s.data.log.melt <- melt( s.data.log, measure.vars=c('rp.delta', 'ap.delta'))
      #data.file = '/home/hslee/test_data/meta_test_data/1442415600000-1443625199999.10002822_4646_3.csv.gz', find.heavy=F, find.pattern=F, find.pattern_extend=T,

      tmp_start.json <- c( summit_flag = sign_flag.b
                           , rp_min = min(s.data.log$rp.delta)
                           , rp_max = max(s.data.log$rp.delta)
                           , rp_med = median(s.data.log$rp.delta)
                           , ap_min = min(s.data.log$ap.delta)
                           , ap_max = max(s.data.log$ap.delta)
                           , ap_med = median(s.data.log$ap.delta)
                           , sample.num = s.list$sum[chosen.e.idx[m_idx]]
                           , min.t = s.list$min.t[chosen.e.idx[m_idx]]
                           , med.t = s.list$med.t[chosen.e.idx[m_idx]]
                           , min.med.rate = s.list$min.med.rate[chosen.e.idx[m_idx]]
                           , lost.sample.num = s.list$lost.sig.num[chosen.e.idx[m_idx]]
                           , lost.sample.rate = s.list$lost.sig.rate[chosen.e.idx[m_idx]])
      
      # end data info
      tmp_end.json <- c( rp_min         = min(e.data.log$rp.delta), 
                         rp_max         = max(e.data.log$rp.delta),
                         rp_med         = median(e.data.log$rp.delta),
                         ap_min         = min(e.data.log$ap.delta), 
                         ap_max         = max(e.data.log$ap.delta),
                         ap_med         = median(e.data.log$ap.delta), 
                         slotNum.zero   = e.data.summary$slot_zeros, 
                         slotNum.one    = e.data.summary$slot_ones, 
                         ZerotoOneratio = e.data.summary$metric.info, 
                         EffTimeOn.med  = e.data.summary$eff.TimeOn.med, 
                         EffTimeOn.min  = e.data.summary$eff.TimeOn.min, 
                         EffTimeOn.max  = e.data.summary$eff.TimeOn.max,
                         EffTimeOn.sd   = e.data.summary$eff.TimeOn.sd, 
                         EffAP_Drop.med = e.data.summary$eff.ap.delta.med, 
                         EffAP_Drop.min = e.data.summary$eff.ap.delta.min, 
                         EffAP_Drop.max = e.data.summary$eff.ap.delta.max, 
                         EffAP_Drop.sd  = e.data.summary$eff.ap.delta.sd, 
                         EffRP_Drop.med = ifelse( is.null(e.data.summary$eff.rp.delta.med), 0, e.data.summary$eff.rp.delta.med ), 
                         EffRP_Drop.min = ifelse( is.null(e.data.summary$eff.rp.delta.min), 0, e.data.summary$eff.rp.delta.min ),
                         EffRP_Drop.max = ifelse( is.null(e.data.summary$eff.rp.delta.max), 0, e.data.summary$eff.rp.delta.max ), 
                         EffRP_Drop.sd  = ifelse( is.null(e.data.summary$eff.rp.delta.sd),  0, e.data.summary$eff.rp.delta.sd  )) 
      
      parameters <- c('genResolution.n'=genResolution.n
                      , 'genEffSize.n'=genEffSize.n
                      , 'staPeriodicity.p'=staPeriodicity.p
                      , 'endEffSlot.p'=endEffSlot.p
                      , 'endConsistency.n'=endConsistency.n
                      , 'min_sigmag.active.n'=min_sigmag.active.n
                      , 'min_sigmag.reactive.n'=min_sigmag.reactive.n
                      , 'min_lossrate.p'=min_lossrate.p
                      , 'clustering.method'=clustering.method
                      , 'main_type'=main_type 
                      , 'sub_type'=sub_type )
      
      json.result[[length(json.result)+1]] <- list(`meta-version` = 1, 
                                                   shape_type = "pattern_scan", 
                                                   rising_edge = tmp_start.json, 
                                                   falling_edge = tmp_end.json, 
                                                   supportRatio = supportRatio, 
                                                   generation_info = list(data_used = list(start = str.t, end = end.t, sampling = 1), 
                                                                          computed = as.character(Sys.time())),
                                                   parameters= parameters)
      
    }
    
    return(json.result)
  } else {
    return(list())
  }
}



#============================================================================================  
# < supplementary functions >
#============================================================================================

get_speed_new <- function(timestamp, value, allowInfSlope = F){
  
  if(class(timestamp)[1] != 'numeric') timestamp <- as.numeric(timestamp)
  time.diff  <- diff(timestamp)
  value.diff <- diff(value)
  slope      <- value.diff / time.diff
  
  if( !allowInfSlope ){
    infinite.slope <- which(!is.finite(slope))
    if( length(infinite.slope) > 0 ){
      finite.slope   <- setdiff( 1:length(slope), infinite.slope )
      if(             1 %in% infinite.slope ) slope[1:min(finite.slope)] <- slope[min(infinite.slope)]
      if( length(slope) %in% infinite.slope ) slope[max(finite.slope):length(slope)] <- slope[max(finite.slope)]
      infinite.slope <- setdiff( infinite.slope, c(1,length(slope)) )
#       idx <-   sapply( infinite.slope, function(x) c( max(finite.slope[ finite.slope < x ]), 
#                                                        min(finite.slope[ finite.slope > x ]) ))
#       slope[ infinite.slope ] <- ( slope[idx[1,]] + slope[idx[2,]] )/2
      split.infSlope <- split( infinite.slope, cumsum(c(-1,diff(infinite.slope)) != 1) )
      idx      <- sapply( split.infSlope, function(x) c( head(x,1)-1, tail(x,1)+1) )
      idx.left <- unlist(mapply( function(minVal,duplicatedItem) rep(minVal,length(duplicatedItem)), 
                                 idx[1,], split.infSlope ))
      idx.rite <- unlist(mapply( function(maxVal,duplicatedItem) rep(maxVal,length(duplicatedItem)), 
                                 idx[2,], split.infSlope ))
      slope[ infinite.slope ] <- ( slope[idx.left] + slope[idx.rite] )/2
    }
  }
  if( length(which(!is.finite(slope))) > 0 ) stop("Algorithm is not correct : get_speed_new")
  return( slope )
}


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

stateCompression <- function( state ){
  
  stateChangePos <- which( head(state,-1) != tail(state,-1) )
  
  state.str <- c(0,stateChangePos) + 1
  state.end <- c(stateChangePos,length(state))
  state.val <- state[ state.str ]
  
  return( data.frame( str = state.str, end = state.end, val = state.val) )
}

compressedPatternMatching <- function( compressedState, targetPattern ){
  
  pattern.length <- length(targetPattern) 
  
  if( targetPattern[1] %in% c(1,2) ){ # rising edge
    
    if( any( targetPattern %in% c(3,4) & shift(targetPattern,-1) %in% c(1,2) ) ){
      # power가 증가 -> 감소 -> 증가로 정의되면 마지막 증가는 패턴에서 제외
      pattern.length <- max(which(targetPattern %in% c(3,4) & shift(targetPattern,-1) %in% c(1,2) ))
    }
    
    matchedPattern.str <- grepl.pattern( compressedState$val, targetPattern )
    matchedPattern.end <- matchedPattern.str + pattern.length - 1
    maxSlope.pos <- matchedPattern.str + which( targetPattern == 1 & shift(targetPattern,-1) %in% c(2,3)) - 1
    minSlope.pos <- matchedPattern.str + which( targetPattern == 3 & shift(targetPattern,-1) %in% c(4,1)) - 1 
    peak.pos <- matchedPattern.str + which( targetPattern %in% c(1,2) & shift(targetPattern,-1) %in% c(3,4)) - 1
    
    result <- data.frame(str = compressedState$str[ matchedPattern.str ], 
                         end = compressedState$end[ matchedPattern.end ])
    
    if( length(maxSlope.pos) != 0 ) result <- cbind( result, data.frame(maxSlope = compressedState$end[ maxSlope.pos ]))
    if( length(minSlope.pos) != 0 ) result <- cbind( result, data.frame(minSlope = compressedState$end[ minSlope.pos ]))
    if( length(peak.pos) == 0 ) peak.pos <- matchedPattern.end
    
    result <- cbind( result, data.frame(peak = compressedState$end[ peak.pos ] + 1)) 
    return( result )
    
  }else if( targetPattern[1] %in% c(3,4) ){ # falling edge 
    
    matchedPattern.str <- grepl.pattern( compressedState$val, targetPattern )
    matchedPattern.end <- matchedPattern.str + pattern.length - 1
    minSlope.pos <- matchedPattern.str + which( targetPattern == 3 & shift(targetPattern,-1) %in% c(4,1,2)) - 1 
    peak.pos <- matchedPattern.str 
    
    return( data.frame(str = compressedState$str[ matchedPattern.str ], 
                       end = compressedState$end[ matchedPattern.end ],
                       minSlope = compressedState$end[ minSlope.pos ],
                       peak = compressedState$str[ peak.pos ]))
  }
}

grepl.pattern <- function(united_state, targetPattern){
  
  n <- length(targetPattern)
  targetPattern.boolean <- sapply( 0:(n-1), function(i) shift( united_state == targetPattern[i+1],-i))
  return( which( rowSums( targetPattern.boolean ) == n ) )
}


DetectPattern_1Hz <- function ( data, position = c("start", "end")
                                , main_type = c("reactive", "active")
                                , sub_type  = c("active", "reactive")
                                , min_mag_main = 3
                                , min_mag_sub = 20
                                , c_factor = 0
                                , debug.mode = F ){
    
  mainLog <- switch( main_type, 
                     'reactive' = data.frame( timestamp = data$timestamp, value = data$reactive_power), 
                     'active'   = data.frame( timestamp = data$timestamp, value = data$active_power) )
  subLog <- switch( sub_type, 
                     'reactive' = data.frame( timestamp = data$timestamp, value = data$reactive_power), 
                     'active'   = data.frame( timestamp = data$timestamp, value = data$active_power) )
  
  # for debug
  if ( nrow(mainLog) != nrow(subLog) ){
    cat("Log mismatch between active and reactive power!", nrow(mainLog), " and ", nrow(subLog), "\n")
    stop("Error")
  }
  
  slope <- get_speed_new(mainLog$timestamp, mainLog$value) 
  state <- scan_state_new(slope) 
  united.state <- stateCompression( state )
  
  edgeType <- match.arg(position)
  if ( edgeType == "start"){
    
    rising.edge.pattern <- list( c(1,2,3,4), c(1,3,4), c(1,2,3,1), c(1,3,1) )
    resultant1 <- ldply( rising.edge.pattern, 
                        function( pattern ){ patternIdx <- compressedPatternMatching( united.state, pattern )
                                             get_patternFeature( mainLog, subLog, slope, patternIdx )} )
    
    rising.edge.pattern <- list( c(1,2), c(1,3) )
    resultant2 <- ldply( rising.edge.pattern, 
                        function( pattern ){ patternIdx <- compressedPatternMatching( united.state, pattern )
                                             get_patternFeature( mainLog, subLog, slope, patternIdx )} )
    resultant2 <- subset( resultant2, !(start.idx %in% resultant1$start.idx))
    if( nrow(resultant2) > 0 ) resultant1 <- rbind.fill( resultant1, resultant2 )
    
    resultant <- resultant1
    if( is.null(resultant) ) stop('There is no proper rising edge')
    resultant <- subset(resultant, (delta >= min_mag_main) & (sub.delta >= min_mag_sub))
    
    if(debug.mode){
      par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
      plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
      abline( v=resultant$start.timestamp, col='red' ) 
      plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
      title(paste('Step1 : rising edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
      mtext("Timestamp", 1, 1, outer=TRUE)
      mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
      par(mfrow=c(1,1))
    }
    return( resultant )
    
  }else if (edgeType == "end"){
    
    falling.edge.pattern <- list( c(3,4), c(3,1), c(3,2) )
    resultant <- ldply( falling.edge.pattern, 
                        function( pattern ){ patternIdx <- compressedPatternMatching( united.state, pattern )
                                             get_patternFeature( mainLog, subLog, slope, patternIdx, c_factor )} )
    if( is.null(resultant) ) stop('There is no proper falling edge')
    resultant <- subset(resultant, (-delta >= min_mag_main) & (-sub.delta >= min_mag_sub))
    
    if(debug.mode){
      par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
      plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
      abline( v=resultant$start.timestamp, col='blue' ) 
      plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
      title(paste('Step1 : falling edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
      mtext("Timestamp", 1, 1, outer=TRUE)
      mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
      par(mfrow=c(1,1))
    }
    
    resultant[,sapply(resultant, is.numeric)] <- abs(resultant[,sapply(resultant, is.numeric)])
    return( resultant )
  }  
}  

#============================================================================================
# h1 : peak - left (main)
# h2 : rite - peak (main)
# delta : rite - left = h1 + h2 (main) 
# sub.delta : rite - left (sub)
get_patternFeature <- function(data, sub.data, speed, patternIndex, c_factor, use_Ap = FALSE, useConsistencyFlag = TRUE){
  
  timestamp <- data$timestamp
  value     <- data$value
  sub.value <- sub.data$value
  
  patternFeature        <- patternIndex[,c('str','end')]
  patternFeature$len    <- patternIndex$end - patternIndex$str + 1
  patternFeature$str.t  <- timestamp[patternIndex$str  ]
  patternFeature$end.t  <- timestamp[patternIndex$end+1]
  
  if('maxSlope' %in% names(patternIndex)) patternFeature$max.slope <- speed[ patternIndex$maxSlope ]
  if('minSlope' %in% names(patternIndex)) patternFeature$min.slope <- speed[ patternIndex$minSlope ]
  
  if('peak' %in% names(patternIndex)){
    peak.value.set           <- value[ patternIndex$peak ]
    patternFeature$h1        <- peak.value.set - value[ patternIndex$str ]
    patternFeature$h2        <- value[ patternIndex$end+1 ] - peak.value.set
  }
  patternFeature$sub.delta <- sub.value[patternIndex$end+1] - sub.value[patternIndex$str]
  patternFeature$delta     <- value[patternIndex$end+1]     - value[patternIndex$str]
  
  if( !missing(c_factor) && c_factor >= 1 ){
    
    patternFeature$org_h2        <-  patternFeature$h2
    patternFeature$org_sub.delta <-  patternFeature$sub.delta
    
    patternFeature <- subset( patternFeature, (str > c_factor*len) & (end < length(value)-c_factor*len) )
    if( nrow(patternFeature) == 0 ) return(NULL)
    
    consistency_flag <- rep(1, nrow(patternFeature))
    if( !use_Ap ){
      consistency_flag[ which(patternFeature$org_h2 > 0 | patternFeature$org_sub.delta > 0) ] <- 0
    }else{
      consistency_flag[ which(patternFeature$org_h2 > 0) ] <- 0
    }
      
    for(c_idx in 1:c_factor) {
      
      left.str <- patternFeature$str - patternFeature$len * c_idx
      left.end <- patternFeature$str - patternFeature$len * (c_idx-1) - 1
      rite.str <- patternFeature$end + patternFeature$len * (c_idx-1) + 2
      rite.end <- patternFeature$end + patternFeature$len * c_idx + 1
      
      tmp_sign_delta <- mapply( function(ls,le,rs,re) mean(value[rs:re]) - mean(value[ls:le]),
                             left.str, left.end, rite.str, rite.end )
      
      tmp_sign_sub.delta <- mapply( function(ls,le,rs,re) mean(sub.value[rs:re]) - mean(sub.value[ls:le]), 
                                    left.str, left.end, rite.str, rite.end)
      if( use_Ap ) tmp_sign_sub.delta <- 0 * tmp_sign_sub.delta
      
      consistency_flag[ which(tmp_sign_delta > 0 & tmp_sign_sub.delta > 0) ] <- 0
      
      patternFeature$h2        <- patternFeature$h2        + tmp_sign_delta
      patternFeature$delta     <- patternFeature$delta     + tmp_sign_delta
      patternFeature$sub.delta <- patternFeature$sub.delta + tmp_sign_sub.delta
      
    }
    
    if( useConsistencyFlag ) patternFeature <- patternFeature[which(consistency_flag == 1),]
    patternFeature$h2        <-  patternFeature$h2        / (c_factor+1)
    patternFeature$delta     <-  patternFeature$delta     / (c_factor+1)
    patternFeature$sub.delta <-  patternFeature$sub.delta / (c_factor+1)  
  }
  names(patternFeature) <- mapvalues( names(patternFeature), warn_missing=F, 
                                      c('str','end','len','str.t','end.t'), 
                                      c('start.idx','end.idx','pattern.size',
                                        'start.timestamp','end.timestamp') )
  
  return(patternFeature)
}



cell.summary <- function( cell, xColName = 'delta', yColName = 'sub.delta' ){
  
  cell  <- unique(cell[ order(cell$start.timestamp), ])
  cell2 <- data.frame( ts = as.numeric(cell$start.timestamp), 
                       main = cell[,xColName], 
                       sub  = cell[,yColName] )
  
  result <- summarize( cell2
                       , sum = length(ts)
                       , min.t = min(diff(ts))
                       , med.t = median(diff(ts))
                       , min.med.rate = min.t/med.t
                       , min.med.rate2 = quantile(diff(ts), .05)/med.t
                       , med.d = median(main)
                       , med.sub.d = median(sub)
                       , lost.sig.num = sum( pmax( round( diff(ts) / med.t ) -1 ,0))
                       , lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num))
  
  if('leftCell.exist'   %in%  names(cell)) result$leftCell.exist   <- unique(cell$leftCell.exist)
  if('bottomCell.exist' %in%  names(cell)) result$bottomCell.exist <- unique(cell$bottomCell.exist)
  
  names(result) <- mapvalues( names(result), c('med.d','med.sub.d'), 
                              c(paste0(xColName,'.med'),paste0(yColName,'.med')), warn_missing=F)
  return(result)
}

cell.summary2 <- function( cell, vec, xColName = 'delta', yColName = 'sub.delta' ){
  
  if( !missing(cell) ){
    
    cell  <- cell[ cell$start.timestamp > min(vec), ]
    cell  <- cell[ order(cell$start.timestamp), ]
    cell2 <- data.frame( ts = as.numeric(cell$start.timestamp), 
                         main = cell[,xColName], 
                         sub  = cell[,yColName] )
    
    
    if( nrow(cell) == 0 ) return( cell.summary2( vec = vec) )
    
    slot_number      <- findInterval( cell2$ts, vec=vec )
    slot_number.rle  <- rle(slot_number)
    slot_Property    <- slot_number.rle$lengths
    slot_ones        <- length( which( slot_Property == 1 ) ) 
    slot_zeros       <- length(vec) - length(slot_number.rle$values)
    effective.rate   <- slot_ones / length(vec) 
    cell2$start.distance <- cell2$ts - as.numeric(vec[ slot_number ])
    
    metric.info <- slot_zeros / slot_ones
    if( slot_ones > 0 ){
      result <- summarize( cell2[ slot_number %in% slot_number.rle$values[ slot_Property == 1 ], ]
                           , Ontime.med    = median(start.distance)
                           , Ontime.min    = min(start.distance)
                           , Ontime.max    = max(start.distance)
                           , Ontime.sd     = sd(start.distance)
                           , powerDrop.med = median(sub)
                           , powerDrop.min = min(sub)
                           , powerDrop.max = max(sub)
                           , powerDrop.sd  = sd(sub)
                           , rpDrop.med    = median(main)
                           , rpDrop.min    = min(main)
                           , rpDrop.max    = max(main)
                           , rpDrop.sd     = sd(main) )
    }else{
      result <- data.frame( matrix(numeric(0), ncol = 12, nrow = 1)) 
    }
    result <- cbind( data.frame('slot_zeros'      = slot_zeros
                                , 'slot_ones'     = slot_ones
                                , 'metric.info'   = metric.info
                                , 'effective.rate' = effective.rate), result )
    
  }else{
    result <- data.frame( matrix(integer(0), ncol = 2,  nrow = 0),
                          matrix(numeric(0), ncol = 14, nrow = 0)) 
  }
  names(result) <- c( "slot_zeros", "slot_ones", "metric.info", 'effective.rate',
                      paste0("eff.TimeOn.",    c("med","min","max", "sd")), 
                      paste0("eff.PowerDrop.", c("med","min","max", "sd")),
                      paste0("eff.rpDrop.",    c("med","min","max", "sd")))
  if( !missing(cell) ){
    if('leftCell.exist'   %in%  names(cell)) result$leftCell.exist   <- unique(cell$leftCell.exist)
    if('bottomCell.exist' %in%  names(cell)) result$bottomCell.exist <- unique(cell$bottomCell.exist)
  }
 
  names(result) <- gsub( 'rpDrop',    xColName, names(result) )
  names(result) <- gsub( 'PowerDrop', yColName, names(result) )
  
  return( result )
}