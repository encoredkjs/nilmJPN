joinMeta <- function( data, meta1, meta2, check.identicalMeta = 1, ratio.threshold = .9, debug.mode = F ){
  
  meta1.signal <- meta2signal(meta1, data, postprocessing=F, show.fig=F)
  meta2.signal <- meta2signal(meta2, data, postprocessing=F, show.fig=F)
  if( !('p' %in% names(meta1.signal)) ) names(meta1.signal) <- mapvalues( names(meta1.signal), c('ap.box','rp.box'), c('p','q'))
  if( !('p' %in% names(meta2.signal)) ) names(meta2.signal) <- mapvalues( names(meta2.signal), c('ap.box','rp.box'), c('p','q'))
  
  if( debug.mode ){
    par(mfrow=c(3,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot( data$timestamp, data$active_power, type='l', ann=FALSE,   xaxt="n")
    plot( meta1.signal$timestamp, meta1.signal$p, type='l', ann=FALSE,   xaxt="n")
    plot( meta2.signal$timestamp, meta2.signal$p, type='l', ann=FALSE,   xaxt="n")
    title(paste('Compare two meta\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Meta2 / Meta2 / Original active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  
  signal.merge <- merge( meta1.signal, meta2.signal, by='timestamp', all = T )
  signal.merge$on.x <- signal.merge$p.x > 0
  signal.merge$on.y <- signal.merge$p.y > 0
  signal.merge$p.z  <- pmin( signal.merge$p.x, signal.merge$p.y )
  
  if( check.identicalMeta == 1 ){
    meta1.signalSummary <- series.to.box.lists(signal.merge$timestamp,signal.merge$p.x,0)[[1]]
    meta2.signalSummary <- series.to.box.lists(signal.merge$timestamp,signal.merge$p.y,0)[[1]]
    meta3.signalSummary <- series.to.box.lists(signal.merge$timestamp,signal.merge$p.z,0)[[1]]
    
    tmp1 <- apply( meta3.signalSummary, 1, 
                   function(row) which((meta1.signalSummary$str <= row['str']) & (row['end'] <= meta1.signalSummary$end)))
    tmp2 <- apply( meta3.signalSummary, 1, 
                   function(row) which((meta2.signalSummary$str <= row['str']) & (row['end'] <= meta2.signalSummary$end)))
    
    true.ratio <- 
      (rep( meta3.signalSummary$duration, 2 ) / meta1.signalSummary$duration[tmp1] > ratio.threshold) &&
      (rep( meta3.signalSummary$duration, 2 ) / meta1.signalSummary$duration[tmp2] > ratio.threshold)

    is.identical <- any( true.ratio )
  }else{
    print(ddply( signal.merge, .(on.x, on.y), summarize, len=length(timestamp) ))
    
    true.ratio   <- nrow(subset( signal.merge, on.x & on.y )) / min( length(which(signal.merge$on.x)),
                                                                     length(which(signal.merge$on.y)) )  
    is.identical <- any( true.ratio > ratio.threshold )
  }
  
  if( is.identical ){
    
    common.name <- setdiff( union(names(meta1), names(meta2)), c('meta-version','generation_info'))
    meta <- mapply(c, meta1[common.name], meta2[common.name], SIMPLIFY=FALSE)
    meta['meta-version'] <- 1
    
    str.t <- as.character(min(data$timestamp))
    end.t <- as.character(max(data$timestamp))
    
    meta[['generation_info']] <- list('data_used' = list( 'start' = str.t, 'end' = end.t, 'sampling' = 1),
                                    'computed' = as.character(Sys.time()))
    
    if( 'cycle' %in% names(meta) ){
      
      signal.merge$p <- pmax( signal.merge$p.x, signal.merge$p.y )
      signalSummary  <- series.to.box.lists(signal.merge$timestamp,signal.merge$p,0)[[1]]
      
      duration.on  <- signalSummary$duration
      duration.on  <- duration.on[!duration.on %in% boxplot.stats(duration.on)$out] # outlier 제거
      duration.off <- signalSummary$str[-1] - signalSummary$end[-nrow(signalSummary)]
      duration.off <- duration.off[ duration.off > 2]
      duration.off <- duration.off[!duration.off %in% boxplot.stats(duration.off)$out] # outlier 제거
      
      on.off <- c('working_time' = mean(duration.on), 
                  'working_time_sigma' = sd(duration.on),
                  'duty_cycle' = mean(duration.on) / (mean(duration.on)+mean(duration.off)))
      meta[['cycle']] <- on.off
    }
    return(meta)
  }else{
    return(list())
  }
}



meta2box_merge <- function( meta1, meta2, data, postprocessing = T, show.fig = F, filename = "meta2box_merge" ){
  
  meta1.signal <- meta2signal(meta1, data, postprocessing=F, show.fig=show.fig)
  meta2.signal <- meta2signal(meta2, data, postprocessing=F, show.fig=show.fig)
  if( !('p' %in% names(meta1.signal)) ) names(meta1.signal) <- mapvalues( names(meta1.signal), c('ap.box','rp.box'), c('p','q'))
  if( !('p' %in% names(meta2.signal)) ) names(meta2.signal) <- mapvalues( names(meta2.signal), c('ap.box','rp.box'), c('p','q'))
  
  signal.merge   <- merge( meta1.signal, meta2.signal, by='timestamp', all = T )
  signal.merge$p <- pmax( signal.merge$p.x, signal.merge$p.y )
  
  if( postprocessing ) box.shape <- post.processing(meta.json, pre.box.shape, show.fig)
  
}

