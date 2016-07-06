meta2PatternScan_1Hz <- function (data, meta.json, extension.p, flexibility.p, consistency.n, postprocessing = T, debug.mode, filename = "meta2pattern" ) {
  
  if(meta.json[['shape_type']] != 'pattern_scan') stop("Meta is not correct : meta2PatternScan_1Hz")
  
  if( debug.mode ) png( paste0(filename,"-%d.png"), width=3000)
  
  result <- meta2PatternScan_summit_1Hz_new(data, meta.json, extension.p, flexibility.p, consistency.n, postprocessing=postprocessing, debug.mode = debug.mode)
  
  if( debug.mode ) dev.off()
  
  return(result)  
}

meta2PatternScan_summit_1Hz_extend <- function( data
                                                , meta
                                                , extension.p = .2
                                                , flexibility.p = .9
                                                , consistency.n = 1
                                                , postprocessing = F
                                                , debug.mode = F
                                                , minsupportRatio = .5
                                                , fillMissingData = F
                                                , filename = 'meta2patternExtend') {
  
  if( debug.mode ) png( paste0(filename,"-%d.png"), width=3000)
  
  if( meta$supportRatio < minsupportRatio ) postprocessing <- F
  
  summit_flag <- as.character(meta$rising_edge['summit_flag'])

  # pattern matching
  InfoExt.p <- extension.p
  m.flex.p  <- flexibility.p
  r.edge    <- meta$rising_edge
  f.edge    <- meta$falling_edge
  
  s.rp_min <- r.edge['rp_min'] *(1-InfoExt.p*sign(r.edge['rp_min']))
  s.rp_max <- r.edge['rp_max'] *(1+InfoExt.p*sign(r.edge['rp_max']))
  s.rp_med <- r.edge['rp_med']
  
  s.ap_min <- r.edge['ap_min'] *(1-InfoExt.p*sign(r.edge['ap_min']))
  s.ap_max <- r.edge['ap_max'] *(1+InfoExt.p*sign(r.edge['ap_max']))
  s.ap_med <- r.edge['ap_med']
    
  s.t_min <- r.edge['min.t'] *(1-InfoExt.p)
  s.t_med <- r.edge['med.t']
  s.t_min.med.rate <- r.edge['min.med.rate']
  
  e.rp_min <- f.edge['rp_min'] *(1-InfoExt.p*sign(f.edge['rp_min']))
  e.rp_max <- f.edge['rp_max'] *(1+InfoExt.p*sign(f.edge['rp_max']))
  e.rp_med <- f.edge['rp_med']
  
  e.ap_min <- f.edge['ap_min'] *(1-InfoExt.p*sign(f.edge['ap_min']))
  e.ap_max <- f.edge['ap_max'] *(1+InfoExt.p*sign(f.edge['ap_max']))
  e.ap_med <- f.edge['ap_med']
  
  e.eff_t_min <- f.edge['EffTimeOn.min'] *(1-InfoExt.p)
  e.eff_t_med <- f.edge['EffTimeOn.med']
  e.eff_t_max <- e.eff_t_med * (1+m.flex.p)
  
  e.eff_ap_med <- f.edge['EffAP_Drop.med']
  e.eff_p_sd   <- f.edge['EffAP_Drop.sd' ]  
  e.eff_Rp_med <- f.edge['EffRP_Drop.med']
  e.eff_Rp_sd  <- f.edge['EffRP_Drop.sd' ]
    
  main_type = switch( summit_flag, '1' = 'reactive', '-1' = '-reactive', '2' = 'active'   )
  sub_type  = switch( summit_flag, '1' = 'active',   '-1' = 'active',    '2' = 'reactive' )
  
  answer.log <- data.frame( timestamp = data$timestamp, p= 0, q= 0)
  
  # search start points
  s.pattern <- DetectPattern_1Hz_new( data = data
                                      , position = "start"
                                      , main_type = main_type
                                      , sub_type = sub_type
                                      , debug.mode = F
                                      , use_Ap = (main_type == 'active'))
  
  if( summit_flag == 2 ){
    s.pattern <- subset( s.pattern, (delta >= s.ap_min) & (delta <= s.ap_max) & 
                           (sub.delta >= s.rp_min) & (sub.delta <= s.rp_max) )
  }else{
    s.pattern <- subset( s.pattern, (delta >= s.rp_min) & (delta <= s.rp_max) & 
                           (sub.delta >= s.ap_min) & (sub.delta <= s.ap_max) )
  }

  
  if( nrow(s.pattern) == 0 ){
    print("Rising edge detect : failed")
    return( answer.log )
  }else if( debug.mode ){
    plot( data$active_power, type='l')
    abline( v = s.pattern$start.idx, col='red' )
  }
    
  s.pattern <- s.pattern[order(s.pattern$start.idx),] 
  
  # examine s.pattern
  while( min(diff(as.numeric(s.pattern$start.timestamp))) < s.t_min ) {
          
    t1.idx <- which.min(diff(as.numeric(s.pattern$start.timestamp)))
    t2.idx <- t1.idx + 1
    t.idx  <- c(t1.idx, t2.idx)
    
    likelihood <- sapply( t.idx, 
                          function(t.idx) sqrt(((s.ap_med - s.pattern$delta[t.idx] )/s.ap_med )^2 +
                                               ((s.rp_med - s.pattern$sub.delta[t1.idx] )/s.rp_med )^2 ))
    
    s.pattern <- s.pattern[ - t.idx[ which.max( likelihood ) ],]
    if (likelihood[1] == likelihood[2]) print("A start point has been removed randomly")
    
  }
  
  # search end points
  e.pattern <- DetectPattern_1Hz_new( data = data
                                      , position = "end"
                                      , main_type = main_type
                                      , sub_type  = sub_type
                                      , c_factor = consistency.n
                                      , debug.mode = F
                                      , use_Ap = (main_type == 'active'))
  
  if( summit_flag != 2 ){
    e.pattern <- subset( e.pattern, (delta >= e.rp_min) & (delta <= e.rp_max) & 
                           (sub.delta >= e.ap_min) & (sub.delta <= e.ap_max) )
  }else{
    e.pattern <- subset( e.pattern, (delta >= e.ap_min) & (delta <= e.ap_max) & 
                           (sub.delta >= e.rp_min) & (sub.delta <= e.rp_max) )
  }

  
  if( nrow(e.pattern) == 0 ){
    print("Falling edge detect : failed")
    return( answer.log )
  }else if( debug.mode ){
    abline( v = e.pattern$start.idx, col='blue' )
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
  }
  
  e.pattern <- e.pattern[ order(e.pattern$start.idx),]
  
  
  # edge selection
  t.duration <- list()
  chosen.ends <- list()
  
  e.pattern          <- subset( e.pattern, start.timestamp >= min(s.pattern$end.timestamp) )
  slot_number        <- findInterval( e.pattern$start.timestamp, vec=s.pattern$end.timestamp )
  e.pattern$duration <- as.numeric(e.pattern$start.timestamp) - as.numeric(s.pattern$end.timestamp[ slot_number ])
  e.pattern          <- subset( e.pattern, (e.eff_t_min < duration) & (duration < e.eff_t_max))
  
  slot_number        <- findInterval( e.pattern$start.timestamp, vec=s.pattern$end.timestamp )
  slot_number.rle    <- rle(slot_number)
  slot_Property      <- slot_number.rle$lengths
  slot_ones          <- length( which( slot_Property == 1 ) ) 
  
  if(main_type == 'active'){
    likelihoodFn <- function(x) abs((e.eff_ap_med  - as.numeric(x$delta) / e.eff_ap_med ))
  }else{
    likelihoodFn <- function(x) sqrt( ((e.eff_Rp_med - as.numeric(x$delta) )/e.eff_Rp_med )^2 +
                                        ((e.eff_ap_med  - as.numeric(x$sub.delta)/e.eff_ap_med )^2 ))
  }
  if( any(slot_Property != 0) ){
    matched.ends <- lapply( 1:length(s.pattern$end.timestamp), function(x){
      e.pattern.sub <- e.pattern[ slot_number == x, ]
      if( nrow(e.pattern.sub) == 0 ) return(NULL)
      calc_likelihood <- likelihoodFn( e.pattern.sub )
      e.pattern.sub[ which.min(calc_likelihood), ]
    })
    tmp_cnt <- length(which(sapply( matched.ends, is.null )))
  }
  resultant.info <- data.frame( start.timestamp = s.pattern$start.timestamp[which(!sapply( matched.ends, is.null ))], 
                                end.timestamp = ldply(matched.ends)$end.timestamp )
  
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
        
    if( any(sapply( matched.ends, is.null )) ){
      resultant.info <- rbind( resultant.info, 
                               data.frame( start.timestamp = s.pattern$start.timestamp[which(sapply( matched.ends, is.null ))], 
                                           end.timestamp = s.pattern$start.timestamp[which(sapply( matched.ends, is.null ))] + e.eff_t_med))
    }
    resultant.info <- resultant.info[order(resultant.info$start.timestamp),]
    
    undetected <- undetected.event( as.numeric(resultant.info$start.timestamp), s.t_med )
    undetected <- as.POSIXct(undetected, origin = "1970-01-01", tz = "Asia/Seoul")
      
    ind <- as.numeric(min(resultant.info$end.timestamp)) - e.eff_t_med:0 - s.t_med
    while( any(ind >= as.numeric(min(data$timestamp))) ){
      str.t <- as.POSIXct(min(ind), origin = "1970-01-01", tz = "Asia/Seoul")
      end.t <- as.POSIXct(max(ind), origin = "1970-01-01", tz = "Asia/Seoul")
      resultant.info <- rbind( resultant.info, 
                               data.frame( start.timestamp = str.t, end.timestamp = end.t ) )
      ind <- ind - s.t_med
    }
    
    ind <- as.numeric(max(resultant.info$start.timestamp)) + 0:e.eff_t_med + s.t_med
    while( any(ind <= as.numeric(max(data$timestamp))) ){
      str.t <- as.POSIXct(min(ind), origin = "1970-01-01", tz = "Asia/Seoul")
      end.t <- as.POSIXct(max(ind), origin = "1970-01-01", tz = "Asia/Seoul")
      resultant.info <- rbind( resultant.info, 
                               data.frame( start.timestamp = str.t, end.timestamp = end.t ) )
      ind <- ind + s.t_med
    }
    
    add.box <- undetected
    if( length(add.box) > 0 ){
      str.t <- undetected
      end.t <- undetected + e.eff_t_med - 1
      resultant.info <- rbind( resultant.info, 
                               data.frame( start.timestamp = str.t, end.timestamp = end.t ) )
    }
  }
  
  if (nrow(resultant.info) != 0) {
    
    # log extraction
    p.pad <- mean( c(s.ap_med, abs(e.eff_ap_med)), na.rm=T )
    q.pad <- mean( c(s.rp_med, abs(e.eff_Rp_med)), na.rm=T )
    
    if(fillMissingData){
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
      mtext("Signal / Active power", 2, 3, outer=TRUE, las=0)
      par(mfrow=c(1,1))
      
      dev.off()
    }
    return(answer.log)
  } else {
    stop("Error: no end point!")
  }
}

generate.PatternScan.meta_1Hz_extend <- function ( data
                                                   , genResolution.n  = 20
                                                   , genEffSize.n     = 15
                                                   , staPeriodicity.p = .1
                                                   , endEffSlot.p     = .1
                                                   , endConsistency.n = 1
                                                   , min_sigmag.active.n = 20
                                                   , min_sigmag.reactive.n = 3
                                                   , min_lossrate.p    = .8
                                                   , clustering.method = 3
                                                   , debug.mode = F
                                                   , use_posRpAp = T
                                                   , use_negRpAp = T
                                                   , use_Ap = T) {
                                                       
  
  # parameters
  s.max_resolution.n <- genResolution.n
  e.max_resolution.n <- genResolution.n
  s.min_gsize.n      <- genEffSize.n
  g.eff_sample.n     <- genEffSize.n 
  e.eff_rate.p       <- endEffSlot.p      # select small -> big
  e.eff_rate_ext.p   <- e.eff_rate.p / 2
  
  # search starting points
  s.pattern.original <- DetectPattern_1Hz_extended( data, "start"
                                                    , min_sigmag.active.n
                                                    , min_sigmag.reactive.n
                                                    , debug.mode = debug.mode
                                                    , use_posRpAp = use_posRpAp
                                                    , use_negRpAp = use_negRpAp
                                                    , use_Ap = use_Ap)
  
  e.pattern.original <- DetectPattern_1Hz_extended( data, "end"
                                                    , min_sigmag.active.n
                                                    , min_sigmag.reactive.n
                                                    , endConsistency.n 
                                                    , debug.mode = debug.mode
                                                    , use_posRpAp = use_posRpAp
                                                    , use_negRpAp = use_negRpAp
                                                    , use_Ap = use_Ap)
  
  json.result <- list()
  
  for( level in 1:2 ){
    
    if( level == 1 &  any(s.pattern.original$type != 'active') ){
      s.pattern <- subset( s.pattern.original, type != 'active')
      e.pattern <- subset( e.pattern.original, type != 'active')
    }else{
      s.pattern <- subset( s.pattern.original, type == 'active')
      e.pattern <- subset( e.pattern.original, type == 'active')
    }
    
    if( nrow(s.pattern) != 0 & nrow(e.pattern) != 0 ){
      
      chosen.e.idx <- c()
      chosen.e.data.list <- list()
      chosen.e.data.summary <- list()
      
      s.info <- edgeCluster( s.pattern
                             , xColName = switch( level, '1' = 'rp.delta', '2' = 'ap.delta' )
                             , yColName = switch( level, '1' = 'ap.delta', '2' = 'ap.delta' )
                             , resolution = s.max_resolution.n
                             , d_factor   = staPeriodicity.p
                             , gp_size    = s.min_gsize.n
                             , clusterMethod = clustering.method
                             , debug.mode = debug.mode
                             , xBuffer_mag = switch( level, '1' = min_sigmag.reactive.n, '2' = min_sigmag.active.n )
                             , yBuffer_mag = switch( level, '1' = min_sigmag.active.n,   '2' = min_sigmag.active.n ))
      
      if (length(s.info) == 0 ||
            nrow( s.list <- subset( getLast(s.info), lost.sig.rate / 100 <= min_lossrate.p)) == 0 ){
        
        print("Start point detection failed")
        
      }else{
        
        s.info <- s.info[s.list$candidate.idx]
        
        # detect corresponding end points for each start group
        print("Detected groups for starting points")
        print(s.list)
        
        for(s_idx in 1:nrow(s.list)){
          
          tmp.s.info <- s.info[[s_idx]]
          
          iter   <- 0 
          tmp_e.eff_rate.p <- e.eff_rate.p
          e.info <- list()
          while( length(e.info) == 0 || 
                   nrow( e.list <- subset(getLast(e.info), metric.info < 1)) == 0 ){
            
            iter <- iter + 1
            
            if( iter > 2 ){
              detection.failed <- T
              print("End point detection failed") 
              break
            } 
            
            if( iter > 1 ){
              print("Reducing 'effctive slot rate' to find end points again")
              tmp_e.eff_rate.p <- tmp_e.eff_rate.p - e.eff_rate_ext.p
            } 
            
            e.info <- edgeCluster( data = subset( e.pattern, type == unique(tmp.s.info$type) )
                                   , xColName = switch( level, '1' = 'rp.delta', '2' = 'ap.delta' )
                                   , yColName = switch( level, '1' = 'ap.delta', '2' = 'ap.delta' )
                                   , resolution = e.max_resolution.n
                                   , d_factor   = tmp_e.eff_rate.p
                                   , vec        = tmp.s.info$end.timestamp
                                   , clusterMethod = clustering.method
                                   , debug.mode = debug.mode
                                   , xBuffer_mag = switch( level, '1' = min_sigmag.reactive.n, '2' = min_sigmag.active.n )
                                   , yBuffer_mag = switch( level, '1' = min_sigmag.active.n,   '2' = min_sigmag.active.n ))
          }
          
          
          if ( iter <= 2 ){
            cat("End point candidate list for appliance ", s_idx, "\n")
            print(e.list)
            
            chooseInterval.1 <- findInterval( nrow(tmp.s.info), c(0,50,100) )
            if( chooseInterval.1 == 1 ){ # 0 <= nrow(s.info[[s_idx]] < 50
              chooseInterval.2 <- findInterval( e.list$metric.info, c(0, .9))
            }else if( chooseInterval.1 == 2 ){ #  50 <= nrow(s.info[[s_idx]] < 100
              chooseInterval.2 <- findInterval( e.list$metric.info, c(0, .6, .9))        
            }else if( chooseInterval.1 == 3 ){ # 100 <= nrow(s.info[[s_idx]]
              chooseInterval.2 <- findInterval( e.list$metric.info, c(0, .2, .6, .9))
            }
            e.list$metricInterval <- chooseInterval.2
            
            ex.e.list <- subset( e.list, metricInterval == min(chooseInterval.2) )
            
            
            rising_ap  <- s.list$ap.delta.med[s_idx]
            rising_rp  <- s.list$rp.delta.med[s_idx]
            falling_ap <- ex.e.list$eff.ap.delta.med
            falling_rp <- ex.e.list$eff.rp.delta.med
            
            tmp_compute <- 
              which.min( sqrt( ((falling_ap + rising_ap)/rising_ap )^2 +
                               ((falling_ap + rising_ap)/rising_ap )^2 ))
            
            chosen.e.idx <- c(chosen.e.idx, s_idx)
            chosen.e.data.list[[length(chosen.e.data.list)+1]]       <- e.info[[ex.e.list$candidate.idx[tmp_compute]]]
            chosen.e.data.summary[[length(chosen.e.data.summary)+1]] <- ex.e.list[tmp_compute, ]
          }       
        }      
      }
      
      if (length(chosen.e.idx) != 0) {
        
        str.t <- as.character(min(data$timestamp))
        end.t <- as.character(max(data$timestamp))
        
        # extract meta data
        for (m_idx in 1:length(chosen.e.idx)) {
          
          s.data.log <- s.info[[chosen.e.idx[m_idx]]]
          e.data.log <- chosen.e.data.list[[m_idx]]
          s.data.summary <- s.list[chosen.e.idx[m_idx],]
          e.data.summary <- chosen.e.data.summary[[m_idx]]
          
          supportRatio <- as.numeric( diff( c( min(s.data.log$start.timestamp),
                                               max(e.data.log$end.timestamp))), units='secs')
          supportRatio <- supportRatio / as.numeric(diff(range(data$timestamp)),units='secs')
          
          # start data info
          s.data.log.melt <- melt( s.data.log, id.vars='cluster', measure.vars = c('ap.delta','rp.delta'))
          s.data.log.cast <- cast( s.data.log.melt, cluster ~ variable, c(min,max,median))
          s.data.log.cast <- setNames( as.numeric(s.data.log.cast), names(s.data.log.cast) )
          names(s.data.log.cast) <- gsub( '.delta', '',    names(s.data.log.cast) )
          names(s.data.log.cast) <- gsub( 'median', 'med', names(s.data.log.cast) )
          s.data.log.cast <- setNames( as.numeric(s.data.log.cast), names(s.data.log.cast) )
          s.data.summary.numeric <- s.data.summary[,c('sum','min.t','med.t','min.med.rate'
                                                              ,'lost.sig.num','lost.sig.rate')]
          s.data.summary.numeric <- setNames( as.numeric(s.data.summary.numeric), names(s.data.summary.numeric) )
            
          tmp_start.json <- c( summit_flag = switch( paste0(unique(s.data.log$type), collapse='//'),
                                                     'reactive/active' = 1, 
                                                     '-reactive/active' = -1,
                                                     'active' = 2 ) 
                               , s.data.log.cast
                               , s.data.summary.numeric )
          
          # end data info
          e.data.log.melt <- melt( e.data.log, id.vars='cluster', measure.vars = c('ap.delta','rp.delta'))
          e.data.log.cast <- cast( e.data.log.melt, cluster ~ variable, c(min,max,median))
          e.data.log.cast <- setNames( as.numeric(e.data.log.cast), names(e.data.log.cast) )
          names(e.data.log.cast) <- gsub( '.delta', '',    names(e.data.log.cast) )
          names(e.data.log.cast) <- gsub( 'median', 'med', names(e.data.log.cast) )
          e.data.log.cast <- setNames( as.numeric(e.data.log.cast), names(e.data.log.cast) )
          
          tmp_end.json <- c( e.data.log.cast, 
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
                             EffRP_Drop.med = ifelse( is.null(e.data.summary$eff.rp.delta.med), NA, e.data.summary$eff.rp.delta.med ), 
                             EffRP_Drop.min = ifelse( is.null(e.data.summary$eff.rp.delta.min), NA, e.data.summary$eff.rp.delta.min ),
                             EffRP_Drop.max = ifelse( is.null(e.data.summary$eff.rp.delta.max), NA, e.data.summary$eff.rp.delta.max ), 
                             EffRP_Drop.sd  = ifelse( is.null(e.data.summary$eff.rp.delta.sd),  NA, e.data.summary$eff.rp.delta.sd  )) 

                             #colwise( function(x) c(min(x),median(x),max(x),sd(x)), )
          parameters <- c('genResolution.n'=genResolution.n
                          , 'genEffSize.n'=genEffSize.n
                          , 'staPeriodicity.p'=staPeriodicity.p
                          , 'endEffSlot.p'=endEffSlot.p
                          , 'endConsistency.n'=endConsistency.n
                          , 'min_sigmag.active.n'=min_sigmag.active.n
                          , 'min_sigmag.reactive.n'=min_sigmag.reactive.n
                          , 'min_lossrate.p'=min_lossrate.p
                          , 'clustering.method'=clustering.method)
          #                      , 'main_type'=main_type 
          #                      , 'sub_type'=sub_type )
          
          json.result[[length(json.result)+1]] <- 
            list(`meta-version` = 1, 
                 shape_type = "pattern_scan", 
                 rising_edge = tmp_start.json, 
                 falling_edge = tmp_end.json, 
                 supportRatio = supportRatio, 
                 generation_info = list(data_used = list(start = str.t, end = end.t, sampling = 1), 
                                        computed = as.character(Sys.time())),
                 parameters= parameters)
        }
      } 
    }
  }
  return(json.result)
}


DetectPattern_1Hz_extended <- function( data, position = c("start", "end")
                                        , min_mag_ap = 20
                                        , min_mag_rp = 20
                                        , c_factor = 0
                                        , debug.mode = F
                                        , use_posRpAp = T
                                        , use_negRpAp = T
                                        , use_Ap = T
                                        ){
  
  data.original   <- data
  result.normal   <- list()
  result.reverse  <- list()
  result.active   <- list()
  # normal direction
  if( use_posRpAp ){
    result.normal <- DetectPattern_1Hz_new( data, position, 'reactive', 'active', c_factor )
    if( position == "start" ){
      result.normal <- subset( result.normal, ( delta >= min_mag_rp) & ( sub.delta >= min_mag_ap) )
    }else if( position == "end" ){
      result.normal <- subset( result.normal, (-delta >= min_mag_rp) & (-sub.delta >= min_mag_ap) )
    }
    
    if( nrow(result.normal) > 0 ){
      result.normal$type <- 'reactive/active'
      names(result.normal) <- mapvalues( names(result.normal), warn_missing=F, 
                                         c('delta','h1','h2','org_h2','sub.delta','org_sub.delta'), 
                                         c('rp.delta','rp.h1','rp.h2','rp.org_h2','ap.delta','ap.org_sub'))
      
      loc <- unlist(apply( result.normal, 1, function(x) x['start.idx']:x['end.idx']), use.names=F)
      data[ loc, c('active_power','reactive_power') ]   <- NA
    }
  }
  
  # reverse direction
  if( use_negRpAp ){
    result.reverse <- DetectPattern_1Hz_new( data, position, '-reactive', 'active', c_factor )
    if( position == "start" ){
      result.reverse <- subset( result.reverse, (-delta >= min_mag_rp) & ( sub.delta >= min_mag_ap) )
    }else if( position == "end" ){
      result.reverse <- subset( result.reverse, ( delta >= min_mag_rp) & (-sub.delta >= min_mag_ap) )
    }
    
    if( nrow(result.reverse) > 0 ){
      result.reverse$type <- '-reactive/active'
      names(result.reverse) <- mapvalues( names(result.reverse), warn_missing=F, 
                                          c('delta','h1','h2','org_h2','sub.delta','org_sub.delta'), 
                                          c('rp.delta','rp.h1','rp.h2','rp.org_h2','ap.delta','ap.org_sub'))
      
      loc <- unlist(apply( result.reverse, 1, function(x) x['start.idx']:x['end.idx']), use.names=F)
      data[ loc, c('active_power','reactive_power') ]   <- NA
    }
  }
    
  # 
  if( use_Ap ){
    result.active <- DetectPattern_1Hz_new( data, position, 'active', 'reactive', c_factor, use_Ap=use_Ap )
    if( position == "start" ){
      result.active <- subset( result.active, ( delta >= min_mag_ap) & ( abs(sub.delta) < min_mag_rp) )
    }else if( position == "end" ){
      result.active <- subset( result.active, (-delta >= min_mag_ap) & ( abs(sub.delta) < min_mag_rp) )
    }
    
    if( nrow(result.active) > 0 ){
      result.active$type <- 'active'
      names(result.active) <- mapvalues( names(result.active), warn_missing=F, 
                                         c('delta','h1','h2','org_h2','sub.delta','org_sub.delta'), 
                                         c('ap.delta','ap.h1','ap.h2','ap.org_h2','rp.delta','rp.org_sub'))
      
      loc <- unlist(apply( result.active, 1, function(x) x['start.idx']:x['end.idx']), use.names=F)
      data[ loc, c('active_power','reactive_power') ]   <- NA
    }
  }
  
  result <- list( result.normal, result.reverse, result.active )
  result <- rbind.fill(result[ sapply(result, is.data.frame) ])
  result <- result[ order(result$start.idx), ]
    
  if(debug.mode){
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data.original$timestamp, y=data.original$active_power, ann=FALSE,   xaxt="n", type='l')
    abline( v=result$start.timestamp, col=factor(result$type))
    plot(x=data.original$timestamp, y=data.original$reactive_power, ann=FALSE, xaxt="n", type='l')
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)),outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  return( result )  
}

DetectPattern_1Hz_new <- function ( data, position = c("start", "end")
                                    , main_type = c("reactive", "active","-reactive", "-active")
                                    , sub_type  = c("active", "reactive","-active", "-reactive")
                                    , c_factor = 0
                                    , debug.mode = F
                                    , use_Ap = F
                                    , useConsistencyFlag = TRUE ){
  library(plyr)

  mainType <- match.arg(main_type)
  mainLog  <- switch(mainType,
                     'reactive'  = {data.frame(timestamp = data$timestamp, value =  data$reactive_power)},
                     'active'    = {data.frame(timestamp = data$timestamp, value =  data$active_power)},
                     '-reactive' = {data.frame(timestamp = data$timestamp, value = -data$reactive_power)},
                     '-active'   = {data.frame(timestamp = data$timestamp, value = -data$active_power)})
  
  subType <- match.arg(sub_type)
  subLog  <- switch(subType,
                    'reactive'  = {data.frame(timestamp = data$timestamp, value =  data$reactive_power)},
                    'active'    = {data.frame(timestamp = data$timestamp, value =  data$active_power)},
                    '-reactive' = {data.frame(timestamp = data$timestamp, value = -data$reactive_power)},
                    '-active'   = {data.frame(timestamp = data$timestamp, value = -data$active_power)})
  
      
  # for debug
  if ( nrow(mainLog) != nrow(subLog) ){
    cat("Log mismatch between active and reactive power!", nrow(mainLog), " and ", nrow(subLog), "\n")
    stop("Error")
  }
  
  slope <- get_speed_new(mainLog$timestamp, mainLog$value) 
  united.state <- stateCompression( scan_state_new(slope)  )
  
  edgeType <- match.arg(position)
  if ( edgeType == "start"){
    edge.pattern <- list( c(1,2,3,4), c(1,3,4), c(1,2,3,1), c(1,3,1) )
  }else if (edgeType == "end"){
    edge.pattern <- list( c(3,4), c(3,1), c(3,2) )
  }
  resultant <- ldply( edge.pattern, 
                      function( pattern ){ patternIdx <- compressedPatternMatching( united.state, pattern )
                                           get_patternFeature( mainLog, subLog, slope, patternIdx, c_factor, use_Ap, useConsistencyFlag )} )

  if ( edgeType == "start"){
    edge.pattern <- list( c(1,2), c(1,3) )
    resultant2 <- ldply( edge.pattern, 
                         function( pattern ){ patternIdx <- compressedPatternMatching( united.state, pattern )
                         get_patternFeature( mainLog, subLog, slope, patternIdx, c_factor, use_Ap, useConsistencyFlag )} )
    resultant2 <- subset( resultant2, !(start.idx %in% resultant$start.idx))
    if( nrow(resultant2) > 0 ) resultant <- rbind.fill( resultant, resultant2 )
  }  
  if( is.null(resultant) ) stop(paste('There is no proper', switch( edgeType , "start"='rising', "end"='falling')))
                                        
  if( substr(mainType,1,1) == '-' ){
    colNames <- names(resultant)[ names(resultant) %in% c('h1','h2','delta','org_h2') ]
    resultant[,colNames] <- resultant[,colNames] * (-1) 
  }
  if( substr(subType,1,1) == '-' ){
    colNames <- names(resultant)[ names(resultant) %in% c('sub.delta','org_sub.delta') ]
    resultant[,colNames] <- resultant[,colNames] * (-1) 
  }
    
  if( mainType == subType ) resultant <- resultant[, -grep("sub", names(resultant))]
    
  if(debug.mode){
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
    abline( v=resultant$start.timestamp, col=switch(edgeType, start = {"red"}, end = {"blue"}))
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
    if( edgeType == 'start' ){
      title(paste('Step1 : rising edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
    }else{
      title(paste('Step1 : falling edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
    }
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  
  return(resultant[ order(resultant$start.timestamp),] )
}