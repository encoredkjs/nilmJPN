speed2state <- function(speed){
  
  idx.nneg.speed <- which( speed >= 0 ) # index of non-negative speed
  idx.neg.speed  <- which( speed <  0 ) # index of     negative speed
  
  state <- character(length = length(speed))
  state[idx.nneg.speed] <- 'R'
  state[idx.neg.speed ] <- 'F'
  state
}

MIN <- function(x) quantile( x, .05, names = F)
MED <- function(x) quantile( x, .50, names = F)
MAX <- function(x) quantile( x, .95, names = F)

calc.periodicity  <- function(x) setNames( MED(x) / MIN(x), 'periodicityRate' )
calc.locality     <- function(x) setNames( MAX(x) / MED(x), 'localityRate'    )
calc.TotalPeriod  <- function(x) diff(x$str)

DetectEvent_1Hz <- function( data
                             , majorColName = 'reactive_power'
                             , minorColName = 'active_power'
                             , ... ){
  
  slope            <- get_speed_new( data$timestamp, data[,majorColName] ) 
  compressedStates <- stateCompression( speed2state( slope ) )
  
  if( head(compressedStates$val,1) != 'R' ) compressedStates <- tail( compressedStates, -1 )
  if( tail(compressedStates$val,1) != 'F' ) compressedStates <- head( compressedStates, -1 )
  
  edgeInfo <- data.frame( str  = subset( compressedStates, val == "R" )$str, 
                          peak = subset( compressedStates, val == "R" )$end + 1, 
                          end  = subset( compressedStates, val == "F" )$end + 1 )
  
  colFeature <- function(colname,str,end,peak){
    strValue   <- data[str,  colname]
    endValue   <- data[end,  colname]
    peakValue  <- data[peak, colname]
    #left.value  <- sapply( edgeInfo$str, function(s) median(data[unique(pmin(pmax(s-c(5:0),1))),colname],na.rm=T))
    #rite.value  <- sapply( edgeInfo$end, function(e) median(data[unique(pmin(pmin(e+c(5:0),nrow(data)))),colname],na.rm=T))
    
    result <- data.frame(
        'str2peak' = peakValue -  strValue
      , 'peak2end' =  endValue - peakValue
      , 'str2end'  =  endValue -  strValue
      #      , 'delta'    = rite.value - left.value
    )
    names(result) <- paste( names(result), colname, sep = '.')
    return(result)
  }
  
  feature  <- lapply( c(majorColName,minorColName), colFeature, str = edgeInfo$str, end = edgeInfo$end, peak = edgeInfo$peak )
  edgeInfo <- cbind( edgeInfo, dplyr::bind_cols(feature) )
  
  edgeInfo <- subset(edgeInfo, eval(parse(text = deparse(...))))
  if( nrow(edgeInfo) == 0 ) return( edgeInfo )
  
  tmp <- which( edgeInfo$end == shift(edgeInfo$str,-1))
  if( length(tmp) > 0 ){
    peak <-
      sapply( c(majorColName, minorColName), function(colname){
        colname <- paste0('str2end.', colname)
        abs( edgeInfo[ tmp, colname ] + edgeInfo[ tmp+1, colname ] ) < 3
      })
    if( (!is.null(dim(peak)) && any(removepeak <- apply( peak, 1, all ))) ||
        is.null(dim(peak)) && any(removepeak <- all( peak) )) edgeInfo <- edgeInfo[ -tmp[removepeak], ]
    #edgeInfo <- edgeInfo[ -(tmp+1), ]
  } 

  tmp <- which( edgeInfo$end == shift(edgeInfo$str,-1))
  if( length(tmp) > 0 ){
  
    contiguousEdge <- split( tmp, cumsum(c(1, diff(tmp) != 1)) )
    contiguousEdge <- lapply( contiguousEdge, function(x) c(x,max(x)+1) )
    
    newStr  <- edgeInfo$str[ sapply( contiguousEdge, function(i) min(i) ) ]
    newEnd  <- edgeInfo$end[ sapply( contiguousEdge, function(i) max(i) ) ]
    newPeak <- edgeInfo$peak[ sapply( contiguousEdge, function(is) is[which.max(data$reactive[edgeInfo$peak[is]])]) ]

    newFeature <- lapply( c(majorColName,minorColName), colFeature, str = newStr, end = newEnd, peak = newPeak )    
    edgeInfo <- rbind( edgeInfo[-c(tmp,tmp+1),], 
                       cbind( data.frame( str = newStr, peak = newPeak, end = newEnd ), dplyr::bind_cols(newFeature) ))
  } 
  
  return( edgeInfo %>% arrange(str) )
}


generate.CyclicBox.meta_extend <- function( data, max.iter = 3
                                            , debug.mode = F
                                            , split.algorithm = 1
                                            , level = 1
                                            , maxLevel = 4
                                            , bin.threshold = 5
                                            , mag.active = 20
                                            , mag.reactive = 3
                                            , known_meta
                                            , useAdaptiveExtensionP = FALSE ){
  
  # projection에 이용
  if( !('timeIdx' %in% names(data)) ) data$timeIdx <- 1:nrow(data)
  
  if( level %/% 3 == 0 ){
    if( split.algorithm == 1 ){
      
      n.hour <- 4
      data.split <- dlply( data, .(format(timestamp,'%F')))
      data.split <- data.split[ sapply( data.split, function(x) length(unique(hour(x$timestamp))) > n.hour ) ]
      hourly.sum <- lapply( data.split, function(x) 
        ddply(x, .(hour=format(timestamp,'%F %H')), summarize, usage = mean(active_power)))
      
      str.t <- sapply( hourly.sum, function(x) x$hour[which.min(stats::filter(x$usage, rep(1,4), sides=2))-1] )
      str.t <- as.POSIXct( str.t, format = '%F %H', tz = attr(data$timestamp,'tzone')) 
      end.t <- str.t + dhours(4) - dseconds(1)
      
      data.list  <- mapply( function(s,e) subset( data, s <= timestamp & timestamp <= e ), 
                            str.t, end.t, SIMPLIFY = F, USE.NAMES = T )
      data.list <- data.list[ sapply( data.list, function(x) length(unique(hour(x$timestamp))) == n.hour ) ]
      
    }else if( split.algorithm == 2 ){
      data.list <- split.data1( data=data,  debug.mode=debug.mode, apply.envelopeDetector = apply.envelopeDetector)
    }
  }else{
    data.list <- list( data )
    names(data.list) <- as.character(min(floor_date(data$timestamp,'day'))) 
  }
  
  if(missing(known_meta)){ known_meta <- list() }
  edge.lists <- list()
  
  if( level %% 2 == 1 ){
    DetectEvent   <- function(data) DetectEvent_1Hz( data, 'reactive_power', 'active_power', 
                                                     substitute( abs(str2end.reactive_power) >= mag.reactive & abs(str2end.active_power) >= mag.active, 
                                                                 list('mag.reactive' = mag.reactive, 'mag.active' = mag.active)))
  }else{
    DetectEvent   <- function(data) DetectEvent_1Hz( data, 'active_power', 'reactive_power', 
                                                     substitute( abs(str2end.reactive_power) <= mag.reactive & abs(str2end.active_power) >= mag.active, 
                                                                 list('mag.reactive' = mag.reactive, 'mag.active' = mag.active)))
  }
  totalEdge   <- DetectEvent( data )
  partialEdge <- ldply( data.list, DetectEvent, .id = 'day' )
  #sapply( data.list, function(x) c(min(x$timeIdx),max(x$timeIdx))  )
  para4clust <- c( 'str2end.reactive_power', 'str2end.active_power', 'str2peak.reactive_power' ) 
  para4meta  <- c( 'str2peak.reactive_power', 'str2end.reactive_power', 
                   'str2peak.active_power', 'str2end.active_power' )
  
  if( nrow(partialEdge) > 0 && !missing(known_meta) && length(known_meta) > 0 ){
    for( i in seq_along(known_meta) ){
      totalEdge    <- meta2edgeDetect( known_meta[[i]], totalEdge, para4meta )[[1]]
      partialEdge  <- rbind.fill( dlply( partialEdge, .(day), 
                                         function(localEdge) meta2edgeDetect( known_meta[[i]], localEdge, para4meta )[[1]] ))
    }
  }
  
  known_para  <- list()
  known_rEdge <- list()
#   while( iter <= max.iter ){
#     
#   }
   for( iter in 1:(max.iter+1) ){
    
    if( iter <= max.iter ){
      para <- search.cyclic.box.meta_extend( df = partialEdge
                                             , listOfColNames = para4clust
                                             , debug.mode = debug.mode
                                             , bin.threshold = bin.threshold
                                             , known_para  = known_para
                                             , known_rEdge = known_rEdge)
    }else{
      para <- list()
    }
    
    
    if( is.null(para) || length(para) == 0 ){
      
      if( level == maxLevel ) return( known_meta )
      
      meta.lists <- generate.CyclicBox.meta_extend( data
                                                    , max.iter = max.iter
                                                    , debug.mode = debug.mode
                                                    , split.algorithm = split.algorithm
                                                    , level = level + 1
                                                    , maxLevel = maxLevel
                                                    , bin.threshold = bin.threshold
                                                    , known_meta = known_meta
                                                    , useAdaptiveExtensionP = useAdaptiveExtensionP )
      return( meta.lists )
    }
    
    print( para[1:3] )
    meta <- from.para.to.JSON( data, para[[1]], para[[2]], para[[3]], para[[1]]['n.pt'] )
    known_para <- append( known_para, para[4] )
    
    nrowOfPartialEdge <- nrow(partialEdge)
    edgeDetected <- meta2edgeDetect( meta, totalEdge, para4meta, useAdaptiveExtensionP = useAdaptiveExtensionP )
    totalEdge    <- edgeDetected[[1]]
    partialEdge  <- partialEdge[ (partialEdge$str + sapply( data.list, function(x) min(x$timeIdx)-1 )[partialEdge$day])
                                 %in% totalEdge$str, ]
    detectedEG   <- c( nrowOfPartialEdge - nrow(partialEdge) ) / 2
    
    if( (!is.null(edgeDetected[[2]]) & !is.null(edgeDetected[[3]])) && 
        (nrow(edgeDetected[[2]])>=bin.threshold & nrow(edgeDetected[[3]])>=bin.threshold)  ){
      
      rEdge.info <- edgeDetected[[2]] %>% summarise_each_(funs(mean,sd,MIN,MAX,median), para4meta)
      fEdge.info <- edgeDetected[[3]] %>% summarise_each_(funs(mean,sd,MIN,MAX,median), para4meta)
      rEdge.info <- merge( rEdge.info, data.frame( n.pt = nrow(edgeDetected[[2]]),
                                                   locality = calc.locality(diff(edgeDetected[[2]]$str) )))
      fEdge.info <- merge( fEdge.info, data.frame( n.pt = nrow(edgeDetected[[2]])))

      replaceStr <- list('MIN'='min','MAX'='max','_mean'='_height','_sd'='_sigma','_median'='_med')
      for( i in seq_along(replaceStr) ){
        names(rEdge.info) <- gsub(names(replaceStr)[i], replaceStr[[i]], names(rEdge.info))
        names(fEdge.info) <- gsub(names(replaceStr)[i], replaceStr[[i]], names(fEdge.info))
      }
      rEdge.info <- setNames(as.numeric(rEdge.info), names(rEdge.info))
      fEdge.info <- setNames(as.numeric(fEdge.info), names(fEdge.info))
      
      on.duration <- edgeDetected[[3]]$end - edgeDetected[[2]]$str
      on.duration <- on.duration[ MIN(on.duration) <= on.duration & on.duration <= MAX(on.duration)]
      tot.duration <- diff(edgeDetected[[2]]$str)
      para[[3]]['working_time']        <- median(on.duration)
      para[[3]]['working_time_sigma']  <-     sd(on.duration)
      para[[3]]['working_time_min']    <-    min(on.duration)
      para[[3]]['working_time_max']    <-    max(on.duration)
      para[[3]]['tot_time']            <- median(tot.duration)
      
      meta <- from.para.to.JSON( data,  rEdge.info[grep('_',names(rEdge.info))]
                                 , fEdge.info[grep('_',names(fEdge.info))], para[[3]], rEdge.info['n.pt']
                                 , rEdge.info['locality'] )
      meta$level <- level
      
      if( round( para[[2]]['n.pt'] / para[[1]]['n.pt'] ) < 2 || para[[2]]['n.pt'] / detectedEG < 2 ){
        known_meta[[length(known_meta)+1]] <- meta
        edge.lists[[length(edge.lists)+1]] <- edgeDetected[2:3]
      }else{
        known_rEdge <- rbind( known_rEdge, para[[5]] )
      }
    }
  }
  return(known_meta)
}


search.cyclic.box.meta_extend <- function( df
                                           , listOfColNames
                                           , bin.threshold = 5
                                           , width.threshold = 50
                                           , debug.mode = F
                                           , max.iter = 20
                                           , periodicityRate.max = 2
                                           , effective.rate.min = .5
                                           , known_para  = NULL
                                           , known_rEdge = NULL){

  if( nrow(df) == 0 ) return(NULL)
  
  # clustering
  df <- edgeClustering( EG = df
                        , listOfColNames = listOfColNames
                        , debug.mode = debug.mode
                        , method = 'original'
                        , max.iter = max.iter, h = .01 )
  
  
  df.subset <- subset( df, cluster %in% as.numeric(names(table(df$cluster))[table(df$cluster) > bin.threshold]) )
  if( nrow(df.subset) == 0 || length(unique(df.subset$cluster)) < 2 ) return(NULL)
  
  if( !is.null(known_rEdge) && length(known_rEdge) > 0 ){
    df.subset.tmp <- df.subset %>% select(-cluster) %>% mutate( rowNumber = 1:nrow(df.subset) )
    df.subset.tmp <- anti_join(df.subset.tmp, select( known_rEdge, -cluster))
    df.subset     <- df.subset[ sort(df.subset.tmp$rowNumber), ]
  }
  
  if( debug.mode ){
    colName4plot <- paste0( 'str2end.', c('active_power', 'reactive_power'))
    plot( df[,colName4plot], col='grey70' )
    points( df.subset[,colName4plot], col = as.numeric(factor(df.subset$cluster)) + 1, pch = 15 )
  }
  
  cluster.info <- df.subset %>% group_by(cluster) %>% 
    summarise_each_(funs(mean,sd,min,max,median), listOfColNames) 
  cluster.info <- merge( cluster.info, ddply( df.subset, .(cluster), summarize, n.pt = length(str) ))
  
  time.diff <- dlply( df.subset, .(cluster), function(xx) 
    unlist(dlply( xx, .(day), function(xxx) diff(sort(xxx$str)) ), use.names = F))
  cluster.info <- merge( cluster.info, ldply( time.diff, calc.periodicity ), all = T )
  
  rising  <- subset( cluster.info, str2end.active_power_mean > 0 & periodicityRate < periodicityRate.max  )
  falling <- subset( cluster.info, str2end.active_power_mean < 0 )
  if( nrow(rising) == 0 || nrow(falling) == 0 ) return( NULL )
  
  if( !is.null(known_para) && length(known_para) > 0 ){
    rising.tmp <- rising %>% select(-cluster) %>% mutate( rowNumber = 1:nrow(rising) )
    rising.tmp <- anti_join( rising.tmp, select(ldply(known_para, .id=NULL ),-cluster) )
    rising     <- rising[ sort(rising.tmp$rowNumber), ]
  }
  
  efftMat <- outer( rising$cluster, falling$cluster, Vectorize( function(r,f){
    rEdge.sub <- subset( df.subset, cluster == r )
    fEdge.sub <- subset( df.subset, cluster == f )
    r.split   <- split( rEdge.sub, rEdge.sub$day, drop = F )
    f.split   <- split( fEdge.sub, fEdge.sub$day, drop = F )
    
    fillIn <- mapply( function(s,e) split(e,findInterval(e$str,v=s$str)), 
                      r.split, f.split, USE.NAMES = F, SIMPLIFY = F )
    effective.rate <- table(sapply( unlist(fillIn, recursive = F), nrow ))
    effective.rate <- as.numeric( effective.rate['1'] / nrow(rEdge.sub) )
    effective.rate
  }))
  
  distMat <- dist( rbind(  rising[,c('str2end.active_power_median','str2end.reactive_power_median')], 
                        -falling[,c('str2end.active_power_median','str2end.reactive_power_median')]))
  distMat <- as.matrix( distMat )[1:nrow(rising),-c(1:nrow(rising)), drop=FALSE]
  distMat[ is.na(efftMat) | effective.rate.min > efftMat | 1 < efftMat ] <- Inf

  rowNumber <- which.max(rising$n.pt)
#   if( !any(is.finite(distMat[rowNumber,])) ){
#     
#   }
   for( i in 1:nrow(rising) ){
    #rowNumber <- which.max(rising$n.pt)
    rowNumber <- order(rising$n.pt, decreasing = T)[i]
    if( any(is.finite(distMat[rowNumber,]))) break
  }
  if( i == nrow(rising) & !any(is.finite(distMat[rowNumber,])) ) return(NULL)
  colNumber <- which.min(distMat[rowNumber,])
  rising    <- rising[ rowNumber, ]
  falling   <- falling[ colNumber,] 
  
  if( abs((rising['str2end.active_power_mean'] + falling['str2end.active_power_mean']) / 
          rising['str2end.active_power_mean']) > 1 ){
    rising <- subset( cluster.info, str2end.active_power_mean > 0 )
    dist <- ddply( rising, .(cluster), function(r) 
      sqrt((r$str2end.active_power_mean + falling$str2end.reactive_power_mean)^2 + 
             (r$str2end.active_power_mean + falling$str2end.reactive_power_mean)^2 ) )
    rising <- rising[which.min(dist$V1),]
  }
  
  if( !any(is.finite(distMat[rowNumber,])) || (nrow(rising) * nrow(falling) == 0) || 
      (sign(rising$str2end.active_power_mean) == sign(falling$str2end.active_power_mean)) ) 
    return(NULL)
  
  rEdge <- subset( df.subset, cluster == rising$cluster  ) 
  fEdge <- subset( df.subset, cluster == falling$cluster ) 
  
  if( debug.mode ){
    plot( df[,colName4plot], col='grey70' )
    points( rEdge[,colName4plot], col='red')
    points( fEdge[,colName4plot], col='blue')
  }
  
  rEdge.split <- split( rEdge, rEdge$day, drop = F )
  fEdge.split <- split( fEdge, fEdge$day, drop = F )
  
  for( i in seq_along(rEdge.split) ){
    rEdge.split[[i]] <- subset( rEdge.split[[i]], str < max(fEdge.split[[i]]$str) )
    fEdge.split[[i]] <- subset( fEdge.split[[i]], str > min(rEdge.split[[i]]$str) )
  }
  
  rEdge.group <- mapply( function(r,f) split.data.frame( r, findInterval( r$str, vec = f$str )), 
                         rEdge.split, fEdge.split, SIMPLIFY = F )
  fEdge.group <- mapply( function(r,f) split.data.frame( f, findInterval( f$str, vec = r$str )), 
                         rEdge.split, fEdge.split, SIMPLIFY = F )
  rEdge.group <- unlist( rEdge.group, recursive = F )
  fEdge.group <- unlist( fEdge.group, recursive = F )
  
  matchingResults <- 
    mapply( function(r,f) as.numeric(outer( f$end, r$str, '-')), rEdge.group, fEdge.group, SIMPLIFY = F )
  slotOne    <- sapply( matchingResults, length ) == 1

  matching.r <- ldply( rEdge.group[ slotOne ], .id = NULL )
  matching.f <- ldply( fEdge.group[ slotOne ], .id = NULL )
  if( is.null(matching.r) || is.null(matching.f) || nrow(matching.r) == 0 || nrow(matching.f) == 0 ) return(NULL)

  onDuration  <- matching.f$str - matching.r$str
  totalPeriod <- setdiff(unlist(dlply( matching.r, .(day), calc.TotalPeriod ),use.names = F),0)
  
  measureVars <- names(matching.r)[ intersect( which(sapply(matching.r, is.numeric)), 
                                               grep('power',names(matching.r)) ) ]
  matching.r  <- matching.r %>% summarise_each_(funs(min, max), measureVars)
  matching.f  <- matching.f %>% summarise_each_(funs(min, max), measureVars)
  
  rEdge.summary  <- cbind(  rising[,setdiff( names(rising), names(matching.r) )], matching.r )
  fEdge.summary  <- cbind( falling[,setdiff(names(falling), names(matching.f) )], matching.f )
  rEdge.summary  <- rEdge.summary[, !(names(rEdge.summary) %in% c('cluster', 'value'))]
  fEdge.summary  <- fEdge.summary[, !(names(fEdge.summary) %in% c('cluster', 'value'))]
  rEdge.summary  <- setNames( as.numeric(rEdge.summary), names(rEdge.summary) )
  fEdge.summary  <- setNames( as.numeric(fEdge.summary), names(fEdge.summary) )
  
  replaceStr <- list('_mean'='_height','_sd'='_sigma','_median'='_med')
  for( i in seq_along(replaceStr) ){
    names(rEdge.summary) <- gsub(names(replaceStr)[i], replaceStr[[i]], names(rEdge.summary))
    names(fEdge.summary) <- gsub(names(replaceStr)[i], replaceStr[[i]], names(fEdge.summary))
  }
  
  duration.summary <- 
    setNames( c(median(onDuration), sd(onDuration), MIN(onDuration), MAX(onDuration)
                , median(onDuration)/ifelse(length(totalPeriod) == 0, 1, median(totalPeriod)))
              , c(paste0('working_time',c('','_sigma','_min','_max')),'duty_cycle') )

  return( list( 'rising' = rEdge.summary, 'falling' = fEdge.summary, 
                'duration' = duration.summary, 'knownPara.r' = rising,
                'rEdge' = rEdge ))
}

meta2edgeDetect <- function( meta
                             , edge
                             , listOfColNames
                             , extension.p = .2
                             , flexibility.p = .15
                             , debug.mode = FALSE
                             , ignoreRP = FALSE
                             , useHardMatching = TRUE
                             , useSoftMatching = FALSE
                             , useAdaptiveExtensionP = FALSE){
  
  if( ignoreRP ) listOfColNames <- listOfColNames[ -grep('reactive_power', listOfColNames) ]
  paraRange  <- sapply( listOfColNames, function(x) outer( x, c('_min','_max'), paste0), simplify = F )
  
  # rising edge condition
  rEdge.para <- lapply(  paraRange, function(x) meta$rising_edge[x] )
  if( useAdaptiveExtensionP &&
      ('str2peak.reactive_power' %in% names(rEdge.para)) && 
      max( max(abs(rEdge.para$str2peak.reactive_power)), 
           max(abs(rEdge.para$str2end.reactive_power)) ) < 3 ){
    
    extension.p.vec <- rep( extension.p, length(rEdge.para) )
    extension.p.vec[ grepl('reactive_power', names(rEdge.para)) ] <- 
      extension.p.vec[ grepl('reactive_power', names(rEdge.para)) ] * 5
    extension.p.vec[ !grepl('reactive_power', names(rEdge.para)) ] <- 
      extension.p.vec[ !grepl('reactive_power', names(rEdge.para)) ] / 2
    
    rEdge.para <- mapply( function(x, weight)  x * c( 1-weight*sign(x[1]), 1+weight*sign(x[2]) ), 
                          rEdge.para, extension.p.vec, SIMPLIFY = FALSE, USE.NAMES = TRUE )
  }else{
    rEdge.para <- lapply( rEdge.para, function(x) x * c( 1-extension.p*sign(x[1]), 1+extension.p*sign(x[2]) ))
  }
  
  # falling edge condition
  fEdge.para <- lapply(  paraRange, function(x) meta$falling_edge[x] )
  if( useAdaptiveExtensionP &&
      ('str2peak.reactive_power' %in% names(fEdge.para)) && 
      max( max(abs(fEdge.para$str2peak.reactive_power)), 
           max(abs(fEdge.para$str2end.reactive_power)) ) < 3 ){
    
    extension.p.vec <- rep( extension.p, length(fEdge.para) )
    extension.p.vec[ grepl('reactive_power', names(fEdge.para)) ] <- 
      extension.p.vec[ grepl('reactive_power', names(fEdge.para)) ] * 5
    extension.p.vec[ !grepl('reactive_power', names(fEdge.para)) ] <- 
      extension.p.vec[ !grepl('reactive_power', names(fEdge.para)) ] / 2
    
    fEdge.para <- mapply( function(x, weight)  x * c( 1-weight*sign(x[1]), 1+weight*sign(x[2]) ), 
                          fEdge.para, extension.p.vec, SIMPLIFY = FALSE, USE.NAMES = TRUE )
  }else{
    fEdge.para <- lapply( fEdge.para, function(x) x * c( 1-extension.p*sign(x[1]), 1+extension.p*sign(x[2]) ))
  }
  # duration condition
  duration.med <- meta$cycle['working_time']
  #duration.min <- meta$cycle['working_time_min'] * (1-flexibility.p)
  #duration.max <- meta$cycle['working_time_max'] * (1+flexibility.p)
  duration.min <- meta$cycle['working_time']
  duration.max <- meta$cycle['working_time']
  duration.min <- meta$cycle['working_time_min']
  duration.max <- meta$cycle['working_time'] * 2 - meta$cycle['working_time_min']
  duration.min <- duration.min * (1-flexibility.p)
  duration.max <- duration.max * (1+flexibility.p)
  
  listOfColNames.r <- listOfColNames
  listOfColNames.f <- listOfColNames
  listOfColNames.f <- listOfColNames.f[-grep('str2peak',listOfColNames.f)]
  
  condition1 <- sapply( listOfColNames.r, function(col){ 
    (rEdge.para[[col]][1] <= edge[,col]) & (edge[,col] <= rEdge.para[[col]][2])})
  condition2 <- sapply( listOfColNames.f, function(col){
    (fEdge.para[[col]][1] <= edge[,col]) & (edge[,col] <= fEdge.para[[col]][2])})
  if( is.matrix(condition1) ){
    condition1 <- apply( condition1, 1, all )
  }else{
    condition1 <- all(condition1)
  }
  if( is.matrix(condition2) ){
    condition2 <- apply( condition2, 1, all )
  }else{
    condition2 <- all( condition2 )
  }
  
  rEdge <- edge[ condition1, ]
  fEdge <- edge[ condition2, ]
  if( nrow(rEdge) == 0 || nrow(fEdge) == 0 )
    return( list(edge, data.frame(), data.frame(), NULL, NULL) )
  
  rEdge0 <- subset( rEdge, str > max(fEdge$end) )
  rEdge  <- subset( rEdge, str < max(fEdge$end) )
  if( nrow(rEdge) == 0 )
    return( list(edge, data.frame(), data.frame(), rEdge0, NULL) )

  fEdge0 <- subset( fEdge, str < min(rEdge$str) )
  fEdge  <- subset( fEdge, str > min(rEdge$str) )
  if( nrow(fEdge) == 0 )
    return( list(edge, data.frame(), data.frame(), rEdge0, fEdge0) )
  
  rEdge.group <- split.data.frame( rEdge, findInterval( rEdge$str, vec = fEdge$str ))
  fEdge.group <- split.data.frame( fEdge, findInterval( fEdge$str, vec = rEdge$str ))
  
  hardMatchingResults <-  mapply( hardMatching, rEdge.group, fEdge.group, SIMPLIFY = F, 
                                  MoreArgs = list( minPulseWidth = duration.min, 
                                                   maxPulseWidth = duration.max, 
                                                   pulseWidth = duration.med ) )
  
  rEdge.match <- rbind.fill(lapply( hardMatchingResults, function(x) x$rEdge ))
  fEdge.match <- rbind.fill(lapply( hardMatchingResults, function(x) x$fEdge ))
  
  rEdge.group <- rbind.fill(lapply( hardMatchingResults, function(x) x$rEdgeRest ))
  fEdge.group <- rbind.fill(lapply( hardMatchingResults, function(x) x$fEdgeRest ))
  
  edge <- subset( edge, !(str %in% rEdge.match$str | str %in% fEdge.match$str) )
  if( !useSoftMatching || is.null(rEdge.match) || 
      nrow(rEdge.match) < 2 || length(undetected.event(rEdge.match$str,total.period = meta$cycle['tot_time'])) == 0 ) 
    return( list(edge, rEdge.match, fEdge.match, rEdge0, fEdge0) )
  
  #return( list(edge, rEdge.match, fEdge.match) )
  
  duration.min <- duration.min * (1-flexibility.p)
  duration.max <- duration.max * (1+flexibility.p)
  
  totEvent.r <- total.event(rEdge.match$str, meta$cycle['tot_time'])
  totEvent.f <- totEvent.r
  totEvent.f$eventTime[  totEvent.f$is.detected ] <- fEdge.match$str
  totEvent.f$eventTime[ !totEvent.f$is.detected ] <- totEvent.f$eventTime[ !totEvent.f$is.detected ] + duration.med
  
  emptySlot <- cbind( which( totEvent.r$is.detected & shift(!totEvent.r$is.detected,-1)), 
                      which(!totEvent.r$is.detected & shift( totEvent.r$is.detected,-1)))
  emptySlot <- cbind( emptySlot, emptySlot[,2] - emptySlot[,1] )
  emptySlot <- alply( emptySlot, 1 )
  
  emptySlotSub.r <- lapply( emptySlot, function(row) 
    seq( totEvent.r$event[row[1]], totEvent.r$event[row[2]+1], length.out = row[3] + 1) )
  emptySlotSub.f <- lapply( emptySlot, function(row) 
    seq( totEvent.f$event[row[1]], totEvent.f$event[row[2]+1], length.out = row[3] + 1) )
  emptySlotSub.r <- rbind.fill.matrix(lapply( emptySlotSub.r, function(x) cbind(head(x,-1), tail(x,-1)) ))
  emptySlotSub.f <- rbind.fill.matrix(lapply( emptySlotSub.f, function(x) cbind(head(x,-1), tail(x,-1)) ))
  
  rEdge.group <- alply( emptySlotSub.r, 1, function(sub) subset( rEdge.group, sub[1] <= str & str < sub[2] ))
  fEdge.group <- alply( emptySlotSub.f, 1, function(sub) subset( fEdge.group, sub[1] <= str & str < sub[2] ))
  
  rEdge.group.extend <- alply( emptySlotSub.r, 1, function(sub) subset( edge, str2end.active_power > 0 & sub[1] <= str & str < sub[2] ))
  fEdge.group.extend <- alply( emptySlotSub.f, 1, function(sub) subset( edge, str2end.active_power < 0 & sub[1] <= str & str < sub[2] ))
  
  softMatchingResults <-  mapply( hardMatching, 
                                  rEdge.group[sapply(rEdge.group,nrow)==1], 
                                  fEdge.group.extend[sapply(rEdge.group,nrow)==1], SIMPLIFY = F, 
                                  MoreArgs = list( minPulseWidth = duration.min, 
                                                   maxPulseWidth = duration.max, 
                                                   pulseWidth = duration.med ) )
  rEdge.match <- rbind( rEdge.match, rbind.fill(lapply( softMatchingResults, function(x) x$rEdge )) )
  fEdge.match <- rbind( fEdge.match, rbind.fill(lapply( softMatchingResults, function(x) x$fEdge )) )
  
  softMatchingResults <-  mapply( hardMatching, 
                                  rEdge.group.extend[sapply(fEdge.group,nrow)==1], 
                                  fEdge.group[sapply(fEdge.group,nrow)==1], SIMPLIFY = F, 
                                  MoreArgs = list( minPulseWidth = duration.min, 
                                                   maxPulseWidth = duration.max, 
                                                   pulseWidth = duration.med ) )
  rEdge.match <- rbind( rEdge.match, rbind.fill(lapply( softMatchingResults, function(x) x$rEdge )) )
  fEdge.match <- rbind( fEdge.match, rbind.fill(lapply( softMatchingResults, function(x) x$fEdge )) )
  
  edge <- subset( edge, !(str %in% rEdge.match$str | str %in% fEdge.match$str) )
  rEdge.match <- rEdge.match[order(rEdge.match$str),]; rEdge.match <- rEdge.match[order(rEdge.match$str),]
  fEdge.match <- fEdge.match[order(fEdge.match$str),]; fEdge.match <- fEdge.match[order(fEdge.match$str),]
  
  return( list(edge, rEdge.match, fEdge.match, rEdge0, fEdge0) )
}


multipleMeta2cyclicBoxExtend <- function( meta.lists
                                          , data
                                          , debug.mode = FALSE
                                          , saveAsFile = FALSE
                                          , listOfColNames = c( 'str2peak.reactive_power', 'str2peak.active_power', 
                                                                'str2end.reactive_power', 'str2end.active_power')
                                          , postProcessing = FALSE
                                          , mag.reactive = 3
                                          , mag.active = 20
                                          , extension.p = .2
                                          , flexibility.p = .15
                                          , ignoreRP = FALSE
                                          , useSoftMatching = TRUE
                                          , keepSingleEdge = FALSE
                                          , filename = 'multiplemeta2cyclic'
                                          , locality.max = 10
                                          , useAdaptiveExtensionP = FALSE ){
  
  if( debug.mode && saveAsFile ) png( paste0(filename,"-%d.png"), width=3000)
  
  level <- unique(sapply( meta.lists, function(x) x$level ))
  
  totalEdge <- data.frame()
  for( l in unique(level %% 2) ){
    
    if( l == 1 ){
      if( ignoreRP ){
        totalEdgeTmp <- DetectEvent_1Hz( data, 'reactive_power', 'active_power', 
                                         substitute( abs(str2end.active_power) >= mag.active, list('mag.active' = mag.active)))
      }else{
        totalEdgeTmp <- DetectEvent_1Hz( data, 'reactive_power', 'active_power', 
                                         substitute( abs(str2end.reactive_power) >= mag.reactive & abs(str2end.active_power) >= mag.active, 
                                                     list('mag.reactive' = mag.reactive, 'mag.active' = mag.active)))
      }
    }else{
      if( ignoreRP ){
        totalEdgeTmp <- DetectEvent_1Hz( data, 'active_power', 'reactive_power', 
                                         substitute( abs(str2end.active_power) >= mag.active, list('mag.active' = mag.active)))
      }else{
        totalEdgeTmp <- DetectEvent_1Hz( data, 'active_power', 'reactive_power', 
                                         substitute( abs(str2end.reactive_power) <= mag.reactive & abs(str2end.active_power) >= mag.active, 
                                                     list('mag.reactive' = mag.reactive, 'mag.active' = mag.active)))
      }
    }
    totalEdge <- rbind( totalEdge, totalEdgeTmp )
  }
  
  totalEdge <- totalEdge[ order(totalEdge$str), ]
  
  detectedEdge <- lapply( meta.lists, function(meta) 
    meta2edgeDetect( meta
                     , totalEdge
                     , listOfColNames = listOfColNames 
                     , extension.p = extension.p
                     , flexibility.p = flexibility.p
                     , ignoreRP = ignoreRP
                     , useSoftMatching = useSoftMatching
                     , useAdaptiveExtensionP = useAdaptiveExtensionP ) )
  #diag(checkCommonEdge) / apply( checkCommonEdge, 1, sum )
  checkCommonEdge <- outer( detectedEdge, detectedEdge, Vectorize(function(x,y){length(intersect(x[[2]]$str,y[[2]]$str))} ) )
  diag(checkCommonEdge) <- 1
  multiplyingEdge <- unique(unlist(apply( checkCommonEdge, 1, function(x) list(which(x!=0)) ), recursive = F))
  nonmultiplyingEdge <- unlist(multiplyingEdge[ sapply( multiplyingEdge, length ) == 1 ])
  multiplyingEdge <- multiplyingEdge[ sapply( multiplyingEdge, length ) != 1 ]
  
  if( length(multiplyingEdge) > 0){
    for( i in seq(multiplyingEdge) ){
      for( j in seq(multiplyingEdge) ){
        if( i < j ){
          if( length(intersect(multiplyingEdge[[i]], multiplyingEdge[[j]])) > 0 ){
            multiplyingEdge[[i]] <- union(multiplyingEdge[[i]],multiplyingEdge[[j]])
            multiplyingEdge[[j]] <- NA
          }
        } 
      }
    }
    multiplyingEdge <- multiplyingEdge[sapply( multiplyingEdge, function(x) unique(!is.na(x)) )]
  }
  
  result       <- list()
  result.order <- NULL
  if( length(nonmultiplyingEdge) > 0 ){
    result       <- lapply( detectedEdge[ nonmultiplyingEdge ], function(x) x[2:3] )
    result.order <- c( result.order, nonmultiplyingEdge )
  } 
  if( (length(nonmultiplyingEdge) > 0) && keepSingleEdge ){
    tmp <- lapply( detectedEdge[ nonmultiplyingEdge ], function(x) x[4:5] )
    exist.singleEdge <- sapply( tmp, function(x) !is.null(x[[1]]) )
    if( any(exist.singleEdge) ){
      for( i in which(exist.singleEdge) ){
        if( tmp[[i]][[1]]$str > nrow(data) - meta.lists[[ nonmultiplyingEdge ]]$cycle['working_time'] * 2 ){
          result[[i]][[1]] <- rbind( result[[i]][[1]], tmp[[i]][[1]] )
          tmp[[i]][[1]]$end <- nrow(data)
          result[[i]][[2]] <- rbind( result[[i]][[2]], tmp[[i]][[1]] )
          result[[i]][[1]] <- result[[i]][[1]][ order(result[[i]][[1]]$str),]
          result[[i]][[2]] <- result[[i]][[2]][ order(result[[i]][[2]]$str),]
          result.order <- c( result.order, nonmultiplyingEdge )
        }
      }
    }
  }
  
  if( length(multiplyingEdge) != 0 ){
    for( i in seq_along(multiplyingEdge) ){
      meta.local <- meta.lists[multiplyingEdge[[i]]]
      shortPeriodic <- order(sapply( meta.local, function(x) x$cycle['working_time'] ))
      
      totalEdgeTmp <- totalEdge
      for( j in shortPeriodic ){
        detectedEdge <-  
          meta2edgeDetect( meta.local[[j]]
                           , totalEdgeTmp
                           , listOfColNames = listOfColNames 
                           , extension.p = extension.p
                           , flexibility.p = flexibility.p
                           , ignoreRP = ignoreRP
                           , useSoftMatching = useSoftMatching )
        result[[length(result)+1]] <- detectedEdge[2:3]
        result.order <- c( result.order, multiplyingEdge[[i]][j] )
        totalEdgeTmp <- detectedEdge[[1]]
      }
    }
  }
  
  signalResult <- vector("list", length(meta.lists))
  for( i in seq_along(result) ){
    result_i <- result[[i]]
    #    result_i[[2]]$end[result_i[[2]]$end == result_i[[1]]$end] <- nrow(data)
    loc <- unlist(mapply( function(r,f) seq(r,f), result_i[[1]]$str, result_i[[2]]$end))
    
    signalTmp <- data.frame( timestamp = data$timestamp, p = 0, q = 0)
    signalTmp$p[loc] <- (meta.lists[[result.order[i]]]$rising_edge['str2end.active_power_med'] - 
                           meta.lists[[result.order[i]]]$falling_edge['str2end.active_power_med'])/2
    signalTmp$q[loc] <- (meta.lists[[result.order[i]]]$rising_edge['str2end.reactive_power_med'] - 
                           meta.lists[[result.order[i]]]$falling_edge['str2end.reactive_power_med'])/2
    #signalResult[[ length(signalResult) + 1 ]] <- signalTmp
    signalResult[[ result.order[i] ]] <- signalTmp
  }
  
  if( debug.mode ){
    par(mfrow=c(length(signalResult)+1,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot( data$timestamp, data$active_power, type='l')
    for( i in seq_along(signalResult) ) plot( signalResult[[i]]$timestamp, signalResult[[i]]$p, type='l')
    
    plot( data$timestamp, data$reactive_power, type='l')
    for( i in seq_along(signalResult) ) plot( signalResult[[i]]$timestamp, signalResult[[i]]$q, type='l')
    par(mfrow=c(1,1))
  }
  
  if( postProcessing ){
    for( i in seq_along(result) ){
      if( 'locality' %in% names(meta.lists[[i]]) && meta.lists[[i]]['locality'] < locality.max ){
        signalResult[[i]] <- post.processing( meta.lists[[i]], signalResult[[i]], debug.mode = F, calcWTbound = F )
        names(signalResult[[i]]) <- c('timestamp','p.residual','q.residual','p','q')
      }else{
        message( 'Skip post-processing :', i, '-th meta is not cyclic => locally cyclic' )
      }
    }
    if( debug.mode ){
      par(mfrow=c(length(signalResult)+1,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
      plot( data$timestamp, data$active_power, type='l')
      for( i in seq_along(signalResult) ){
        if('locality' %in% names(meta.lists[[i]]) && meta.lists[[i]]['locality'] < locality.max ){
          plot( signalResult[[i]]$timestamp, signalResult[[i]]$p, type='l')
        }else{
          plot( signalResult[[i]]$timestamp, signalResult[[i]]$p, type='l', col='grey70')
        }
      } 
    }
  }
  
  
  if( debug.mode && saveAsFile ) dev.off()
  
  return(signalResult)
}

# minPulseWidth와 maxPulseWidth 사이에 rEdge와 fEdge가 matching이 되는 경우에 대해 pairing
hardMatching <- function( rEdgeCandidates, fEdgeCandidates, minPulseWidth, maxPulseWidth, pulseWidth ){
  
  rEdgeCandidates0 <- subset( rEdgeCandidates, str > max(fEdgeCandidates$end) )
  rEdgeCandidates  <- subset( rEdgeCandidates, str < max(fEdgeCandidates$end) )
  if( nrow(rEdgeCandidates) == 0 ) 
    return( list(rEdgeRest=rEdgeCandidates, fEdgeRest=fEdgeCandidates, rEdge=NULL, fEdge=NULL) )
  
  fEdgeCandidates0 <- subset( fEdgeCandidates, str < min(rEdgeCandidates$str) )
  fEdgeCandidates  <- subset( fEdgeCandidates, str > min(rEdgeCandidates$str) )
  if( nrow(fEdgeCandidates) == 0 ) 
    return( list(rEdgeRest=rEdgeCandidates, fEdgeRest=fEdgeCandidates, rEdge=NULL, fEdge=NULL) )
  
  pulseWidthTable <- outer( fEdgeCandidates$end, rEdgeCandidates$str, '-' )
  pulseWidthTable[ minPulseWidth > pulseWidthTable ] <- NA
  pulseWidthTable[ maxPulseWidth < pulseWidthTable ] <- NA
  if( all(is.na(pulseWidthTable)) ) 
    return( list(rEdgeRest=rEdgeCandidates, fEdgeRest=fEdgeCandidates, rEdge=NULL, fEdge=NULL) )
  
  bestFittingFairs <- abs( pulseWidthTable - pulseWidth )
  bestFittingFairs <- which( bestFittingFairs == min(bestFittingFairs, na.rm=T), arr.ind = T )
  
  rEdge <- rEdgeCandidates[bestFittingFairs[,2],]
  fEdge <- fEdgeCandidates[bestFittingFairs[,1],]
  
  rEdgeRest <- rEdgeCandidates[-bestFittingFairs[,2],]
  fEdgeRest <- fEdgeCandidates[-bestFittingFairs[,1],]
  
  return( list(rEdgeRest=rEdgeRest, fEdgeRest=fEdgeRest, rEdge=rEdge, fEdge=fEdge) )
}

edgeClustering <- function( EG, listOfColNames, debug.mode = F, method = 'ms', max.iter = 20, ... ){
  
  if( all(! names(EG) %in% listOfColNames) ) return( rep(1,nrow(EG)))
  EG.subset <- EG[ , names(EG) %in% listOfColNames ]
  
  if( method == 'ms' ){
    EG$cluster <- LPCM::ms( EG.subset, ... )$cluster.label
  }
  if( method == 'original' ){
    
    findMaxGap  <- function(power) mean(sort(power)[which.max(diff(sort(power)))+c(0,1)])
    findCluster <- function(df,colName){
      if(nrow(df)<2) return(rep(1,nrow(df)))
      return(findInterval(df[,colName], vec = findMaxGap(df[,colName])))
    }  
    
    clusterIdx <- data.frame(matrix(1,nrow(EG),length(listOfColNames)))
    names(clusterIdx) <- gsub('X','cluster',names(clusterIdx))
    EG <- cbind( EG, clusterIdx ) 
    EG.list <- list( EG )
    
    for( iter in 1:max.iter ){
      
      check.periodicity <-  lapply( EG.list, function(edge){
        unlist(dlply(edge, .(day), function(x) diff(x$str)), recursive=F, use.names=F)
      })
      check.periodicity[sapply( check.periodicity, length ) == 0] <- 0
      need.more.division <- sapply( check.periodicity, calc.periodicity ) > 2
      need.more.division[ is.na(need.more.division) ] <- F
      need.more.division[ sapply( EG.list, function(x){
        nrow(x) > 1 && sign(max(x$str2end.active_power)) != sign(min(x$str2end.active_power))
      }) ] <- T
      
      convergedDf <- EG.list[!need.more.division]
      EG.list <- lapply( EG.list[need.more.division], function(x){
        for( i in seq_along(listOfColNames) ){
          x.split <- dlply(x, names(clusterIdx))
          cluster.new <- lapply( x.split, function(xx) findCluster(xx,listOfColNames[[i]]))
          cluster.new <- unlist( cluster.new, use.names = F)
          x <- rbind.fill(x.split)
          x[,paste0('cluster',i)] <- cluster.new
        }
        return( dlply( x, names(clusterIdx) )  )
      })
      EG.list <- unlist( EG.list, recursive = F )
      EG.list <- append( EG.list, convergedDf )
      
    }
    names(EG.list) <- seq_along(EG.list)
    EG <- ldply( EG.list, .id = 'cluster', stringsAsFactors = F)
    EG <- EG[, ! names(EG) %in% c(names(clusterIdx))]
    if(debug.mode){
      if( ncol(EG.subset) == 2 ){
        points( EG[,listOfColNames], col = as.numeric(EG$cluster) + 1 )
      } 
      if( ncol(EG.subset)  > 2 ){
        pairs( EG[,listOfColNames], col = as.numeric(EG$cluster) + 1 )
      }
    } 
  }
  return(EG)
}

