# doNILM <- function(meta, log, DRAW = FALSE) {
#   
#   app.meta <- fromJSON(toJSON(meta, pretty = TRUE))
#   
#   result <- lapply(app.meta, function(meta) {
#     shape_type = meta[["shape_type"]]
#     if (shape_type == "cyclic_box") {
#       temp <- meta2cyclic_box(meta.json = meta, data = log, 
#                               postprocessing = T, show.fig = DRAW)
#     }else if (shape_type == "high_power") {
#       temp <- meta2heavy_load(meta.json = meta, data = log, show.fig = DRAW)
#     }else {
#       temp <- data.frame(timestamp = log$timestamp, ap.box = 0, rp.box = 0)
#     }
#     return(temp)
#   })
#   
#   for (i in 1:length(result)) {
#     if (length(names(result[[i]])) == 5 ) {
#       names(result[[i]]) <- c('timestamp', 'p.reverse', 'q.reverse', 'p', 'q')
#       result[[i]] <- result[[i]][c('timestamp', 'p', 'q')]
#     }
#   }
#   return(result)
# }

data2each.signal <- function( filename  ){
  
  x <- read.csv( file = filename, header = F )
  names(x)=c('timestamp','voltage','current','active_power','reactive_power','apparent_power','power_factor')
  x$timestamp = as.POSIXct(x$timestamp,tz="Asia/Seoul",origin="1970-01-01")
  
  y = x[!duplicated(x),]
  
  if( difftime( y$timestamp[ max(which(y$active_power > 10)) ],
                y$timestamp[ min(which(y$active_power > 10)) ], units='hours') < 4 )
    stop('Not enough data : cyclic')
  
  if( length(which(y$active_power < quantile( y$active_power, .05) + 10 )) / length(y$active_power) > .5 )
    stop('Less than the base amount : cyclic')
  
  if( length(which(y$active_power < 0 )) / length(y$active_power) > .25 )
    stop('Negative active power : cyclic')
  
  if( length(which(y$active_power < 10 )) / length(y$active_power) > .8 )
    stop('No appliance : cyclic')
  
  #============================================================================================
  #
  #                                   Meta info 찾기
  #
  #============================================================================================
  
  #find cyclic boxes
  meta.list.cyclic = list()
  meta.list.cyclic <- tryCatch({
    generate.CyclicBox.meta(y,max.iter=3,debug.mode=T,split.algorithm=1)
  }, error = function(err) {list()})
  
  #find heavy loads
  meta.list.high = list()
  meta.list.high <- tryCatch({
    generate.HeavyLoad.meta(y,max.iter=3,debug.mode=T)
  }, error = function(err) {list()})
  
  
  #find pattern matching results
  meta.list.pattern = list()
  meta.list.pattern <- tryCatch({
    generate.PatternScan.meta_1Hz(data = y,
                                  genResolution.n = 250, genEffSize.n = 15, 
                                  staPeriodicity.p = 0.2, endEffSlot.p = 0.1)
  }, error = function(err) {list()})
  
  #merge list
  meta.list = list()
  meta.list = append(meta.list,meta.list.cyclic)
  meta.list = append(meta.list,meta.list.high)
  meta.list = append(meta.list,meta.list.pattern)
  
  if(length(meta.list)==0) return('There is no appliance')
  
  names(meta.list)=seq(0,length(meta.list)-1)
  meta.list.json <- toJSON( meta.list,pretty=TRUE )
  
  #============================================================================================
  #
  #                                   Box shape 찾기
  #
  #============================================================================================
  
  app.meta = fromJSON(meta.list.json)
  result <- lapply( app.meta, function(meta){ 
    shape_type = meta[['shape_type']] 
    if( shape_type == 'cyclic_box' ){ 
      temp <- meta2cyclic_box(meta.json=meta, data=y, postprocessing=T, show.fig=T )
    }else if( shape_type == 'high_power' ){
      temp <- meta2heavy_load(meta.json=meta, data=y,show.fig=T)
    }else if( shape_type == 'pattern_scan' ){
      temp <- meta2PatternScan_1Hz(data = y, meta = meta, extension.p =0.0, flexibility.p = 0.9)
    }else{ temp <- data.frame(timestamp = y$timestamp, ap.box = 0, rp.box = 0)}})
  
  return(result)
}


meta2signal <- function( meta, data, postprocessing=T, show.fig=T, extension.p =0.2, flexibility.p = .9, 
                         relApDiff.max = .15, relRpDiff.max = .2, measurementErrorBound = 10 ){
  library(plyr)
  
  sampling <- meta$generation_info$data_used$sampling
  data.original <- data
  ts <- data.original$timestamp
  
  if( sampling != 1 ){
    resolutionMin <- as.integer(sampling / 60)
    data.original$timestamp <- floor_date(data.original$timestamp,'hour') + 
      dminutes((minute(data.original$timestamp) %/% resolutionMin) * resolutionMin)
    data   <- smoothing( data, resolutionMin )
  } 
  
  shape_type = meta[['shape_type']]
  if( length(shape_type) == 1 ){
    if( shape_type == 'cyclic_box' ){ 
      result <- meta2cyclic_box(meta.json=meta, data=data, postprocessing=postprocessing, show.fig=show.fig )
    }else if( shape_type == 'high_power' ){
      result <- meta2heavy_load(meta.json=meta, data=data,show.fig=show.fig)
    }else if( shape_type == 'pattern_scan' ){
      result <- meta2PatternScan_1Hz(data = data, meta = meta, extension.p =extension.p, 
                                     flexibility.p = flexibility.p, consistency.n=1, 
                                     postprocessing=postprocessing, debug.mode=show.fig)
    }else if( shape_type == 'pattern_scan_heavy' ){
      result <- meta2PatternScanHeavy_1Hz(data = data, meta = meta, extension.p = extension.p, 
                                          relApDiff.max = relApDiff.max, relRpDiff.max = relRpDiff.max, 
                                          measurementErrorBound = measurementErrorBound, 
                                          debug.mode = show.fig)
    }else{ result <- data.frame(timestamp = data$timestamp, ap.box = 0, rp.box = 0)}
  }else if( length(shape_type) == 2 ){
    results <- lapply( shape_type, function(type){
      meta.local <- meta
      meta.local[['shape_type']] <- type
      signal <- meta2signal( meta.local, data, postprocessing=F, show.fig=F, extension.p =extension.p, flexibility.p = flexibility.p )
      return(signal[,c('timestamp','p','q')])
    })
    
    if( !identical( results[[1]]$timestamp, results[[2]]$timestamp ) ){
      print('timestamp does not match')
      return(list())
    }
    
    result <- data.frame( timestamp = results[[1]]$timestamp, 
                          p = apply( sapply( results, function(x) x$p ), 1, max ), 
                          q = apply( sapply( results, function(x) x$q ), 1, max )) 
    after.postProcessing <- post.processing( meta, result, debug.mode=show.fig )
    return( after.postProcessing )
  }else{
    print("Did not implement yet")
    return(list())
  }
  
  if( sampling != 1 ){
    result <- merge( data.original, result, all.x = T )
    result$timestamp <- ts
  }
  
  return(result)  
}

data2meta <- function( data ){
  
  library(RJSONIO)
  library(reshape)
  
  if( difftime( data$timestamp[ max(which(data$active_power > 10)) ],
                data$timestamp[ min(which(data$active_power > 10)) ], units='hours') < 4 )
    print('Not enough data : cyclic')
  
  if( length(which(data$active_power < quantile( data$active_power, .05) + 10 )) / length(data$active_power) > .5 )
    print('Less than the base amount : cyclic')
  
  if( length(which(data$active_power < 0 )) / length(data$active_power) > .25 )
    print('Negative active power : cyclic')
  
  if( length(which(data$active_power < 10 )) / length(data$active_power) > .8 )
    print('No appliance : cyclic')
  
  #============================================================================================
  #
  #                                   Meta info 찾기
  #
  #============================================================================================
  
  #find cyclic boxes
  meta.list.cyclic <- list()
  meta.list.cyclic <- tryCatch(
    expr  = generate.CyclicBox.meta(data,max.iter=5,debug.mode=F,split.algorithm=1),
    error = function(err) {
      print(err$message)
      list()
    })
  
  #find heavy loads
  meta.list.high <- list()
  meta.list.high <- tryCatch(
    expr  = generate.HeavyLoad.meta( data, max.iter = 5, debug.mode = F ),
    error = function(err) {
      print(err$message)
      list()
    })
  
  #find heavy loads
  meta.list.patternHigh <- list()
  meta.list.patternHigh <- tryCatch(
    expr  = generate.PatternScanHeavy.meta_1Hz_extend( data
                                                       , relApDiff.max = .15
                                                       , relRpDiff.max = .20
                                                       , min_sigmag.active.n = 800
                                                       , min_sigmag.reactive.n = 3
                                                       , timeSpan = ehours(1)
                                                       , measurementErrorBound = 10
                                                       , debug.mode = F )[[3]],
    error = function(err) {
      print(err$message)
      list()
    })
  
  #find pattern matching results
  meta.list.pattern1 <- list()
  meta.list.pattern1 <- tryCatch(
    expr  = generate.PatternScan.meta_1Hz( data = data, 
                                           genResolution.n = 1000, genEffSize.n = 15, 
                                           staPeriodicity.p = 0.1, endEffSlot.p = 0.1, 
                                           endConsistency.n = 1,
                                           clustering.method = 1, debug.mode = F),
    error = function(err) {
      print(err$message)
      list()
    })
  
  meta.list.pattern2 <- list()
  meta.list.pattern2 <- tryCatch(
    expr  = generate.PatternScan.meta_1Hz( data = data, 
                                           genResolution.n = 1000, genEffSize.n = 15, 
                                           staPeriodicity.p = 0.1, endEffSlot.p = 0.15, 
                                           endConsistency.n = 1,
                                           clustering.method = 1, debug.mode = F),
    error = function(err) {
      print(err$message)
      list()
    })
  
  meta.list.pattern3 <- list()
  meta.list.pattern3 <- tryCatch(
    expr  = generate.PatternScan.meta_1Hz( data = data, 
                                           genResolution.n = 20, genEffSize.n = 15, 
                                           staPeriodicity.p = 0.1, endEffSlot.p = 0.1, 
                                           endConsistency.n = 1,
                                           clustering.method = 3, debug.mode = F),
    error = function(err) {
      print(err$message)
      list()
    })
  
  meta.list.pattern4 <- list()
  meta.list.pattern4 <- tryCatch(
    expr  = generate.PatternScan.meta_1Hz( data = data, 
                                           genResolution.n = 20, genEffSize.n = 15, 
                                           staPeriodicity.p = 0.1, endEffSlot.p = 0.15, 
                                           endConsistency.n = 1,
                                           clustering.method = 3, debug.mode = F),
    error = function(err) {
      print(err$message)
      list()
    })
  
  meta.list.residual <- list()
  meta.list.residual <- tryCatch(
    expr  = generate.PatternScan.meta_summit_1Hz_Residual( data = data
                                                           , Resolution = 1000
                                                           , EffSize = 15
                                                           , Periodicity = 0.1
                                                           , Thres_Height = 20
                                                           , Thres_Delta = 10
                                                           , sign_flag = 2),
    error = function(err) {
      print(err$message)
      list()
    })
  
  meta.list.StandbyPower <- list()
  meta.list.StandbyPower <- tryCatch(
    expr = generate.StandbyPower.meta(),
    error = function(err) {
      print(err$message)
      list()
    })
  
  #merge list
  meta.list = list()
  meta.list = append(meta.list,meta.list.cyclic)
  meta.list = append(meta.list,meta.list.high)
  meta.list = append(meta.list,meta.list.patternHigh)
  meta.list = append(meta.list,meta.list.pattern1)
  meta.list = append(meta.list,meta.list.pattern2)
  meta.list = append(meta.list,meta.list.pattern3)
  meta.list = append(meta.list,meta.list.pattern4)
  
  if(length(meta.list)==0) return('There is no appliance')
  
  names(meta.list)=seq(0,length(meta.list)-1)
  meta.list.json <- toJSON( meta.list, pretty=TRUE )
}