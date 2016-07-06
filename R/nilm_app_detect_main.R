app.detect.main <- function( data.file,
                              find.heavy = FALSE, find.pattern = FALSE, find.pattern_extend = FALSE, find.patternHigh = TRUE,
                              find.rice = TRUE, find.standby = TRUE, find.cyclic = TRUE, find.ac = TRUE, find.washer = TRUE,
                              find.tv = TRUE,
                              check.data = TRUE, debug.mode = FALSE, printLog = TRUE) {
  ###########################

  #y = SithETL::loadCompactPowerData(data.file, sampling_rate = 15)

  if (is.character(data.file)) {
    x <- PreprocessNHz(data.file, sampling_rate = 15)
    names(x) <- c('timestamp', 'active_power', 'reactive_power')
  } else if(is.data.frame(data.file)){
    x <- data.file %>% select(timestamp, active_power, reactive_power)
  }
  x$timestamp = as.numeric(x$timestamp)
  x$timestamp = as.POSIXct(x$timestamp, tz="Asia/Seoul", origin="1970-01-01")

  y <- x
  if (check.data) {
    if (difftime(y$timestamp[max(which(y$active_power > 10))],
                 y$timestamp[min(which(y$active_power > 10))], units = 'hours') < 4){
      message('Not enough data : cyclic')
      find.cyclic <- FALSE
    }

    if (length(which(y$active_power < quantile( y$active_power, .05) + 10)) / length(y$active_power) > .5){
      message('Less than the base amount : cyclic')
      find.cyclic <- FALSE
    }

    if (length(which(y$active_power < 0 )) / length(y$active_power) > .25){
      message('Negative active power : cyclic')
      find.cyclic <- FALSE
    }

    if (length(which(y$active_power < 10 )) / length(y$active_power) > .8 ){
      message('No appliance : cyclic')
      find.cyclic <- FALSE
    }
  }

  #
  # Cyclic pattern
  #
  meta.list.cyclic <- list()
  if (find.cyclic) {
    meta.list.cyclic <- tryCatch(
      {
        #generate.CyclicBox.meta(y,max.iter=3,debug.mode=F,split.algorithm=1)
        generate.CyclicBox.meta_extend(data = y, maxLevel = 2)
      }, error = function(err) {
        list()
      })
  }


  #
  # Pattern scan
  #
  meta.list.pattern <- list()
  if ( find.pattern ) {
    meta.list.pattern <- tryCatch({
      generate.PatternScan.meta_1Hz(
        data = y, genResolution.n = 20, genEffSize.n = 15, staPeriodicity.p = 0.1,
        endEffSlot.p = 0.1, endConsistency.n = 1, clustering.method = 3, debug.mode = F)},
      error = function(err) {list()})

    #
    # Reduce resolution
    #
    if ( length(meta.list.pattern) == 0 ) {
      resolutionMin = 1 # default 1-min
      y.smoothing=smoothing(y,resolutionMin=resolutionMin)
      meta.list.pattern <- tryCatch({
        generate.PatternScan.meta_1Hz(
          data = y.smoothing, genResolution.n = 20, genEffSize.n = 15, staPeriodicity.p = 0.1,
          endEffSlot.p = 0.1, endConsistency.n = 1, clustering.method = 3, debug.mode = F)},
        error = function(err) {list()})
      meta.list.pattern <- lapply(meta.list.pattern, function(x) {
        x$generation_info$data_used$sampling <- resolutionMin * 60
        return(x)
      })
    }

    #
    # Reduce resolution one more time
    #
    if (length(meta.list.pattern) == 0) {
      resolutionMin = 3 # default 1-min
      y.smoothing=smoothing(y,resolutionMin=resolutionMin)
      meta.list.pattern <- tryCatch({
        generate.PatternScan.meta_1Hz(
          data = y.smoothing, genResolution.n = 20, genEffSize.n = 15, staPeriodicity.p = 0.1,
          endEffSlot.p = 0.1, endConsistency.n = 1, clustering.method = 3, debug.mode = F)
      }, error = function(err) {
        list()
      })
      meta.list.pattern <- lapply(meta.list.pattern, function(x) {
        x$generation_info$data_used$sampling <- resolutionMin * 60
        return(x)
      })
    }
  }

  #
  #find heavy loads
  #
  meta.list.high = list()
  if ( find.heavy ) {
    meta.list.high <- tryCatch({
      high.power.detect( y, max.iter=3,debug.mode = FALSE)
    }, error = function(err) {
      list()
    })
  }




  meta.list.pattern_extend = list()
  if ( find.pattern_extend ) {
    meta.list.pattern_extend <- tryCatch({
      generate.PatternScan.meta_1Hz_extend(
        data = y, genResolution.n = 20, genEffSize.n = 15, staPeriodicity.p = 0.1,
        endEffSlot.p = 0.1, endConsistency.n = 1, clustering.method = 3, debug.mode = F)
    }, error = function(err) {
      list()
    })

    #
    # Reduce resolution
    #
    if ( length(meta.list.pattern_extend)==0 ) {
      resolutionMin = 1 # default 1-min
      y.smoothing <- smoothing(y, resolutionMin = resolutionMin)
      meta.list.pattern_extend <- tryCatch({
        generate.PatternScan.meta_1Hz_extend(
          data = y.smoothing, genResolution.n = 20, genEffSize.n = 15, staPeriodicity.p = 0.1,
          endEffSlot.p = 0.1, endConsistency.n = 1, clustering.method = 3, debug.mode = F)
      }, error = function(err) {
        list()
      })
      meta.list.pattern_extend <- lapply(meta.list.pattern_extend, function(x) {
        x$generation_info$data_used$sampling <- resolutionMin * 60
        return(x)
      })
    }

    #
    # Reduce resolution one more time
    #
    if (length(meta.list.pattern_extend) == 0) {
      resolutionMin = 3 # default 1-min
      y.smoothing <- smoothing(y, resolutionMin = resolutionMin)
      meta.list.pattern_extend <- tryCatch({
        generate.PatternScan.meta_1Hz_extend(
          data = y.smoothing, genResolution.n = 20, genEffSize.n = 15, staPeriodicity.p = 0.1,
          endEffSlot.p = 0.1, endConsistency.n = 1, clustering.method = 3, debug.mode = F)
      }, error = function(err) {
        list()
      })
      meta.list.pattern_extend <- lapply(meta.list.pattern_extend, function(x) {
        x$generation_info$data_used$sampling <- resolutionMin * 60
        return(x)
      })
    }
  }

#  print( paste( "runtime of part 2: ", ( proc.time() - ptm)[ 3], "secs"))
#  ptm <- proc.time()


  #
  # TV
  #
  meta.list.tv <- list()
  if (find.tv) {
    meta.list.tv <- tryCatch(
      {
        list(generate.meta.tv( X = y ))
      }, error = function(err) {
        list()
      })
  }

  #
  # rice cooker
  #
  ptm <- proc.time()
  meta.list.warm <- list()
  if (find.rice) {
    meta.list.warm <- tryCatch(
      {
        # 1hz version
        # generate.PatternScan.meta.ricecooker_1Hz(data = y, eff_size = 10, periodicity =0.1, thres_ap.h = 15, thres_rp.delta = 20)
        # 15hz
        generate.PatternScan.meta.ricecooker_15Hz(data = y)
      }, error = function(err) {
        list()
      })
  }
#  print( paste( "point3:  15hz run completed", (proc.time()- ptm)[ 3]))
  print(paste( "ricecookr:", proc.time() - ptm))

  #
  # Washer
  #
  meta.list.washer = list()
  if (find.washer) {
    meta.list.washer <- tryCatch(
      {
        list(generate.meta.washer(X = y, debug = F, samplingRate = 15))
      }, error = function(err) {
        list()
      })
  }

  ###########
  ### 1 Hz
  ###########

  x2 <- convert_15Hz_to_1Hz_from_Raw(x)

  # names(x2) <- c( 'timestamp','voltage','current','active_power','reactive_power','apparent_power','power_factor')
  names(x2) <- c( 'timestamp', 'active_power', 'reactive_power' )


#  print( paste( "1hz data conversion completed", (proc.time()- ptm)[ 3]))
#  ptm <- proc.time()


  y <- x2


  #
  #find pattern heavy load
  #
  meta.list.patternHigh <- list()
  if (find.patternHigh){
    meta.list.patternHigh <- tryCatch(expr = generate.PatternScanHeavy.meta_1Hz_extend(
      y,relApDiff.max = 0.15, relRpDiff.max = 0.2, min_sigmag.active.n = 800,
      min_sigmag.reactive.n = 3, timeSpan = dhours(1), measurementErrorBound = 10,
      debug.mode = F)[[3]], error = function(err) {
        print(err$message)
        list()
      })
  }



  #
  # Air-conditioner
  #
  meta.list.ac = list()
  if (find.ac) {
    meta.list.ac <- tryCatch(
      {
        #list( generate.ac.meta(y) )
        generate.ac.meta(y, sampling = 360)
      }, error = function(err) {
        list()
      })
  }


  #
  # standby power
  #
  meta.list.standby = list()
  if (find.standby) {
    meta.list.standby <- tryCatch(
      {
        generate.StandbyPower.meta(y)
      }, error = function(err) {
        list()
      })
  }

  #
  # merge list
  #
  meta.list <- list()
  meta.list <- append(meta.list, meta.list.standby)
  if (printLog && length(meta.list.standby) > 1){
    print("Detected :: meta.list.standby")
    print(toJSON(meta.list))
  }

  if (length(meta.list.cyclic) > 0) {

    for (k in seq(meta.list.cyclic)) {
      rising_ap <- abs(meta.list.cyclic[[k]][['rising_edge']][['str2end.active_power_med'  ]])
      rising_rp <- abs(meta.list.cyclic[[k]][['rising_edge']][['str2end.reactive_power_med']])

      if ((abs(rising_ap) < 100 & abs(rising_rp) < 10) |
           any( is.na( meta.list.cyclic[[k]][['cycle']] ))) {
        #pass
      } else {
        meta.list <- append(meta.list, list(meta.list.cyclic[[k]]))
        if (printLog) {
          print("Detected :: meta.list.cyclic")
          print(toJSON(meta.list))
        }
      }
    }
  }

  if (length(meta.list.warm) > 0) {
    meta.list <- append(meta.list, list(meta.list.warm))
    if (printLog) {
      print("Detected :: meta.list.warm")
      print(toJSON(meta.list))
    }
  }

  if (length(meta.list.ac) > 0) {
    meta.list <- append(meta.list, meta.list.ac)
    if (printLog){
      print("Detected :: meta.list.ac")
      print(toJSON(meta.list))
    }
  }

  if (length(meta.list.washer) > 0) {
    meta.list <- append(meta.list, meta.list.washer)
    if (printLog){
      print("Detected :: meta.list.washer")
      print(toJSON(meta.list))
    }
  }


  if (length(meta.list.tv) > 0) {
    meta.list <- append(meta.list, meta.list.tv)
    if (printLog) {
      print("Detected :: meta.list.tv")
      print(toJSON(meta.list))
    }
  }


  for (k in seq(meta.list.pattern)) {
    rising_ap  <- abs(meta.list.pattern[[k]][['rising_edge']][['ap_med']])
    falling_ap <- abs(meta.list.pattern[[k]][['falling_edge']][['EffAP_Drop.med']])
    rising_rp  <- abs(meta.list.pattern[[k]][['rising_edge']][['rp_med']])
    falling_rp <- abs(meta.list.pattern[[k]][['falling_edge']][['EffRP_Drop.med']])
    sampling   <- meta.list.pattern[[k]][['generation_info']][['data_used']][['sampling']]

    time_on <- meta.list.pattern[[k]][['falling_edge']][['EffTimeOn.med']]
    if ((falling_ap-rising_ap)/rising_ap > 0.15
        || (rising_ap-falling_ap)/rising_ap > 0.5
        || abs(rising_rp-falling_rp)/rising_rp > 0.5
        || (min(rising_rp,falling_rp) < 15 & time_on > (3600*3) )
        || time_on < 10) {

    } else {
      meta.list = append(meta.list,list(meta.list.pattern[[k]]))
      if (printLog){
        print("Detected :: meta.list.pattern")
        print(toJSON(meta.list))
      }
    }
  }

  for (k in seq(meta.list.high)) {
    rising_ap  <- meta.list.high[[k]][['rising_edge']][['ap_height']]
    falling_ap <- -1*meta.list.high[[k]][['falling_edge']][['ap_height']]

    if (rising_ap < 400
        || abs(rising_ap-falling_ap) / rising_ap > 0.25) {

    } else {
      meta.list = append(meta.list,list(meta.list.high[[k]]))
      if (printLog){
        print("Detected :: meta.list.high")
        print(toJSON(meta.list))
      }
    }
  }

  for (k in seq(meta.list.pattern_extend)) {
    rising_ap  <- abs(meta.list.pattern_extend[[k]][['rising_edge']][['ap_med']])
    falling_ap <- abs(meta.list.pattern_extend[[k]][['falling_edge']][['ap_med']])
    rising_rp  <- abs(meta.list.pattern_extend[[k]][['rising_edge']][['rp_med']])
    falling_rp <- abs(meta.list.pattern_extend[[k]][['falling_edge']][['rp_med']])
    sampling   <- abs(meta.list.pattern_extend[[k]][['generation_info']][['data_used']][['sampling']])
    time_on    <- meta.list.pattern_extend[[k]][['falling_edge']][['EffTimeOn.med']]

    if ((falling_ap-rising_ap)/rising_ap > 0.15
        || (rising_ap-falling_ap)/rising_ap > 0.5
        || (rising_rp-falling_rp)/rising_rp > 0.5
        #|| (min(rising_rp,falling_rp) < 15 & time_on > (3600*3) )
        || time_on > (3600*3)
        || time_on < 10
        || (rising_ap<40 & rising_rp<1) ) {

    } else {
      meta.list = append(meta.list, list(meta.list.pattern_extend[[k]]))
      if (printLog){
        print("Detected :: meta.list.pattern_extend")
        print(toJSON(meta.list))
      }

    }
  }

  if (length(meta.list.patternHigh) > 0) {
    meta.list <- append(meta.list, meta.list.patternHigh)
    if (printLog){
      print("Detected :: meta.list.patternHigh")
      print(toJSON(meta.list))
    }
  }

  if (length(meta.list) > 0) {
    names(meta.list) = seq(0, length(meta.list)-1)
    meta.list.json <- toJSON(meta.list, pretty = TRUE)
    if (printLog) {
      print("Detected :: convert to JSON")
      print(toJSON(meta.list))
    }
  } else {
    meta.list.json <- I(toJSON(''))
  }
  return(meta.list.json)
}
