
hourlyF = function(df, sampling = 15, col = "p"){
  out <- as.list(rep(0, 24))
  names(out) <- as.character(0:23)

  df <- df[!is.na(df$timestamp),]
  out.temp <- dlply(df, .(hour = lubridate::hour(timestamp)),
                    function(x){ sum(x[[col]], na.rm = TRUE) / 3600 / sampling})
  # to put hourly into the right place
  out[match(names(out.temp), names(out))] <- out.temp
  return(out)
}

app.usage.main <- function(data.file, app.meta.json, find.heavy = F,
                           find.heavy.pattern = T, find.pattern = T, find.rice = T,
                           find.cyclic = T, find.ac = F, find.washer=T, find.standby = T,
                           find.tv = T, show.fig = F){

  # #########################################################################
  #
  # data.file & app.meta.json should be defiend
  #
  # ##########################################################################

  if (is.character(data.file)){
    # x = file2data( data.file)
    #y = SithETL::loadCompactPowerData(data.file, sampling_rate = 15)

    y <- PreprocessNHz(data.file, sampling_rate = 15)

    names(y) <- c('timestamp', 'active_power','reactive_power')

    time.diff <- c(diff(y$timestamp), 999)
    y <- y[time.diff != 0,]
    # from file2data function (removed because this takes awhile)
    # time.diff <- diff(x$timestamp )
    # for (i in unique(x$timestamp[time.diff == 0])){
    #  x <- rbind(subset(x, timestamp != i),
    #             data.frame(lapply(subset(x, timestamp == i), function(col) median(col))))
    # }

    y$timestamp <- as.POSIXct(y$timestamp, tz = "Asia/Seoul", origin = "1970-01-01")

  } else if (is.data.frame(data.file)) {
    y = data.file
  }

  #
  #
  # ---------------------- Start Main Code ----------------------------------
  #
  #

  app.meta.all <- fromJSON(app.meta.json)
  app.usage <- list()
  app.usage.res <- list()
  app.usage.res[['daily' ]] <- sum(y$active_power) / 3600 / 15.
  app.usage.res[['hourly']] <- hourlyF(y, sampling = 15, col="active_power")

  if (length(app.meta.all) > 0) {

    app.meta.cyclic <- app.meta.all[sapply(app.meta.all, function(x) x$shape_type == "cyclic_box")]
    app.meta.noncyclic.15hz <- app.meta.all[sapply(app.meta.all, function(x) x$shape_type %in%
                                                     c("TV", "ricecooker_pattern_scan", "WashingMachine"))]
    app.meta.noncyclic.1hz  <- app.meta.all[sapply(app.meta.all, function(x) x$shape_type %in%
                                                     c('StandbyPower', "pattern_scan_heavy", "HMM"))]

    if (find.cyclic & length(app.meta.cyclic) > 0) {
      app.usage.cyclic <- multipleMeta2cyclicBoxExtend(meta = app.meta.cyclic, data = y, postProcessing = T)
      names(app.usage.cyclic) <- names(app.meta.cyclic)
      print(seq_along(app.meta.cyclic))
      for (iapp in names(app.usage.cyclic)){
        if (length(app.usage.cyclic[[iapp]]) != 0){
          temp <- app.usage.cyclic[[iapp]]
        } else {
          temp <- x
          temp$p <- 0
        }
        app.usage[[iapp]] = list()
        app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm = TRUE) / 3600 / 15
        app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp)
        app.usage[[iapp]][['shape_type']] <- paste("cyclic", iapp)

        app.usage.res[['daily']] = max(app.usage.res[['daily']] - app.usage[[iapp]][['daily']], 0)
        app.usage.res[['hourly']] <- as.list(
          tryCatch(sapply(as.character(0:23), function(x)
            max(0, app.usage.res[['hourly']][[x]] - app.usage[[iapp]][['hourly']][[x]])),
                   error = function(e) {
                     resCatch = rep(0, 24)
                     names(resCatch) <- 0:23
                     return(resCatch)
                   }))
        #j = j + 1
      }
    }

    for (iapp in names(app.meta.noncyclic.15hz)){
      print(iapp)
      app.meta = app.meta.noncyclic.15hz[[iapp]]
      shape_type = app.meta[['shape_type']]

      if (!((shape_type == "TV" & find.tv) |
            (shape_type == "ricecooker_pattern_scan" & find.rice)|
            (shape_type == "WashingMachine" & find.washer))) {
        #pass
      } else {
        app.usage[[iapp]] = list()

        if (shape_type == "TV" & find.tv) {
          temp <- tryCatch({meta2tv(meta = app.meta, data = y)},
                           error = function(err) { data.frame(timestamp = y$timestamp, p = 0)})
          app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm = TRUE) / 3600 / 15.
          app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp, sampling = 15)
          app.usage[[iapp]][['shape_type']] <- shape_type

        } else if (shape_type == "ricecooker_pattern_scan" & find.rice) {
          temp <- tryCatch({ meta2PatternScan.ricecooker_15Hz(data = y, app.meta) },
                           error = function(err) {data.frame(timestamp = y$timestamp, p = 0) })
          if (show.fig) {print(head(temp))}
          app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm = TRUE) / 3600. # make sure output format of ricecooker_15Hz
          app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp, sampling = 1)
          app.usage[[iapp]][['shape_type']] <- shape_type
        } else if (shape_type == "WashingMachine" & find.washer) {
          temp <- tryCatch(predict.washer(meta = app.meta, data = y, debug = show.fig),
                           error = function(err) {
                             data.frame(timestamp = y$timestamp, p = 0)
                            }
                           )
          if (show.fig) { print(head(temp)) }
          app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm = TRUE) / 3600 / 15.
          app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp, sampling = 15)
          app.usage[[iapp]][['shape_type']] <- shape_type
        }

        app.usage.res[['daily' ]] <- max(0, app.usage.res[['daily']] - app.usage[[iapp]][['daily']])
        app.usage.res[['hourly']] <- as.list(
          tryCatch(sapply(as.character(0:23), function(x)
            max(0, app.usage.res[['hourly']][[x]] - app.usage[[iapp]][['hourly']][[x]])),
                   error = function(e) {
                                        resCatch = rep(0, 24)
                                        names(resCatch) <- 0:23
                                        return(resCatch)
                                      }))
      }
    }

    ##########################
    #  Downscale algorithm
    ##########################

    x2 <- convert_15Hz_to_1Hz_from_Raw(y)
    x2 <- x2[order(x2$timestamp), ]
    y  <- x2

    for (iapp in names(app.meta.noncyclic.1hz)) {
      print(iapp)
      app.meta <- app.meta.noncyclic.1hz[[iapp]]
      shape_type <- app.meta[['shape_type']]

      if (!((shape_type == 'StandbyPower' & find.standby) |
            (shape_type == "HMM" & find.ac)  |
            (shape_type == "pattern_scan_heavy" & find.heavy.pattern))) {
        #pass
      } else {

        app.usage[[iapp]] <- list()

        if (shape_type == 'StandbyPower' & find.standby) {
          #base_power <- tryCatch({calculateStandbyPower(y)/nrow(y) }, error = function(err) { 0 })
          base_power <- app.meta[['value']]
          temp <- data.frame(timestamp = y$timestamp, p = base_power)
          if (show.fig) { print(head(temp)) }
          app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm = TRUE) / 3600.
          app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp, sampling = 1)
          app.usage[[iapp]][['shape_type']] <- shape_type

        } else if (shape_type == "HMM" & find.ac) {
          temp <- tryCatch({ meta2ac(data = y, app.meta, hsmm = T) },
                           error = function(err) { data.frame(timestamp = y$timestamp, p = 0) })
          if (show.fig) { print(head(temp)) }
          app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm=TRUE) / 3600.
          app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp, sampling = 1)
          app.usage[[iapp]][['shape_type']] <- shape_type

        } else if (shape_type == "pattern_scan_heavy" & find.heavy.pattern) {
          temp <- tryCatch({ meta2PatternScanHeavy_1Hz(data = y, meta = app.meta, extension.p = 0.2,
                                                       relApDiff.max = 0.15, relRpDiff.max = 0.2,
                                                       measurementErrorBound = 10,
                                                       debug.mode = show.fig) }, error = function(err) { data.frame(timestamp=y$timestamp, p=0) })
          if (show.fig) { print(head(temp)) }
          app.usage[[iapp]][['daily'     ]] <- sum(temp$p, na.rm = TRUE) / 3600.
          app.usage[[iapp]][['hourly'    ]] <- hourlyF(temp, sampling = 1)
          app.usage[[iapp]][['shape_type']] <- shape_type

        } 

        app.usage.res[['daily']]  <- max(app.usage.res[['daily']] - app.usage[[iapp]][['daily']], 0)
        app.usage.res[['hourly']] <- as.list(
          tryCatch(sapply(as.character(0:23), function(x) max(0, app.usage.res[['hourly']][[x]] - app.usage[[iapp]][['hourly']][[x]])),
                   error = function(e) {
                     resCatch = rep(0, 24)
                     names(resCatch) <- 0:23
                     return(resCatch)
                     }))
      }
    }
  }
  app.usage[['999']] <- app.usage.res
  app.usage[['999']][['shape_type']] <- "unknown"
  return(app.usage)
}

# obsolete
file2data <- function(data.file){

  x <- data.table::fread(data.file, showProgress = F, colClasses = rep("numeric", 3), header=F, sep=',')
  names(x) <- c('timestamp', 'active_power', 'reactive_power')

  time.diff <- diff(x$timestamp)
  for (i in unique(x$timestamp[time.diff == 0])){
    x <- rbind(subset(x, timestamp != i),
               data.frame(lapply(subset(x, timestamp == i), function(col) median(col))))
  }

  x <- x[order(x$timestamp), ]
  x$timestamp <- as.POSIXct(x$timestamp, tz="Asia/Seoul", origin="1970-01-01")

  return(x)
}
