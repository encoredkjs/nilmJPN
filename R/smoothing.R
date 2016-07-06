smoothing <- function( data, resolutionMin ){
  data$timestamp <- floor_date(data$timestamp,'hour') + 
    dminutes((minute(data$timestamp) %/% resolutionMin) * resolutionMin)
  data <- ddply( data, .(timestamp), summarize, 
                 active_power = median(active_power), 
                 reactive_power = median(reactive_power))
  return(data)
}

data2LowResolutionMeta <- function(data,resolutionMin){
  library(RJSONIO)
  data <- smoothing( data, resolutionMin )
  meta <- data2meta( data )
  meta <- fromJSON(meta)
  meta <- lapply( meta, function(x){x$generation_info$data_used$sampling <- resolutionMin * 60; return(x)})
  return( toJSON(meta) )
}

LowResolutionMeta2signal <- function( meta, data, postprocessing=T, show.fig=T, extension.p =0.2, flexibility.p = 0.9 ){
  
  sampling <- meta$generation_info$data_used$sampling
  data.original <- data
  ts <- data.original$timestamp
  
  if( sampling != 1 ){
    resolutionMin <- as.integer(sampling / 60)
    data.original$timestamp <- floor_date(data.original$timestamp,'hour') + 
      dminutes((minute(data.original$timestamp) %/% resolutionMin) * resolutionMin)
    data   <- smoothing( data, resolutionMin )
  } 
  
  signal <- meta2signal( meta, data, postprocessing, show.fig, extension.p, flexibility.p)
  data.original <- merge( data.original, signal, all.x = T )
  data.original$timestamp <- ts
  
  return( data.original[,c('timestamp','p','q')] )
}

