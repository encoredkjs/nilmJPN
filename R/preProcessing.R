preProcessing <- function(data, min_mag_h1 = 50, winSize = 15, measurementErrorBound = 10){

  rPattern <- DetectPattern_1Hz( data
                                 , position = "start"
                                 , main_type = 'active'
                                 , sub_type  = 'active'
                                 , min_mag_main = -measurementErrorBound
                                 , min_mag_sub  = -measurementErrorBound
                                 , debug.mode = F )

  rPattern  <- rPattern[ order(rPattern$start.idx), ]
  rPattern1 <- subset( rPattern, (abs(delta) < measurementErrorBound) & (h1 > min_mag_h1) )

  peakySignal.idx <- unlist( apply( rPattern1, 1, function(row)
    as.numeric(row['start.idx']):(as.numeric(row['end.idx'])+1) ), use.names=F )

  peaky.str <- data$active_power[ rPattern1$start.idx   ]
  peaky.end <- data$active_power[ rPattern1$end.idx + 1 ]
  peaky.len <- rPattern1$pattern.size + 1

  underPeak <- unlist( mapply( function(s,e,n) seq(s,e,length.out=n), peaky.str, peaky.end, peaky.len ), use.names=F )

  if( length(peakySignal.idx) != length(underPeak) ){
    stop("Error : dimention does not match (preProcessing)")
  }

  data$originalActive <- data$active_power
  data$active_power[peakySignal.idx] <- underPeak

  rPattern2 <- subset( rPattern, h1 > min_mag_h1 )
  rPattern2 <- subset( rPattern2, !(start.idx %in% rPattern1$start.idx))

  peakySignal.idx <- rbind( rPattern2$start.idx,
                            sapply( rPattern2$end.idx, function(x) pmin(x+c(0:winSize),nrow(data))) )

  tmp <- apply( peakySignal.idx, 2, function(x){ signal.diff <- data$active_power[x[1]] - data$active_power[x[-1]]
                                                 if( any(abs(signal.diff) < measurementErrorBound, na.rm=T)){
                                                   loc <- which.min(abs(signal.diff))
                                                   str <- x[1]
                                                   end <- x[-1][loc]
                                                   return(c(str,end))
                                                 }else return(NULL) })
  tmp <- tmp[!sapply( tmp, is.null )]

  peakySignal.idx <- unlist( lapply( tmp, function(row) as.numeric(row[1]):(as.numeric(row[2])) ), use.names=F )
  peaky.str <- data$active_power[sapply( tmp, function(x) x[1] )]
  peaky.end <- data$active_power[sapply( tmp, function(x) x[2] )]
  peaky.len <- sapply( tmp, function(x) x[2] - x[1] + 1 )
  underPeak <- unlist( mapply( function(s,e,n) seq(s,e,length.out=n), peaky.str, peaky.end, peaky.len ), use.names=F )

  if( length(peakySignal.idx) != length(underPeak) ){
    stop("Error : dimention does not match (preProcessing)")
  }
  data$active_power[peakySignal.idx] <- pmin( underPeak, data$active_power[peakySignal.idx] )

  return(data)
}

convert_15Hz_to_1Hz_from_Raw <- function( data, sampling_rate_raw = 15 ){

  myindex <- seq( 1, nrow(data), by = sampling_rate_raw)

  data.out <- data.frame(
    timestamp      = data$timestamp[ myindex],
    active_power   = sapply( myindex, function(x) mean( data$active_power[  seq(x, x+sampling_rate_raw -1)],
                                                       na.rm = TRUE)),
    reactive_power = sapply( myindex, function(x) mean( data$reactive_power[seq(x, x+sampling_rate_raw -1)],
                                                       na.rm = TRUE)))
  data.out
}

PreprocessNHz <- function( data.file, sampling_rate = 15){
  # data.file : filename
  data_raw <- read.csv(data.file, header = F)

  data <- data.frame( timestamp =
                        as.vector( sapply( data_raw[, 1],
                                           function(x) x + (1 / sampling_rate * (1:sampling_rate - 1)))))

  data$active_power   <- as.vector( t( data_raw[ , 1                 + 1:sampling_rate]))
  data$reactive_power <- as.vector( t( data_raw[ , sampling_rate + 1 + 1:sampling_rate]))

  # sort
  data = data[order(data$timestamp), ]
  options(digits.secs = 3)
  data
}

extractActivation = function( df, app = NULL, total = F, maxGap = 30*60, minP = 200, margin = 20 ){
  if( is.null( appliance ) ){ dataAppliance = subset( df, appliance == df$appliance[1] )
  }else dataAppliance = subset( df, appliance == app )

  if( nrow( dataAppliance ) == 0 ) stop( "extractActivation: no data" )

  dataAppliance$act = dataAppliance$active_power > minP

  myrle = rle( dataAppliance$act )
  dataAppliance$rle = do.call( "c", lapply( myrle$length, function(x) rep( x, each = x ) ) )

  mywhich = which( myrle$lengths > maxGap & !myrle$values )

  if( length( mywhich ) == 0 ){ # all activation
    out = list( dataAppliance )
  }else if( myrle$lengths[1] == nrow(dataAppliance)){
    out = list()
  }else{

    temp2 = sapply( mywhich, function(x){
      out = c()
      if( x == 1 ) out[1] = 1 else out[1] = sum( myrle$lengths[1:(x-1)] )
      out[2] = out[1] + myrle$length[x]
      out
    } )

    if( temp2[1] == 1 ) temp2 = temp2[-1] else temp2 = c( 1, temp2 )
    if( temp2[length(temp2)] == nrow( dataAppliance ) ) temp2 = temp2[ -length(temp2) ] else temp2 = c( temp2, nrow(dataAppliance) )

    temp3 = matrix( temp2, nrow = 2 )

    out = lapply( 1:ncol(temp3), function(x) dataAppliance[ max( 1, temp3[1,x] - margin):min( temp3[2,x] + margin, nrow(dataAppliance) ),] )
  }

  # output total data too
  if( total ){
    dataTotal = subset( df, appliance == "총량" )
    totalExtract = function( oneData ){

      startIndex = which.min( abs( dataTotal$timestamp - min( oneData$timestamp ) ) )
      endIndex = which.min( abs( dataTotal$timestamp - max( oneData$timestamp )  ) )

      dataTotal[startIndex:endIndex,]
    }

    out = lapply( out, function(x) rbind(x[,-c(ncol(x)-0:1)], totalExtract(x)) )
  }

  out
}

