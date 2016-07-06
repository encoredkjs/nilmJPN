#' @title edgeDetector
#' @author (author) prof. / (interpreter) choid
#' @description interpreted version of MATLAB edgeDetector

#source(file="from.Wonjong.Rhee/edgeMatching/edgeDetectorCleaner.R")
FLAT    <-  0
RISING  <-  1
FALLING <- -1
UNKNOWN <- -99


edgeDetector <- function(data, DEBUG_MODE = F) {
  ##
  # EDGEDETECTOR Summary of this function goes here
  #
  # TODO: Rising edge should be modeled as an object (or strucutre or class
  # or a container depending on language in use)
    
  ## Handle sign errors of raw data. This can happen due to installation
  #  problem.
#  if( mean(data$ap) < 0) data$ap = -data$ap 
#  if( mean(data$rp) < 0) data$rp = -data$rp 
  
  # TODO: The following should be handled in data loader.
  # TODO: No need to doulbe check it here in the future.
  #   if (length(dataTs) ~= length(dataAp)) || (length(dataTs) ~= length(dataRp))
  # 		error('Length mismatch among dataTp, dataAp, dataRp.');
  # 	end
  
  ## Define global constants (these are just nominal data; dirty coding in matlab)
  # global FLAT;
  # global RISING;
  # global FALLING;
  # global UNKNOWN;
  
  # Determine r/p/f values based on diff
  # Diff is defined as the difference between current sample and the next
  # sample (i.e. value(next) - value(current)). Because of the definition,
  # the last sample is not evaluated and remains at the default value UNKNOWN.
  SLOPE_THRESHOLD <- 3.0        # Min slope to determine rising/flat/falling segments.
  # This is defined as an absolute value of neighboring two samples.
  # If needed, convert it to reflect the value of time difference.
  MAX_ALLOWED_TIME_DIFF_THRESHOLD <- 2.0  # Tag a sample only if time difference is less than this value.
  if( as.numeric( min(diff(data$ts)), units='secs' ) > MAX_ALLOWED_TIME_DIFF_THRESHOLD ){
    MAX_ALLOWED_TIME_DIFF_THRESHOLD <- as.numeric( min(diff(data$ts)), units='secs' ) * 2
  }
  # Does not need to be the same as sampling period
  
  ap.diff <- diff(data$ap)
  rp.diff <- diff(data$rp)
  
  ## Per sample tags
  DATA_LENGTH <- nrow(data)
  tagsAp <- rep(UNKNOWN, DATA_LENGTH)
  tagsRp <- rep(UNKNOWN, DATA_LENGTH)
  
  tagsAp[which( ap.diff <= -SLOPE_THRESHOLD )]    <- FALLING
  tagsAp[which( ap.diff >=  SLOPE_THRESHOLD )]    <- RISING
  tagsAp[which( abs(ap.diff) < SLOPE_THRESHOLD )] <- FLAT
  
  tagsRp[which( rp.diff <= -SLOPE_THRESHOLD )]    <- FALLING
  tagsRp[which( rp.diff >=  SLOPE_THRESHOLD )]    <- RISING
  tagsRp[which( abs(rp.diff) < SLOPE_THRESHOLD )] <- FLAT

  for( i in which( diff(data$ts) > MAX_ALLOWED_TIME_DIFF_THRESHOLD ) ){
    tagsAp[i] <- tagsAp[i - 1]
    tagsRp[i] <- tagsRp[i - 1]
  } 
  
#   tagsAp[ which( diff(data$ts) > MAX_ALLOWED_TIME_DIFF_THRESHOLD ) ] <- 
#     tagsAp[ which( diff(data$ts) > MAX_ALLOWED_TIME_DIFF_THRESHOLD ) - 1 ]
#   tagsRp[ which( diff(data$ts) > MAX_ALLOWED_TIME_DIFF_THRESHOLD ) ] <- 
#     tagsRp[ which( diff(data$ts) > MAX_ALLOWED_TIME_DIFF_THRESHOLD ) - 1 ]
  
  
  # %% DEBUG_MODE = true? Then plot a figure.
  # if DEBUG_MODE
  # figure(1)
  # plot(dataTs, dataAp,'b--'); hold on
  # plot(dataTs, dataRp,'c--');
  # hold off
  # 
  # figure(2)
  # plot(dataTs, dataAp,'b--'); hold on
  # plot(dataTs(tagsAp==RISING), dataAp(tagsAp==RISING),'r^')
  # plot(dataTs(tagsAp==FALLING), dataAp(tagsAp==FALLING),'kv')
  # plot(dataTs(tagsAp==FLAT), dataAp(tagsAp==FLAT),'go')
  # plot(dataTs(tagsAp==UNKNOWN), dataAp(tagsAp==UNKNOWN),'y*')
  # 
  # plot(dataTs, dataRp,'c--');
  # plot(dataTs(tagsRp==RISING), dataRp(tagsRp==RISING),'m^')
  # plot(dataTs(tagsRp==FALLING), dataRp(tagsRp==FALLING),'kv')
  # plot(dataTs(tagsRp==FLAT), dataRp(tagsRp==FLAT),'go')
  # plot(dataTs(tagsRp==UNKNOWN), dataRp(tagsRp==UNKNOWN),'y*')    
  # hold off;
  # end
  
  ## Examine samples in sequence to take care of exceptions.
  # Update tagsAP and tagsRP based on the results.
  # Exception 1: Height of a RISING sequence is not large enough. Convert the
  # sequence to FLAT.
  # Exception 2: Height of a FALLING sequence is not large enough. Convert the
  # sequence to FLAT.
  # Exception 3: Handle one-off FLAT sample between RISING and FALLING.
  # Convert to either RISING or FALLING depending on the slope.
  
  tagsAp <- edgeDetectorCleaner(tagsAp, data$ap, DEBUG_MODE)
  tagsRp <- edgeDetectorCleaner(tagsRp, data$rp, DEBUG_MODE)
  
  ## Sequential-analysis to detect rising-edges, falling edges 
  # Generate signaure vectors as a by-product.
  # Create a structure for summarizing section information.
  # rising : RISING - FALLING - FLAT signatures
  # falling : FALLING - RISING - FLAT signatures
  
  risingEdgeList  <- edgeDetector_sub(data, tagsAp, tagsRp, 'rising' )
  fallingEdgeList <- edgeDetector_sub(data, tagsAp, tagsRp, 'falling')	
  
  return(list(risingEdgeList, fallingEdgeList))
}



## subfunction
edgeDetector_sub <- function(data, tagsAp, tagsRp, mode = 'rising') {
  ## mode -> 'rising'(default) or 'falling'
  
  # Note: RISING, FALLING, FLAT, UNKNOWN are defined as global 
  # parameters in edgeDetector.m.
  
  ## mode setup
  # rising edge detector : detect rising -> falling -> flat 
  # falling edge detector : detect falling -> rising -> flat
  
  if( mode == 'falling' ) edge.seq <- c('first'=FALLING, 'second'=RISING,  'third'=FLAT)
  if( mode == 'rising'  ) edge.seq <- c('first'=RISING,  'second'=FALLING, 'third'=FLAT)
    
  ## Sequential-analysis to detect edges 
  # Generate signaure vectors as a by-product.
  # Create a structure for summarizing section information.
  
  # Initially set matrix (structure) height to be max and cut it later.
  # (For MATLAB memory mangement)
  # NUM_ATTRIBUTES <- 12
  # lengthTags <- length(tagsAp)
  # 	edgeList <- data.frame(matrix(0, ncol = NUM_ATTRIBUTES, nrow = lengthTags))
  # numEdges <- 0
  
  ## A risingEdge is defined as a sequence of RISING-FALLING-FLAT. FALLING is
  ## optional (RISING-FALLING). In the output, we include RISING-UNKNOWN as
  ## well for analysis. (TODO: Consider removing RISING-UNKNOWN)
  ## Other sequence starting with RISING (e.g. RISING-FALLING-RISING, RISING-UNKNOWN)
  ## is not considered as a risingEdge.
  # TODO: Currently, all the decisions are based on tagsAp only. Consider
  # TODO: tagsRp for boundary decisions.
  
  tagsAp.rle <- rle(tagsAp)
  str.ind <- c(0,cumsum(tagsAp.rle$length[-length(tagsAp.rle$length)])) + 1
  end.ind <- cumsum(tagsAp.rle$length) + 1
  tagsAp.rle[['Ap.height']] <- data$ap[end.ind] - data$ap[str.ind]
  tagsAp.rle[['Rp.height']] <- data$rp[end.ind] - data$rp[str.ind]

  tagsAp.rle[['steadyStateHeightAp']] <- tagsAp.rle[['Ap.height']]
  tagsAp.rle[['steadyStateHeightRp']] <- tagsAp.rle[['Rp.height']]
  for( k in which( tagsAp.rle$values == edge.seq['first'] )){
    if( tagsAp.rle$values[k+1] == edge.seq['second'] ){
      tagsAp.rle[['steadyStateHeightAp']][k:(k+1)] <- sum(tagsAp.rle[['Ap.height']][k:(k+1)])
      tagsAp.rle[['steadyStateHeightRp']][k:(k+1)] <- sum(tagsAp.rle[['Rp.height']][k:(k+1)])
    }
  }
  tagsAp.rle <- data.frame(unclass(tagsAp.rle))
  
  calc.edge.height.length <- function(k){
    tmp <- c( str.ind[k], 
              tagsAp.rle[['steadyStateHeightAp']][k], 
              tagsAp.rle[['steadyStateHeightRp']][k],
              tagsAp.rle$lengths[k:(k+2)],
              tagsAp.rle$Ap.height[k:(k+2)],
              tagsAp.rle$Rp.height[k:(k+2)] )                
    
    if( (k+1) <= nrow(tagsAp.rle) ){
      if( tagsAp.rle$values[k+1] != edge.seq['second'] ){
        tmp[c(6,9,12)] <- tmp[c(5,8,11)]
        tmp[c(5,8,11)] <- 0
      }
    }else tmp[c(5:6,8:9,11:12)] <- 0
    
    if( (k+2) <= nrow(tagsAp.rle) ){
      if( (tagsAp.rle$values[k+1] == edge.seq['second']) & 
            (tagsAp.rle$values[k+2] != edge.seq['third'])  ) tmp[c(6,9,12)] <- 0 
    }else tmp[c(6,9,12)] <- 0
    
    return(tmp)
  }
  
  edgeList <- lapply( which( tagsAp.rle$values == edge.seq['first']), calc.edge.height.length )
  edgeList <- ldply(edgeList)	
  
  return(edgeList)
  
  # This part may not be needed for R
  # Truncate
  # edgeList = edgeList(1:numEdges, :); 
}