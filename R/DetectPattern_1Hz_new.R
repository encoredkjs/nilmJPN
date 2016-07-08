DetectPattern_1Hz_new <- function ( data, position = c("start", "end")
                                    , main_type = c("reactive", "active","-reactive", "-active")
                                    , sub_type  = c("active", "reactive","-active", "-reactive")
                                    , c_factor = 0
                                    , debug.mode = FALSE
                                    , use_Ap = FALSE
                                    , useConsistencyFlag = TRUE ){

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
      title(paste('Step1 : rising edge detection result\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
    }else{
      title(paste('Step1 : falling edge detection result\n', min(data$timestamp), '--', max(data$timestamp)), 
            outer=TRUE)
    }
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }
  
  return(resultant[ order(resultant$start.timestamp),] )
}
