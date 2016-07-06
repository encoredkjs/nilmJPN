#library( lubridate )
#library( plyr )
#library( dplyr )
#library( depmixS4 )
#library( Rsolnp )
#library( hsmm )
#library( tvd ) # move into function?


meta2ac <- function( data, 
                     meta.json = NULL, dimension = 2, sampling = 1,
                     par.free = F, min.diff = 800, min.length = 600,
                     tvd = F, hsmm = T, tuning = c( .8, .9, 0 ), fast = F ){
  
  data2 = data[seq(1, nrow(data), by = sampling), ]
  
  # RUN
  hmmout = suppressWarnings( 
             hmmNoMeta( data = data2, par.free = par.free, hsmm = hsmm, 
                        dimension = dimension, min.diff = min.diff, 
                        tvd = tvd, fast = fast ))

  # post
  hmmout2 = state_post( hmmout, hsmm = hsmm, min.length = min.length/sampling )
  
  state2pq( hmmout2, originalData = data, tuning = tuning )
}


generate.ac.meta <- function( data, 
                             byMonth = T, sampling = 60, dimension = 1, 
                             par.free = F, min.diff = 500, 
                             min.length = 300, tvd = T, hsmm = F){
  out = NULL
  
  # sampling
  data = data[seq(1,nrow(data), by=sampling), ]

  # meta generation rule
  myCondition = TRUE
  if ( byMonth ){
    mymonth = unique( lubridate::month( data$timestamp ) )
    myCondition = any( mymonth %in% c(6, 7, 8)) # meta generate only these months
  } 
  if ( myCondition ){
    
    # RUN
    hmmout = suppressWarnings( 
      hmmNoMeta( data = data, 
                 par.free = par.free, hsmm = hsmm, dimension = dimension, 
                 min.diff = min.diff, tvd = tvd  ))
    
    # meta record (1/2)
    if ( hsmm ){ 
          out$method = "hsmm" 
    }else out$method = "depmixS4"  
    
    temp = segmentsSummary( hmmout$state )
    out$segment_num_hmm    = temp$num
    out$segment_mediansec_hmm = temp$median * sampling    
    
    
    # record (2/2)
    out$shape_type = "HMM"
    out$'meta-version' = 2
    out$generation_info = 
      list( data_used = list( start = min( data$timestamp ), 
                              end = max( data$timestamp ),
                              sampling = sampling), 
            computed = Sys.time() )
    # out$init.param    = init.param
  }  
  if ( is.null( out ) ){ 
    list()
  }else{
    list( out ) 
  }
}


hmmNoMeta <- function( data, par.free = F, dimension = 1, min.diff = 500, sampling = 1, tvd = F, hsmm = F, fast=F ){
  # data = testData; par.free = F; dimension = 1; min.diff = 500; windows = 60; sampling = 1; tvd= T; hsmm=F
  data = data[ order( data$timestamp ), ]
  data = data[ seq( 1, nrow(data), by = sampling ), ]
  
  
  mydate = unique( floor_date( data$timestamp, "day" ) )
  out = c()
  
  if ( tvd ){
    data$active_power =  tvd1d( data$active_power, lambda = 2000 )
  }   
  for ( i in mydate ){
    tempData = data[ floor_date( data$timestamp, "day" ) == i, ]
    
    # hmm run condition 
    if ( max( tempData$active_power, na.rm=T ) - min( tempData$active_power, na.rm = T ) < min.diff ){
      tempData$state <- 1 
    }else{
      ##############
      # hmm call
      ##############  
      tryCatch({
        meta.json <- list( prior = c(.99, .01), transition = c( .99, .01, .01, .99 ), 
                           APEmissionDiff = c( 850, 600 ), RPEmissionDiff = c( -10, 100 ) )
        
        tempData  = data2ac_hmm( meta.json = meta.json, data = tempData, par.free = par.free, 
                                 dimension = dimension, hsmm= hsmm, fast = fast )
      }, error = function(e){
        tempData$state <- 1
      })
    }
    
    # original state
    out = rbind( out, tempData )
    
    ## parameter update : need to refine
  }
  out
}

####

state2pq <- function( data, originalData = NULL, min.diff = 400, tuning = c( .8, .9, 0.6 )  ){
  # data = hmmout2; min.diff = 500; tuning = c( .9, .99, .6 )
  
  # center predict for all state call 
  
  data$p <- tuning[3]
  data$q <- 0
  
  # floor
  if ( length( unique( data$state)) == 2 ){ 
    data$p[ data$state == 2 ] <- 
    pmax( tuning[3], ( data$active_power[data$state==2]  -
                        median( data$active_power[data$state==1], na.rm = T ) ) * tuning[1] )
  }
  
  data$state[ data$p < min.diff ] <- 1 
  data$p[     data$p < min.diff ] <- tuning[3]
  
  
  # ceiling by quantile
  pvec = data$p[data$state == 2 ]
  if ( length( pvec ) > 0 ) data$p <- pmin( data$p, quantile( pvec, tuning[2], na.rm=T ) ) 
  
  # q
  data$q[ data$state == 2 ] <- data$reactive_power[data$state==2] -
    median( data$reactive_power[data$state==1], na.rm=T )
  
  # upscale
  if ( is.null( originalData ) ) originalData = data
  originalData$p = approx( data$timestamp, data$p, xout = originalData$timestamp, rule = 2 )$y
  originalData$q = approx( data$timestamp, data$q, xout = originalData$timestamp, rule = 2 )$y
  
  originalData
}



segmentsSummary <- function( vec, value = 2 ){
  myrle     = rle( vec )
  mylengths = myrle$lengths[ myrle$values == value ] 
  
  list( num = length( mylengths ), median = median( mylengths ) )
}

notReliableF <- function( data ){
  durationRate = mean( data$state == 2, na.rm = T )
  ( durationRate < 30/60/60/24 | durationRate > 0.5 )
}

overshootF <- function( data, rawdata, windows = 10 ){
  data$diffState = c( NA, diff( data$state ) )

  rawdata$diff   = c( NA, abs( diff( rawdata$reactive_power )))

  myts = data$timestamp[data$diffState == 1]

  data.frame( timestamp = myts, 
              diff = sapply( myts, 
                             function(x) max( rawdata$diff[ (rawdata$timestamp >=  x - windows)  &
                                               (rawdata$timestamp <=  x + windows)    ], na.rm=T ) ) )
}


data2ac_hmm <- function( meta.json = NULL, data, hsmm = F, par.free = F, fast = F, 
                        dimension=2, min.diff = 400 ){
    if ( hsmm == T ){ 
      data$state = detectHSMM( data = data, meta = NULL )$hsmmv
    }else{ 
      depmixout = detectDepmixS4( data = data, AC = meta.json, dimension = dimension, 
                                  par.free = par.free, fast = fast, hsmm = hsmm )
      data$state = depmixout$state
    }
    data    
}

#####
# newer version of detectAndpostprocess
detectDepmixS4 <- function( data, AC, 
                           dimension = 2, transition = F, par.free = T, hsmm = F, fast = F ){
  if ( par.free ){
    out = detect( data = data, AC=AC, dimension=dimension, transition = transition, APFixed = c(F,F,F,F) )
  }else if ( fast ){
    out = detect( data = data, AC=AC, dimension=dimension, transition = transition, APFixed = c(T,T,T,T),
                  RPFixed = c( T,T,T,T ), fixTransition = T )
  }else{
    # fix the difference rather than means themselves
    MeanFixed = c( T, F )
    SdFixed   = c( F, F )
    transFix  = F
    
    out = detectC( data=data, AC=AC, dimension=dimension, transition=transition, MeanFixed=MeanFixed, SdFixed=SdFixed, 
                   transFix = transFix )
  }
  state = out$fm@posterior$state
  list( state = state, mod = out$mod, fm = out$fm, data = data )  # df with only state
}


state_post <- function( data, min.diff = 400, min.length = 60, hsmm = F ){
  # state ordering 
  tempData = setACto2( data )
  
  if ( hsmm ) tempData = smooth.tvdf( tempData )
  
  # exclude Short    
  tempData$state = excludeShort( tempData$state, min.length = min.length )
  
  tempData
}



smooth.tvdf <- function( data, lambda = 1000 ){
  data$state = round( tvd1d( data$state, lambda=lambda ) )
  class( data ) <- c( class( data), "commonModels" )
  data
}


filterAC <- function( data, windows = 5, minO = 400 ){
  data$stateDiff = c( NA, diff( data$state )           )
  data$AbsRDiff  = c( NA, diff( data$reactive_power )  )
    
  startIndex = which( data$stateDiff ==  1 )
  endIndex   = which( data$stateDiff == -1 ) - 1

  if ( startIndex[1] > endIndex[1] ) endIndex = endIndex[-1]

  if ( startIndex[ length(startIndex) ] > endIndex[ length(endIndex) ] ) endIndex <- c( endIndex, nrow(data) ) 
  
  indexDF = data.frame( startIndex=startIndex, endIndex = endIndex )
  
  
  tempqq = sapply( 1:nrow(indexDF), 
                   function(x){
                     max( data$AbsRDiff[ max(1, 
                       indexDF[x,'startIndex'] - windows):min(nrow(data), indexDF[x,'startIndex'] + windows )], na.rm=T ) 
                     })
  indexDF_X = indexDF[ tempqq < minO, ]

  if ( nrow( indexDF_X ) != 0 ){
    for( i in 1:nrow( indexDF_X ) ){
      data$state[ indexDF_X[i,"startIndex"]:indexDF_X[i,"endIndex"]  ] <- 1   
    }
  } 
  data[,!( colnames( data ) %in% c("stateDiff", "AbsRDiff"))]
}


setACto2 <- function( data ){
  mymed = sapply( c(1,2), function(x) median( data$active_power[data$state == x], na.rm=T ) )
    if ( mymed[1] > mymed[2] ) data$state <- 3 - data$state
  
  data
}




parUpdate <- function( hmmout, AC, weight = 0.5, min.diff = 350 ){
  durationRate = mean( hmmout$state == 2, na.rm = T )
  
  newAC = list()
  
  if ( durationRate < 30/60/60/24 | durationRate > 0.5 ){ # no update condition
    newAC = AC
  }else{
    temp <- hmmout[ ,c( "active_power", "reactive_power", "state" )] 
    
    myEst <- ddply( temp, "state", summarise, meanA = mean( active_power ), sdA = sd( active_power ), 
                    meanR = mean( reactive_power ), sdR = sd( reactive_power )  )
    
    myDiff = myEst[2,-1] - myEst[1,-1]
    
    if ( myDiff[1] < min.diff ){
      newAC = AC
    }else{
      newPars = list( APEmissionDiff = unlist( myDiff )[1:2], 
                      RPEmissionDiff = unlist( myDiff )[3:4]  )
      
      newAC = lapply( c("APEmissionDiff", "RPEmissionDiff"), 
                      function(x) AC[[x]]*( 1-weight ) + newPars[[x]]*weight )
      names( newAC ) <- c("APEmissionDiff", "RPEmissionDiff")
    }
  }
  newAC
}


detectHSMM <- function( data, meta = NULL ){
  # data = evalData4; meta = NULL; min.diff=400
  if ( is.null( meta ) ){
    meta = list()
    meta$prior = c(.90, .10 )
    meta$transition = c( .999, .001, .01, .99 )
    meta$APEmissionDiff= c( 850, 500 )
  }
  
  # emission
  temp = subset( data, active_power < quantile( active_power, .75 ) )
  bgParam = c( mean( temp$active_power,   na.rm=T ), sd( temp$active_power,   na.rm=T ) )
  
  pardf2 = c( 0, 0, meta$APEmissionDiff ) + rep( unlist( bgParam[1:2] ), 2 )  
  
  pipar  <- meta$prior
  tpmpar <- t( matrix( meta$transition, nrow=2 ) )
  # rdpar  <- c( 100, .3 )
  odpar  <- list( mean = pardf2[c(1,3)], var = pardf2[c(2,4)]^2 )
  
  # run length distribution 
#  rd = "nbinom"; rdpar = list( r = c(1E5, 1E5), pi = c(0.1, 0.1) )
#    rd = "log";  rdpar = list( p = c(.4, .4) ) 
    rd = "pois"; rdpar = list( lambda = c( 1e-6, 1e-3) )
  
  #  myhsmm  = hsmm(         data$active_power, od="norm", rd = rd, rd.par = rdpar, 
  #                          pi.par = pipar, tpm.par = tpmpar, od.par = odpar       ) 
  myhsmmv = hsmm.viterbi( data$active_power, od="norm", rd = rd, rd.par = rdpar, 
                         pi.par = pipar, tpm.par = tpmpar, od.par = odpar        )

  data.frame( hsmmv = myhsmmv$path )  
}

detect <- function( data, AC, 
                   dimension = 2, transition = F, APFixed = c( T, F, T, F ), 
                   RPFixed = c( F, F, F, F ), fixPrior = F, fixTransition = F ){
  tempData = data
  pardf1 = data.frame( par = c( AC$prior, AC$transition ), fixed = c(rep(fixPrior, 2), 
                       rep( fixTransition, 4 + 4*transition) ))
  
  # emission
  temp = subset( tempData, active_power < quantile( active_power, .75 ) )
  bgParam = c( meanA = mean( temp$active_power,   na.rm=T ), sdA = sd( temp$active_power,   na.rm=T ), 
               meanR = mean( temp$reactive_power, na.rm=T ), sdR = sd( temp$reactive_power, na.rm=T ) )
    
  pardf2 = data.frame( par = c( 0, 0, AC$APEmissionDiff ) + 
                         rep( unlist( bgParam[1:2] ), 2 ), fixed = APFixed )
  
  if ( dimension==2 ){
    pardf3 = data.frame( par = c( 0, 0, AC$RPEmissionDiff ) + 
                           rep( unlist( bgParam[3:4] ), 2 ), fixed = RPFixed )
    pardf2 = rbind( pardf2, pardf3 )[c(1, 2, 5, 6, 3, 4, 7, 8), ]
  }
  
  pardf = rbind( pardf1, pardf2 )
  
  # transition modeling
  transitionForm = as.formula( "~1" )
  if ( transition ){
    tempData$diff_reactive_power = c( diff( tempData$reactive_power ), NA )
    tempData$diff_reactive_power = max( 0, tempData$reactive_power )
    transitionForm = as.formula( "~ diff_reactive_power" )
  }

  if ( dimension == 1 ){
    mod <- depmix( response = active_power ~ 1,  
                   transition = transitionForm, data = tempData, nstates = 2 )
  }else{
    mod <- depmix( response = list( active_power ~ 1, reactive_power ~ 1 ), 
                   family = list( gaussian(), gaussian() ), 
                   transition = transitionForm, data = tempData, nstates = 2 )
  }
  newpars = pardf$par # getpars( mod )

  mod = setpars( mod, newpars )
  
  fm = depmixS4::fit( mod, fixed = pardf$fixed )
  
  list( fm = fm, data = tempData, mod = mod, computed = Sys.time() )  
}

# detectC uses conrows argument for more flexible modeling compared to detect
detectC <- function( data, AC, dimension = 2, transition = F, MeanFixed = c( T, F ), SdFixed = c( F, F ), transFix = F ){

  tempData = data
  
  # background
  pardf1 = data.frame( par = c( AC$prior, AC$transition ), fixed = c( F, F, seq( transFix, along.with=AC$transition ) ))
  
  # emission
  temp = subset( tempData, active_power < quantile( active_power, .75 ) )
  bgParam = c( meanA = mean( temp$active_power,   na.rm=T ), sdA = sd( temp$active_power,   na.rm=T ), 
               meanR = mean( temp$reactive_power, na.rm=T ), sdR = sd( temp$reactive_power, na.rm=T ) )
  
  pardf2 = data.frame( par = c( 0, 0, AC$APEmissionDiff ) + 
                         rep( unlist( bgParam[1:2] ), 2 ), fixed = NA )
  
  if ( dimension==2 ){
    pardf3 = data.frame( par = c( 0, 0, AC$RPEmissionDiff ) + 
                           rep( unlist( bgParam[3:4] ), 2 ), fixed = NA )
    pardf2 = rbind( pardf2, pardf3 )[c(1,2,5,6,3,4,7,8),]
  }
  pardf = rbind( pardf1, pardf2 )

  # edit 
  
  transitionForm = as.formula( " ~ 1" )
  
  if ( dimension == 1 ){
    conrq = rbind( c( rep(0, nrow(pardf1)), -1, 0, 1, 0 ),
                   c( rep(0, nrow(pardf1)),  0,-1, 0, 1 ) )
  
    conr = conrq[ c( MeanFixed[1], SdFixed[1] ), , drop=F ]
    conr.value = AC$APEmissionDiff[c(MeanFixed[1], SdFixed[1])] 

    mod <- depmix( response = active_power ~ 1,  
                   transition = transitionForm, data = tempData, nstates = 2 )
  }else{
    conrq = rbind(  c( rep(0, nrow(pardf1)), -1, 0, 0, 0, 1, 0, 0, 0 ),
                    c( rep(0, nrow(pardf1)),  0, 0, -1, 0, 0, 0, 1, 0 ),
                    c( rep(0, nrow(pardf1)),  0, -1, 0, 0, 0, 1, 0, 0 ),
                    c( rep(0, nrow(pardf1)),  0, 0, 0, -1, 0, 0, 0, 1 )  )

    conr = conrq[ c( MeanFixed, SdFixed ), ,drop=F ]
    conr.value = c( AC$APEmissionDiff[1], AC$RPEmissionDiff[1], 
                    AC$APEmissionDiff[2], AC$RPEmissionDiff[2]
                                        )[ c( MeanFixed, SdFixed) ] 
    mod <- depmix( response = list( active_power~1, reactive_power~1), 
                   family = list( gaussian(), gaussian() ), 
                   transition = transitionForm, data = tempData, nstates = 2 )
  }
  mod = setpars( mod, pardf$par )
  
  if ( all( MeanFixed == c(F,F) ) ){
    fm = depmixS4::fit( mod ) 
  }else fm = depmixS4::fit( mod, conrows = conr, 
                      conrows.upper = conr.value, conrows.lower = conr.value ) 
  
  list( fm = fm, data = tempData, mod = mod, computed = Sys.time() )
}

# list
conrowsF3 <- function( numappl, sd = TRUE ){
  #statesVec = c( 2, 2, 2 )
  # numappl =2 ; sd= F
  statesVec = rep( 2, numappl )
  temp = conrowFRaw( statesVec )
  
  mat.mean = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( x, 0 ) ) ) 
  mat.sd   = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( 0, x ) ) )
  
  mat1 = c(); mat3 = c()
  for( i in 1:nrow(mat.mean)){
    mat1 = rbind( mat1, rbind( mat.mean[i,], mat.sd[i,] ) )
    mat3 = rbind( mat3, rbind( c( 1-sum(mat.mean[i,]), 0 ), c( 0, 1-sum( mat.mean[i, ]))))
  }
  mat2 =  diag(-1, nrow(mat1))
  mat.new = cbind( mat1, mat2, mat3 )
  
  mat.dummy = matrix( 0, nrow = nrow(mat.new), ncol = ncol(mat.new)/2 + (ncol(mat.new)/2)^2 )
  out = cbind( mat.dummy, mat.new  )
  if ( sd == F ) out = out[seq(1, nrow(out), by=2),, drop=F]
  
  out
}

# active_power and reactive_power

conrowsF4 <- function( numappl, sd = T, dim = 2 ){
  # statesVec = c( 2, 2, 2 )
  statesVec = rep( 2, numappl )
  
  nstate = 2^length(statesVec)
  temp = conrowFRaw( statesVec )
  
  mat.a.mean = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( x, 0 ))) 
  mat.a.sd   = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( 0, x )))
  
  if ( dim == 1 ){  
    mat.a.mean = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( x, 0 ) ) ) 
    mat.a.sd   = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( 0, x ) ) )
    
    mat1 = c(); mat3 = c()
    for( i in 1:nrow(mat.a.mean) ){
      mat1 = rbind( mat1, mat.a.mean[i,], mat.a.sd[i,] )
      mat3 = rbind( mat3, c( 1-sum( mat.a.mean[i,] ), 0 ), c( 0, 1-sum( mat.a.sd[i,] ) ) )
    }
    
  }else{
    mat.a.mean = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( x, 0, 0, 0 ) )) 
    mat.a.sd   = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( 0, x, 0, 0 ) ))
    mat.r.mean = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( 0, 0, x, 0 ) )) 
    mat.r.sd   = do.call( cbind, lapply( data.frame( temp ), function(x) cbind( 0, 0, 0, x ) ))
    
    mat1 = c()
    mat3 = c()
    for( i in 1:nrow(mat.a.mean) ){
      mat1 = rbind( mat1, rbind( mat.a.mean[i,], mat.a.sd[i,], 
                                 mat.r.mean[i,], mat.r.sd[i,]  ) )
      mat3 = rbind( mat3, c( 1 - sum( mat.a.mean[i,] ), 0, 0, 0 ), 
                    c( 0, 1 - sum(mat.a.sd[i,] ), 0, 0 ),
                    c( 0, 0, 1 - sum(mat.r.mean[i,]), 0 ),
                    c( 0, 0, 0, 1 - sum(mat.r.sd[i,] ) ))
    }
  }
  
  mat2 =  diag( -1, nrow(mat1) )
  mat.new = cbind( mat1, mat2, mat3 )
  
  mat.dummy = matrix( 0, nrow = nrow(mat.new), ncol = nstate + (nstate)^2 )
  out = cbind( mat.dummy, mat.new )


  if ( sd == F ) out = out[seq(1, nrow(out), by=2),]
  
  out
}

conrowFRaw <- function( stateVec ){
  # stateVec = c( 3, 2, 2, 2 )
  nstates = prod( stateVec )
  
  myExpandGrid = expand.grid( lapply( stateVec , function(x) 0:(x-1) ) )
  
  temp = t(do.call( rbind, lapply( myExpandGrid, function(x) ( sapply( x, function(y) as.integer( intToBits(y) )[ 1:max(x) ] )))))
  subset( temp, rowSums( temp ) > 1 )
}

excludeShort <- function( statevec, min.length = 120 ){
  # hmmfit = temp; transition=F; dimension=2; min.diff= 350
  myrle = rle( statevec )
  mylength = do.call( "c", lapply( myrle$lengths, function(x)  rep( x, each = x ) ))
  statevec[ mylength < min.length ] <- 1
  statevec
}

myStateF <- function( vec, nstateVarg ){
  # vec = fm@posterior$state; nstateVarg = nstateV
  innerf <- function( x ){
    nstateVargm1 = nstateVarg - 1
    totalNstate = prod( nstateVarg )
    cumsum  = cumsum( nstateVarg )
    cumsum2 = cumsum( nstateVargm1 ) 
    
    out = rep( 0, length( nstateVarg ))
    
    if ( x <= sum( nstateVargm1 ) ){
      out[as.numeric( cut( x, breaks = c( 0, cumsum2 ) ))] <- 1
      
    }else if ( x < totalNstate ){
      #temp = conrowF(nstateVarg)[x - sum( nstateVargm1 ), 1:(sum(nstateVargm1))]
      
      temp = conrowFRaw(nstateVarg)[ x - sum(nstateVargm1) , ]
      
      for( i in 1:length(nstateVarg) ){ # i = 1
        temp2 = temp[ (cumsum2[i]-nstateVargm1[i]+1):cumsum2[i] ]
        if ( any( temp2==1) )        out[i] <- which( temp2 ==1 )
      }
    }
    out
  }
  ldply( vec, innerf )
}
