generate.PatternScanHeavy.meta_1Hz_extend <- function ( data
                                                        , relApDiff.max = .15
                                                        , relRpDiff.max = .20
                                                        , min_sigmag.active.n = 800
                                                        , min_sigmag.reactive.n = 3
                                                        , timeSpan = dhours(1)
                                                        , measurementErrorBound = 10
                                                        , debug.mode = F
                                                        , bin.threshold = 5) {


  # search starting points
  s.pattern.original <- DetectPattern_1Hz_extended( data, "start"
                                                    , min_sigmag.active.n
                                                    , min_sigmag.reactive.n
                                                    , debug.mode = F )
  message( "Searching starting points : DONE" )
  
  # search ending points
  e.pattern.original <- DetectPattern_1Hz_extended( data, "end"
                                                    , min_sigmag.active.n * (1-relApDiff.max)
                                                    , min_sigmag.reactive.n
                                                    , 0
                                                    , debug.mode = F )
  message( "Searching ending points : DONE" )
  
  base.value       <- min_sigmag.active.n + min(data$active_power)
  risingEdge.diff  <- diff(s.pattern.original$start.timestamp) > timeSpan
  risingEdge.drop  <- mapply( function(r1,r2) any( data$active_power[r1:r2] < base.value ),
                              head(s.pattern.original$end.idx,-1),
                              tail(s.pattern.original$end.idx,-1) )

  # rising edge 간격이 timeSpan 보다 크거나, 중간에 baseline까지 떨어지면 다른 group으로 간주
  # 동일 group 내의 edge들만 matching 적용
  if( nrow(s.pattern.original) == 1 ){
    risingEdge.group = list(s.pattern.original)
  }else{
    risingEdge.group <- mapply( function(i,j) s.pattern.original[i:j,], SIMPLIFY=F,
                                c(1, which(risingEdge.diff | risingEdge.drop) + 1 ),
                                c(which(risingEdge.diff | risingEdge.drop), nrow(s.pattern.original)) )
  }

  # rising edge의 그룹을 중심으로 falling edge의 그룹을 지정
  tmp <- findInterval( e.pattern.original$start.timestamp,
                       vec = rbind.fill(lapply( risingEdge.group, function(x) x[1,] ))$start.timestamp )
  fallingEdge.group <- lapply( seq_along(risingEdge.group),
                               function(i) e.pattern.original[which(i == tmp),] )

  # 동일 그룹 내의 rising/falling edge에 대해 matching condition 지정
  # Step 1 : 그룹 내에 동일 갯수의 r/f edge가 있고, 순서대로 값이 거의 일치하면 적합하다고 판단
  is.matched <- mapply( function(r,f) ( nrow(r) == nrow(f) &&
                                          all(abs((r$ap.delta + f$ap.delta) / r$ap.delta) < relApDiff.max) &&
                                          all(abs((r$rp.delta + f$rp.delta) / r$rp.delta) < relRpDiff.max)),
                        risingEdge.group, fallingEdge.group )

  rEdge <- rbind.fill( risingEdge.group[ is.matched ]);  risingEdge.group[ is.matched ] <- NULL
  fEdge <- rbind.fill(fallingEdge.group[ is.matched ]); fallingEdge.group[ is.matched ] <- NULL

  # Step 2
  matching.a <- mapply( function(r,f) list( outer( r$ap.delta, f$ap.delta, function(x,y) abs((x+y)/x) ),
                                            outer( r$rp.delta, f$rp.delta, function(x,y) abs((x+y)/x) ),
                                            outer( r$start.timestamp, f$start.timestamp, '-' )),
                        risingEdge.group, fallingEdge.group, SIMPLIFY=F )
  matching.b <- lapply( matching.a,
                        function(x) x[[1]] < relApDiff.max & x[[2]] < relRpDiff.max & x[[3]] < 0 )


  for( i in seq_along(risingEdge.group) ){

    matching.table <- matching.b[[i]]
    rEdge.no       <- nrow(matching.table)
    matching.edge  <- NULL

    if( !is.null(rEdge.no) && ncol(matching.table) > 0 ){
      for( row in 1:rEdge.no ){

        matching.log    <- matching.table[row,]
        rEdge.candidate <- risingEdge.group[[i]][ row, ]
        fEdge.candidate <- fallingEdge.group[[i]][ which(matching.log), ]

        if( nrow(fEdge.candidate) != 0 ){

          for( j in 1:nrow(fEdge.candidate) ){
            period <- interval( rEdge.candidate$end.timestamp,
                                fEdge.candidate$start.timestamp[j])
            if( any(data$active_power[ data$timestamp %within% period ] <
                      abs(rEdge.candidate$ap.delta) / 2) ){
              col <- which(matching.log)[j]
              matching.log[col:length(matching.log)] <- F
              break
            }else if( abs( data$active_power[ rEdge.candidate$start.idx ] -
                             data$active_power[fEdge.candidate$end.idx[j] ] ) < measurementErrorBound ){
              col <- which(matching.log)[j]
              if( col + 1 <= length(matching.log) ) matching.log[(col+1):length(matching.log)] <- F
              break
            }
          }
          fEdge.candidate <- fallingEdge.group[[i]][ which(matching.log), ]

          if( nrow(rEdge.candidate) == 1 & nrow(fEdge.candidate) == 1 ){
            rEdge <- rbind( rEdge, rEdge.candidate )
            fEdge <- rbind( fEdge, fEdge.candidate )
            matching.edge <- rbind( matching.edge, c(row,which(matching.log)))
            matching.table[,which(matching.log)] <- F
          }
        }
      }
      if( !is.null(matching.edge) && nrow(matching.edge) > 0 ){
        risingEdge.group[[i]]  <-  risingEdge.group[[i]][ - matching.edge[,1], ]
        fallingEdge.group[[i]] <- fallingEdge.group[[i]][ - matching.edge[,2], ]
      }
    }
  }

  meta <- search.PatternScanHeavy.meta( data, rEdge, fEdge, bin.threshold = bin.threshold, debug.mode = debug.mode )

  return(list( rEdge, fEdge, meta ))
}

search.PatternScanHeavy.meta <- function(data, rEdge, fEdge, bin.threshold = 5, width.threshold = 50, debug.mode = F){

  if(nrow(rEdge) != nrow(fEdge)) stop("Dimension does not match")

  max.iter <- 10
  names(rEdge) <- paste0('r.',names(rEdge))
  names(fEdge) <- paste0('f.',names(fEdge))
  df <- cbind( rEdge, fEdge )
  df <- mutate( df, cluster = 1, ap.re = 1, rp.re = 1, ap.fe = 1, rp.fe = 1 )
  df.list <- list(df)

  for( iter in 1:max.iter ){
    df.list <- lapply( df.list, function(df){
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, ap.re = find.cluster(r.ap.delta))
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, rp.re = find.cluster(r.rp.delta))
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, ap.fe = find.cluster(f.ap.delta))
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, rp.fe = find.cluster(f.rp.delta))
      df.sublist <- dlply( df, .(ap.re,rp.re,ap.fe,rp.fe) )
    })
    #df.list[ which(lapply( df.list, nrow ) <= 1) ] <- NULL
    df.list <- unlist( df.list, recursive = F )
    df.list <- lapply( seq_along(df.list), function(i){df.list[[i]]$cluster <- i; return(df.list[[i]])} )
    #for( i in 1:length(df.list) ) df.list[[i]]$cluster <- i
  }

  # 빈도수가 일정값 이상인 경우만 meta info 만들기
  df.subset <- rbind.fill( df.list[ sapply( df.list, nrow ) >= bin.threshold ] )
  if( is.null(df.subset) ) return(NULL)

  df.subset$duration <- df.subset$f.end.timestamp - df.subset$r.start.timestamp
  df.melt <- melt( df.subset, id.vars='cluster', measure.vars=c('r.ap.delta', 'r.rp.delta',
                                                                'f.ap.delta', 'f.rp.delta',
                                                                'duration'))
  df.melt$variable <- gsub('.delta','',df.melt$variable)
  cluster.info <- cast( df.melt, cluster ~ variable, c(min,median,max) )
  cluster.info <- merge( cluster.info, ddply( df.melt, .(cluster), summarize, n.pt = length(variable)/5) )

  if(debug.mode){
    plot( df$r.ap.delta, df$r.rp.delta )
    for( i in 1:nrow(cluster.info) ){
      localClusterInfo <- cluster.info[i,]
      rect( localClusterInfo$r.ap_min,
            localClusterInfo$r.rp_min,
            localClusterInfo$r.ap_max,
            localClusterInfo$r.rp_max, border = "blue")
    }
  }

  json.result <- list()
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  for( i in 1:nrow(cluster.info) ){
    localCluster.r <- cluster.info[i, grep('r.',names(cluster.info), fixed=T)]
    localCluster.f <- cluster.info[i, grep('f.',names(cluster.info), fixed=T)]
    localCluster.d <- cluster.info[i, grep('duration',names(cluster.info), fixed=T)]
    localCluster.n <- cluster.info[i, 'n.pt']
    names(localCluster.r) <- gsub('r.','',names(localCluster.r), fixed=T)
    names(localCluster.f) <- gsub('f.','',names(localCluster.f), fixed=T)

    localCluster.r <- setNames( as.numeric(localCluster.r), names(localCluster.r) )
    localCluster.f <- setNames( as.numeric(localCluster.f), names(localCluster.f) )
    localCluster.d <- setNames( as.numeric(localCluster.d), names(localCluster.d) )
    localCluster.n <- as.numeric(localCluster.n)

    type.table <- table(subset( df.subset, cluster == cluster.info[i,'cluster'] )$r.type)
    type       <- names(which.max(type.table))
    type       <- switch( type, 'reactive/active' = 1, '-reactive/active' = -1, 'active' = 2)

    json.result[[length(json.result)+1]] <-
      list(`meta-version` = 1,
           shape_type = "pattern_scan_heavy",
           summit_flag = type,
           rising_edge  = localCluster.r,
           falling_edge = localCluster.f,
           duration     = localCluster.d,
           n.pt         = localCluster.n,
           generation_info = list(data_used = list(start = str.t, end = end.t, sampling = 1),
                                  computed = as.character(Sys.time())))
  }
  return(json.result)
}


meta2PatternScanHeavy_1Hz <- function (data
                                       , meta
                                       , extension.p = .2
                                       , relApDiff.max = .15
                                       , relRpDiff.max = .20
                                       , measurementErrorBound = 10
                                       , debug.mode
                                       , filename = 'meta2patternHeavy'
                                       , ignoreRP = F) {

  if( debug.mode ) png( paste0(filename,"-%d.png"), width=3000)

  # pattern matching
  r.edge    <- setNames( as.numeric(meta$rising_edge),  names(meta$rising_edge)  )
  f.edge    <- setNames( as.numeric(meta$falling_edge), names(meta$falling_edge) )
  duration  <- meta$duration
  answer.log <- data.frame( timestamp = data$timestamp, p= 0, q= 0)

  s.rp_min <- r.edge['rp_min'] * (1-extension.p*sign(r.edge['rp_min']))
  s.rp_max <- r.edge['rp_max'] * (1+extension.p*sign(r.edge['rp_max']))
  s.rp_med <- r.edge['rp_median']

  s.ap_min <- r.edge['ap_min'] * (1-extension.p*sign(r.edge['ap_min']))
  s.ap_max <- r.edge['ap_max'] * (1+extension.p*sign(r.edge['ap_max']))
  s.ap_med <- r.edge['ap_median']

  e.rp_min <- f.edge['rp_min'] * (1-extension.p*sign(f.edge['rp_min']))
  e.rp_max <- f.edge['rp_max'] * (1+extension.p*sign(f.edge['rp_max']))

  e.ap_min <- f.edge['ap_min'] * (1-extension.p*sign(f.edge['ap_min']))
  e.ap_max <- f.edge['ap_max'] * (1+extension.p*sign(f.edge['ap_max']))
  e.ap_med <- f.edge['ap_median']

  main_type <- switch( as.character(meta$summit_flag), '1' = 'reactive', '-1' = '-reactive', '2' = 'active'   )
  sub_type  <- switch( as.character(meta$summit_flag), '1' = 'active',   '-1' = 'active',    '2' = 'reactive' )

  # search start points
  s.pattern <- DetectPattern_1Hz_new( data = data
                                      , position   = "start"
                                      , main_type  = main_type
                                      , sub_type   = sub_type
                                      , debug.mode = F )

  if( as.character(meta$summit_flag) == 2 ){
    s.pattern <- subset( s.pattern, ((delta >= s.ap_min) & (delta <= s.ap_max)) )
    if( !ignoreRP ){
      s.pattern <- subset( s.pattern, (sub.delta >= s.rp_min) & (sub.delta <= s.rp_max) )
    }else{
      s.pattern <- subset( s.pattern, abs(sub.delta) < s.rp_max * 2 )
    }
  }else{
    s.pattern <- subset( s.pattern, (sub.delta >= s.ap_min) & (sub.delta <= s.ap_max))
    if( !ignoreRP ){
      s.pattern <- subset( s.pattern, (delta >= s.rp_min) & (delta <= s.rp_max) )
    }else{
      s.pattern <- subset( s.pattern, abs(delta) < s.rp_max * 2 )
    }
  }
  s.pattern <- s.pattern[order(s.pattern$start.timestamp),]

  # search end points
  e.pattern <- DetectPattern_1Hz_new( data = data
                                      , position   = "end"
                                      , main_type  = main_type
                                      , sub_type   = sub_type
                                      , c_factor   = 0
                                      , debug.mode = F )

  if( as.character(meta$summit_flag) == 2 ){
    e.pattern <- subset( e.pattern, ((delta >= e.ap_min) & (delta <= e.ap_max)) )
    if( !ignoreRP ){
      e.pattern <- subset( e.pattern, (sub.delta >= e.rp_min) & (sub.delta <= e.rp_max) )
    }else{
      e.pattern <- subset( e.pattern, abs(sub.delta) < e.rp_max * 2 )
    }
  }else{
    e.pattern <- subset( e.pattern, (sub.delta >= e.ap_min) & (sub.delta <= e.ap_max))
    if( !ignoreRP ){
      e.pattern <- subset( e.pattern, (delta >= e.rp_min) & (delta <= e.rp_max) )
    }else{
      e.pattern <- subset( e.pattern, abs(delta) < e.rp_max * 2 )
    }
  }
  e.pattern <- e.pattern[ order(e.pattern$start.timestamp),]

  names(s.pattern)[ names(s.pattern) == 'delta' ] <-
    paste0(switch( as.character(meta$summit_flag), '1' = 'rp', '-1' = 'rp', '2' = 'ap'),'.delta')
  names(s.pattern)[ names(s.pattern) == 'sub.delta' ] <-
    paste0(switch( as.character(meta$summit_flag), '1' = 'ap', '-1' = 'ap', '2' = 'rp'),'.delta')
  names(e.pattern)[ names(e.pattern) == 'delta' ] <-
    paste0(switch( as.character(meta$summit_flag), '1' = 'rp', '-1' = 'rp', '2' = 'ap'),'.delta')
  names(e.pattern)[ names(e.pattern) == 'sub.delta' ] <-
    paste0(switch( as.character(meta$summit_flag), '1' = 'ap', '-1' = 'ap', '2' = 'rp'),'.delta')


  if(debug.mode){
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
    if( nrow(s.pattern) > 0 ) abline( v=s.pattern$start.timestamp, col='red'  )
    if( nrow(e.pattern) > 0 ) abline( v=e.pattern$start.timestamp, col='blue' )
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)),
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }

  if( nrow(s.pattern) == 0 || nrow(e.pattern) == 0 ) return(answer.log)

  s.pattern0 <- subset( s.pattern, start.timestamp  > max(e.pattern$end.timestamp) )
  s.pattern  <- subset( s.pattern, start.timestamp <= max(e.pattern$end.timestamp) )

  e.pattern0 <- subset( e.pattern, start.timestamp  < min(s.pattern$end.timestamp) )
  e.pattern  <- subset( e.pattern, start.timestamp >= min(s.pattern$end.timestamp) )

  risingEdge.group  <- split.data.frame( s.pattern, findInterval( s.pattern$start.idx, vec = e.pattern$start.idx ))
  fallingEdge.group <- split.data.frame( e.pattern, findInterval( e.pattern$start.idx, vec = s.pattern$start.idx ))

   if( nrow(s.pattern0) > 0 ){
     risingEdge.group <- append( risingEdge.group, dlply(s.pattern0, .(start.idx )) )
     tmp <-
       rep( list(subset( e.pattern, start.timestamp > min(s.pattern0$start.timestamp) )), nrow(s.pattern0) )
     fallingEdge.group <- append( fallingEdge.group, tmp )
   }

   if( nrow(e.pattern0) > 0 ){
     fallingEdge.group <- append( dlply(e.pattern0, .(start.idx)), fallingEdge.group )
     tmp <-
       rep( list(subset( s.pattern, start.timestamp < max(e.pattern0$start.timestamp) )), nrow(e.pattern0) )
     risingEdge.group <- append( tmp, risingEdge.group )
   }


  # 동일 그룹 내의 rising/falling edge에 대해 matching condition 지정
  # Step 1 : 그룹 내에 동일 갯수의 r/f edge가 있고, 순서대로 값이 거의 일치하면 적합하다고 판단
  is.matched <- mapply( function(r,f) ( nrow(r) == nrow(f) &&
                                          all(abs((r$ap.delta + f$ap.delta) / r$ap.delta) < relApDiff.max) &&
                                          all(abs((r$rp.delta + f$rp.delta) / r$rp.delta) < relRpDiff.max)),
                        risingEdge.group, fallingEdge.group )

  rEdge <- rbind.fill( risingEdge.group[ is.matched ]);  risingEdge.group[ is.matched ] <- NULL
  fEdge <- rbind.fill(fallingEdge.group[ is.matched ]); fallingEdge.group[ is.matched ] <- NULL

  # Step 2 & 3
  for( step in 2:3 ){

    if( step == 2 ){
      matching.a <- mapply( function(r,f) list( outer( r$ap.delta, f$ap.delta, function(x,y) abs((x+y)/x) ),
                                                outer( r$rp.delta, f$rp.delta, function(x,y) abs((x+y)/x) ),
                                                outer( r$start.timestamp, f$start.timestamp, '-' )),
                            risingEdge.group, fallingEdge.group, SIMPLIFY=F )
      matching.b <- lapply( matching.a,
                            function(x) x[[1]] < relApDiff.max & x[[2]] < relRpDiff.max & x[[3]] < 0 )
    }else if( step == 3 ){
      matching.a <- mapply( function(r,f) list( outer( r$ap.delta, f$ap.delta, function(x,y) abs((x+y)/x) ),
                                                outer( r$start.timestamp, f$start.timestamp, '-' )),
                            risingEdge.group, fallingEdge.group, SIMPLIFY=F )
      matching.b <- lapply( matching.a,
                            function(x) x[[1]] < relApDiff.max & x[[2]] < 0 )
    }


    for( i in seq_along(risingEdge.group) ){

      matching.table <- matching.b[[i]]
      rEdge.no       <- nrow(matching.table)
      matching.edge  <- NULL

      if( !is.null(rEdge.no) ){
        for( row in 1:rEdge.no ){

          matching.log    <- matching.table[row,]
          rEdge.candidate <- risingEdge.group[[i]][ row, ]
          fEdge.candidate <- fallingEdge.group[[i]][ which(matching.log), ]

          if( nrow(fEdge.candidate) != 0 ){

            for( j in 1:nrow(fEdge.candidate) ){
              period <- interval( rEdge.candidate$end.timestamp,
                                  fEdge.candidate$start.timestamp[j])
              if( any(data$active_power[ data$timestamp %within% period ] <
                        abs(rEdge.candidate$ap.delta) / 2) ){
                col <- which(matching.log)[j]
                matching.log[col:length(matching.log)] <- F
                break
              }else if( abs( data$active_power[ rEdge.candidate$start.idx ] -
                               data$active_power[fEdge.candidate$end.idx[j] ] ) < measurementErrorBound ){
                col <- which(matching.log)[j]
                if( col + 1 <= length(matching.log) ) matching.log[(col+1):length(matching.log)] <- F
                break
              }
            }
            fEdge.candidate <- fallingEdge.group[[i]][ which(matching.log), ]

            if( nrow(rEdge.candidate) == 1 & nrow(fEdge.candidate) == 1 ){
              rEdge <- rbind( rEdge, rEdge.candidate )
              fEdge <- rbind( fEdge, fEdge.candidate )
              matching.edge <- rbind( matching.edge, c(row,which(matching.log)))
              matching.table[,which(matching.log)] <- F
            }
          }
        }
        if( !is.null(matching.edge) && nrow(matching.edge) > 0 ){
          risingEdge.group[[i]]  <-  risingEdge.group[[i]][ - matching.edge[,1], ]
          fallingEdge.group[[i]] <- fallingEdge.group[[i]][ - matching.edge[,2], ]
        }
      }
    }
  }

  if( length(risingEdge.group) > 0 && nrow(risingEdge.group <- rbind.fill(risingEdge.group)) > 0 ){
    for( i in 1:nrow(risingEdge.group) ){
      halfEdge <- ( data$active_power[ data$timestamp == risingEdge.group$start.timestamp[i] ] +
                      data$active_power[ data$timestamp == risingEdge.group$end.timestamp[i] ] ) / 2
      localData <- subset( data, timestamp %within% interval( risingEdge.group$end.timestamp[i],
                                                              risingEdge.group$end.timestamp[i] + dhours(1) ))
      idx <- min(which( (localData$active_power > halfEdge) & shift((localData$active_power < halfEdge), -1) ))
      fEdge <- rbind.fill( fEdge, data.frame( end.timestamp = localData$timestamp[idx]
                                              , ap.delta = -risingEdge.group$ap.delta[i]
                                              , rp.delta = -risingEdge.group$rp.delta[i]) )
      rEdge <- rbind.fill( rEdge, risingEdge.group[i,])
    }
  }

  if( length(fallingEdge.group) > 0 && nrow(fallingEdge.group <- rbind.fill(fallingEdge.group)) > 0 ){
    for( i in 1:nrow(fallingEdge.group) ){
      halfEdge <- ( data$active_power[ data$timestamp == fallingEdge.group$start.timestamp[i] ] +
                      data$active_power[ data$timestamp == fallingEdge.group$end.timestamp[i] ] ) / 2
      localData <- subset( data, timestamp %within% interval( fallingEdge.group$end.timestamp[i] - dhours(1),
                                                              fallingEdge.group$end.timestamp[i] ))
      idx <- max(which( (localData$active_power < halfEdge) & shift((localData$active_power > halfEdge), -1) ))
      rEdge <- rbind.fill( rEdge, data.frame( start.timestamp = localData$timestamp[idx]
                                              , ap.delta = -fallingEdge.group$ap.delta[i]
                                              , rp.delta = -fallingEdge.group$rp.delta[i]) )
      fEdge <- rbind.fill( fEdge, fallingEdge.group[i,])
    }
  }


  if(debug.mode){
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
    if( nrow(rEdge) > 0 ) abline( v=rEdge$start.timestamp, col='red'  )
    if( nrow(fEdge) > 0 ) abline( v=fEdge$end.timestamp, col='blue' )
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, type='l')
    title(paste('Step2 : Edge Matching 결과\n', min(data$timestamp), '--', max(data$timestamp)),
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
  }

  if (nrow(rEdge) != 0) {

    # log extraction
    p.pad <- ( s.ap_med + median(rEdge$ap.delta) ) / 2
    q.pad <- ( s.rp_med + median(rEdge$rp.delta) ) / 2

    signalOn <- unlist( mapply( function(str,end){which(answer.log$timestamp %within% interval(str, end))},
                                rEdge$start.timestamp, fEdge$end.timestamp ) )
    answer.log$p[signalOn] <- p.pad
    answer.log$q[signalOn] <- q.pad


    if(debug.mode){
      str.t <- min(answer.log$timestamp[ answer.log$p != 0 ]) - eminutes(10)
      end.t <- max(answer.log$timestamp[ answer.log$p != 0 ]) + eminutes(10)
      par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
      plot(x=data$timestamp[ data$timestamp %within% interval(str.t,end.t)],
           y=data$active_power[ data$timestamp %within% interval(str.t,end.t)],
           ann=FALSE,   xaxt="n", type='l')
      plot(x=answer.log$timestamp[ answer.log$timestamp %within% interval(str.t,end.t)],
           y=answer.log$p[ answer.log$timestamp %within% interval(str.t,end.t)],
           ann=FALSE, type='l')
      title(paste('Step3 : consumption approximation\n', str.t, '--', end.t),
            outer=TRUE)
      mtext("Timestamp", 1, 1, outer=TRUE)
      mtext("approx / Active power", 2, 3, outer=TRUE, las=0)
      par(mfrow=c(1,1))
    }
    if( debug.mode ) dev.off()
    return(answer.log)
  } else {
    if( debug.mode ) dev.off()
    stop("Error: no matching point!")
  }
}
