# find.split.pt <- function(power){
#   if(length(power) < 2) return(-10^6)
#   split.interval <- sort(power)[which.max(diff(sort(power))) + c(0,1)]
#   if( diff(split.interval) < width.threshold ) return(-10^6)
#   split.pt       <- mean(split.interval)
#   return( split.pt )
# }
# find.split.pt2 <- function(power, threshold = 50){
#   if(length(power) < 2) return(Inf)
#   split.interval <- data.frame( left = sort(power)[ which(diff(sort(power)) > threshold)     ], 
#                                 rite = sort(power)[ which(diff(sort(power)) > threshold) + 1 ] )
#   split.pt <- ( split.interval$left + split.interval$rite ) / 2
#   return( c( -Inf, sort(split.pt), Inf ) )
# }
find.branch.pt <- function( power, threshold = 50, method = 2 ){
  if( length(power) < 2 ) return( c(-Inf,Inf) )
  power.sort  <- sort(power)
  
  if( method == 1 ){
    pt.ind <- which.max(diff(power.sort))
    branch.pt.left <- power.sort[ pt.ind     ]
    branch.pt.rite <- power.sort[ pt.ind + 1 ]
    if( branch.pt.left + threshold > branch.pt.rite ) return( c(-Inf,Inf) )
  }else if( method == 2 ){
    if( all(diff(power.sort) <= threshold) ) return( c( -Inf, Inf ) )
    pt.ind <- which(diff(power.sort) > threshold)
    branch.pt.left <- power.sort[ pt.ind     ]
    branch.pt.rite <- power.sort[ pt.ind + 1 ]
  }
  branch.pt.cent <- ( branch.pt.left + branch.pt.rite ) / 2

  return( c( -Inf, sort(branch.pt.cent), Inf ) )
}

find.cluster <- function(power, threshold = 50, method = 2) 
  findInterval( power, find.branch.pt(power, threshold, method), rightmost.closed = T )

search.heavy.load.meta <- function(df, bin.threshold = 5, width.threshold = 50, debug.mode = F){
  
  max.iter <- 10
  df <- mutate( df, cluster = 1, ap.re = 1, rp.re = 1, ap.fe = 1, rp.fe = 1 )
  df.list <- list(df)
  
  for( iter in 1:max.iter ){
    df.list <- lapply( df.list, function(df){
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, ap.re = find.cluster(ap.rising,  100)) 
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, rp.re = find.cluster(rp.rising,   50)) 
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, ap.fe = find.cluster(ap.falling, 100)) 
      df <- ddply(df, .(ap.re,rp.re,ap.fe,rp.fe), mutate, rp.fe = find.cluster(ap.falling,  50)) 
      df.sublist <- dlply( df, .(ap.re,rp.re,ap.fe,rp.fe) ) 
    })
    #df.list[ which(lapply( df.list, nrow ) <= 1) ] <- NULL
    df.list <- unlist( df.list, recursive = F )
    df.list <- lapply( seq_along(df.list), function(i){df.list[[i]]$cluster <- i; return(df.list[[i]])} )
    #for( i in 1:length(df.list) ) df.list[[i]]$cluster <- i
  }
  
  # 빈도수가 일정값 이상인 경우만 meta info 만들기
  df.subset <- rbind.fill( df.list[ sapply( df.list, nrow ) > bin.threshold | 
                                      sapply( df.list, function(x) max(difftime(x$end,x$str,units='hour'))) > .5 ] )
  if( is.null(df.subset) ) return(NULL)
  
  df.melt <- melt( df.subset, id.vars='cluster', measure.vars=c('ap.rising', 'rp.rising', 
                                                                'ap.falling', 'rp.falling'))
  cluster.info <- cast( df.melt, cluster ~ variable, c(mean,sd) )
  cluster.info <- merge( cluster.info, ddply( df.melt, .(cluster), summarize, n.pt = length(variable)/4) )
  
  if( all(is.na(cluster.info$ap.rising_sd)) ) return(NULL)
  cluster.info <- dlply( subset(cluster.info, !is.na(ap.rising_sd)), .(cluster))
  cluster.info <- lapply( cluster.info, function(x){
    list( r.edge = c( 'ap_height' = x$ap.rising_mean,  'ap_sigma' = x$ap.rising_sd,
                      'rp_height' = x$rp.rising_mean,  'rp_sigma' = x$rp.rising_sd ),
          f.edge = c( 'ap_height' = x$ap.falling_mean, 'ap_sigma' = x$ap.falling_sd,
                      'rp_height' = x$rp.falling_mean, 'rp_sigma' = x$rp.falling_sd),
          n.pt   = x$n.pt)})
    
  
  if(debug.mode){
    plot( df$ap.rising, df$rp.rising )
    for( i in seq_along(cluster.info) ){
      r.edge <- cluster.info[[i]]$r.edge
      f.edge <- cluster.info[[i]]$f.edge
      x1 = qnorm( 0.01, r.edge['ap_height'], r.edge['ap_sigma'] )
      x2 = qnorm( 0.99, r.edge['ap_height'], r.edge['ap_sigma'] )
      y1 = qnorm( 0.01, r.edge['rp_height'], r.edge['rp_sigma'] )
      y2 = qnorm( 0.99, r.edge['rp_height'], r.edge['rp_sigma'] )
      rect( x1, y1, x2, y2, border = "blue")
    }
    points( df$ap.rising[  df$str %in% df.subset$str ], df$rp.rising[  df$str %in% df.subset$str ], 
            col='red')
  }  
  return(cluster.info)
}

  #   max.iter <- 10
  #   for( iter in 1:max.iter ){
  #     split.pt.new <- unlist(daply( df, .(cluster), summarize, find.split.pt(ap.rising, 50) ))
  #     split.pt.new <- split.pt.new[ split.pt.new > -10^6 ]
  #     if( length(split.pt.new) == 0 ) break
  #     split.pt <- sort(c( split.pt, split.pt.new ))
  #     print(split.pt.new)
  #     df$cluster <- findInterval( df$ap.rising, vec = c( min(df$ap.rising), split.pt, max(df$ap.rising) ), 
  #                                 rightmost.closed=T)
  #   }
  #   abline( v = split.pt, col='red' )
  #   ap.split.pt <- c(0,split.pt,10000)
  #   
  #   df$ap.cluster <- df$cluster
  #   df$cluster <- 1
  #   for( cluster in unique(df$ap.cluster) ){
  #     split.pt <- numeric(0)
  #     df.subset <- df[ df$ap.cluster == cluster, ]
  #     for( iter in 1:max.iter ){
  #       split.pt.new <- unlist(daply( df.subset, .(cluster), summarize, find.split.pt(rp.rising, 50) ))
  #       split.pt.new <- split.pt.new[ split.pt.new > -10^6 ]
  #       if( length(split.pt.new) == 0 ) break
  #       split.pt <- sort(c( split.pt, split.pt.new ))
  #       df.subset$cluster <- findInterval( df.subset$rp.rising, 
  #                                          vec = c( min(df.subset$rp.rising), split.pt, 
  #                                                   max(df.subset$rp.rising) ), 
  #                                          rightmost.closed=T)
  #     }
  #     df$rp.cluster[ df$ap.cluster == cluster ] <- df.subset$cluster
  #     if( length(split.pt) > 0)
  #     segments(ap.split.pt[cluster], split.pt, ap.split.pt[cluster+1], y1 = split.pt, col='blue' )
  #     print(split.pt)
  #   }
  #   
  #   cluster.info <- ddply( df, .(ap.cluster,rp.cluster), summarize,
  #                          ap.mu    = mean(ap.rising), rp.mu    = mean(rp.rising),
  #                          ap.sigma = sd(ap.rising),   rp.sigma = sd(rp.rising), n.pt = length(ap.rising) )
  #   cluster.info <- cluster.info[!is.na(cluster.info$ap.sigma),]

search.meta.info.old <- function(a, condition = T, debug.mode = T){
  
  require(ggplot2)
  if( debug.mode  ){ fig <- ggplot(a, aes(x=active,y=reactive)) + geom_point() }
  
  # step 1
  break.pts      <- seq( min(a$active[a$active>0]), max(a$active[a$active>0]), length.out = 30 )
  active.cut     <- cut( a$active[a$active>50], break.pts, include.lowest = T )
  prior.interval <- rev(break.pts)[ which.max( rev(table( active.cut )) ) + c(1,0) ]
  prior.interval <- search.interval( prior.interval, break.pts, a )
  
  which.on  <- which( findInterval( a$active, prior.interval ) == 1 ) 
  series <- a$active[ which.on ]
  r.edge <- c('ap_height' = mean(series), 'ap_sigma' = sd(series))
  if( any(is.na(r.edge)) ) r.edge[is.na(r.edge)] <- 1
  
  if( debug.mode  ){ fig <- fig + annotate( xmin=prior.interval[1], xmax=prior.interval[2], 
                                            ymin=-Inf, ymax=Inf, fill='blue', alpha=.05,'rect') }
  
  # step 2
  rp.sign <- sign(a$reactive[ which.on ])
  rp.sign <- as.numeric(names(which.max(table(rp.sign))))
  if( rp.sign != 0 ){
    
    break.pts  <- seq( min(a$reactive[a$reactive *rp.sign >0]), 
                       max(a$reactive[a$reactive *rp.sign >0]), length.out = 30 )
    if( max(break.pts) - min(break.pts) < 10 ){
      series <- a$reactive[a$reactive *rp.sign >0]
    }else{
      reactive.cut <- cut( a$reactive[which.on], break.pts, include.lowest = T )
      prior.interval <- break.pts[ which.max( table( reactive.cut ) ) + c(0,1) ]
      
      series <- a$reactive[ which( findInterval( a$reactive, prior.interval ) == 1 ) ]
    }
    r.edge <- c( r.edge, c('rp_height' = mean(series), 'rp_sigma' = sd(series)))
  }else{ r.edge <- c( r.edge, c('rp_height' = 0, 'rp_sigma' = 0)) }      
  
  if( debug.mode  ){ fig <- fig + annotate( ymin=prior.interval[1], ymax=prior.interval[2], 
                                            xmin=-Inf, xmax=Inf, fill='blue', alpha=.05,'rect') }
  
  if( condition ){
    # step 3
    break.pts  <- seq( min(a$active[a$active<0]), max(a$active[a$active<0]), length.out = 30 )
    active.cut <- cut( a$active[a$active< -50], break.pts, include.lowest = T )
    prior.interval <- outer( which( table( active.cut ) %in%  tail( sort(table(active.cut)), 3) ), 
                             c(0,1), function(i,j) break.pts[i+j])
    prior.interval <- prior.interval[which.min( abs(apply( prior.interval, 1, mean ) + r.edge['ap_height'] )),]
    prior.interval <- search.interval( prior.interval, break.pts, a )
    
    which.off  <- which( findInterval( a$active, prior.interval ) == 1 ) 
    series <- a$active[ which.off ]
    f.edge <- c('ap_height' = mean(series), 'ap_sigma' = sd(series))
    if( any(is.na(f.edge)) ) f.edge[is.na(f.edge)] <- 1
    
    if( debug.mode  ){ fig <- fig + annotate( xmin=prior.interval[1], xmax=prior.interval[2], 
                                              ymin=-Inf, ymax=Inf, fill='blue', alpha=.05,'rect') }
    
    # step 4
    if( rp.sign != 0 ){
      
      break.pts  <- seq( min(a$reactive[a$reactive *rp.sign <0]), 
                         max(a$reactive[a$reactive *rp.sign <0]), length.out = 30 )
      
      if( max(break.pts) - min(break.pts) < 10 ){
        series <- a$reactive[a$reactive *rp.sign <0]
      }else{
        reactive.cut <- cut( a$reactive[which.off], break.pts, include.lowest = T )
        prior.interval <- break.pts[ which.max( table( reactive.cut ) ) + c(0,1) ]
        
        series <- a$reactive[ which( findInterval( a$reactive, prior.interval ) == 1 ) ]
      }
      f.edge <- c( f.edge, c('rp_height' = mean(series), 'rp_sigma' = sd(series)))
      if( is.na(f.edge['rp_sigma']) ) f.edge['rp_sigma'] <- 1
      if( debug.mode  ){ fig <- fig + annotate( ymin=prior.interval[1], ymax=prior.interval[2], 
                                                xmin=-Inf, xmax=Inf, fill='blue', alpha=.05,'rect'); print(fig) }
      
    }else{ f.edge <- c( f.edge, c('rp_height' = 0, 'rp_sigma' = 0)) }
  }
  
  # step 5
  which.on  <- which.on[ which(which.on < max(which.off)) ]
  duration.tmp <- sapply( which.on, function(x){min(which.off[x<which.off])} ) - which.on 
  
  if( length(duration.tmp) == 0 ){
    series <- NA
  }else if( length(duration.tmp) == 1 ){
    series <- break.pts
  }else{
    break.pts  <- seq( min(duration.tmp), max(duration.tmp), length.out = 10 )
    active.cut <- cut( duration.tmp, break.pts, include.lowest = T )
    prior.interval <- break.pts[ which.max( table( active.cut ) ) + c(0,1) ]
    
    series <- duration.tmp[ which( findInterval( duration.tmp, prior.interval ) == 1 ) ]
  }
  
  on.off <- c('working_time' = mean(series), 'working_time_sigma' = sd(series))
  if( is.na(on.off['working_time_sigma']) ) on.off['working_time_sigma'] <- 10
  return( list( r.edge, f.edge, on.off ))
}

  
search.cyclic.box.meta <- function(df, bin.threshold = 5, width.threshold = 50, 
                                   debug.mode = F, condition=T, level = 1, 
                                   use.periodicity = T, use.minThre = T){
  library(reshape)
  if(debug.mode) plot( df$active, df$reactive )
    
  if( level == 1 ){ 
    df.original <- df
    df[ which(rowSums(df == 0)==1), c('active','reactive') ] <- 0
  } 
  if( use.minThre ) df[ (35 > abs(df$active)) | (abs(df$active) > 600), c('active','reactive')] <- 0
  
  ap.cluster <- function( ap ){ threshold <- ifelse(length(ap) <= bin.threshold, width.threshold, 
                                                    pmax( ifelse(any(ap==0),30,10), 
                                                          max(diff(sort(diff(sort(ap)))))))
                                find.cluster( ap,threshold=threshold )} 
  rp.cluster <- function( rp ){ threshold <- ifelse(length(rp) <= bin.threshold, width.threshold, 
                                                    pmax( ifelse(any(rp==0),20,10), 
                                                         max(diff(sort(diff(sort(rp)))))))
                                find.cluster( rp,threshold=threshold )} 
  
  max.iter <- 10
  df     <- mutate( df, cluster = 1, ap.re = 1, rp.re = 1, ap.fe = 1, rp.fe = 1 )
  df$idx <- seq_along(df$active)
  
  df.list <- list(df[ rowSums(df[,c('active','reactive')]==0) < 2, ])

  for( iter in 1:max.iter ){
    df.list <- lapply( df.list, function(x){
      #if( all(x$active != 0))
        x <- ddply(x, .(ap.re, rp.re), function(x.sub) {x.sub$ap.re <- ap.cluster(x.sub$active);   return(x.sub)})
      #if( all(x$reactive != 0))
        x <- ddply(x, .(ap.re, rp.re), function(x.sub) {x.sub$rp.re <- rp.cluster(x.sub$reactive); return(x.sub)})
      x.sublist <- dlply( x, .(ap.re,rp.re) ) 
    })
    df.list <- unlist( df.list, recursive = F )
    df.list[ sapply( df.list, nrow ) <= bin.threshold ] <- NULL
    if( length(df.list) == 0 ){
      if( level > 0 ){
        result <- search.cyclic.box.meta(df.original, bin.threshold = bin.threshold, width.threshold = width.threshold, 
                                         debug.mode = debug.mode, level = level - 1)
        return(result)
      }
    }
    df.list <- lapply( seq_along(df.list), function(i){df.list[[i]]$cluster <- i; return(df.list[[i]])} )
  }
  
  # 빈도수가 일정값 이상인 경우만 meta info 만들기
  df.subset <- rbind.fill( df.list[ sapply( df.list, nrow ) > bin.threshold ] )
  names(df.subset)[1] <- 'day'
  check.periodicity <- dlply( df.subset, .(cluster), function(df.cluster) 
    unlist(dlply(df.cluster, .(day), function(x) diff(x$idx)), recursive=F, use.names=F))
  check.periodicity[sapply( check.periodicity, length ) == 0] <- 0
  calc.periodicity <- function(x){ quantile(x,.5) / quantile(x,.05) } 
  need.more.division <- sapply( check.periodicity, calc.periodicity ) > 2
  need.more.division[ is.na(need.more.division) ] <- F
  
  if( any(need.more.division) ){
    for( i in names(check.periodicity)[ need.more.division ] ){
      local.df <- subset( df.subset, cluster == i )
      ap.diff <- diff(sort(local.df$active))
      local.df$division <- local.df$active > sum(sort(local.df$active)[ which.max(ap.diff) + c(0,1) ]) / 2
      tmp <- daply( local.df, .(division), function(x) calc.periodicity(diff(x$idx)) )
      if( max(tmp, na.rm=T) < calc.periodicity( diff(local.df$idx) )){
        df.subset$cluster[ df.subset$idx %in% local.df$idx[local.df$division] ] <- max(df.subset$cluster) + 1
      }
    }
    df.list <- dlply( df.subset, .(cluster) )
  }
  
  
  if( is.null(df.subset) ){
    if( level == 1 ){
      result <- search.cyclic.box.meta(df.original, bin.threshold = bin.threshold, width.threshold = width.threshold, 
                             debug.mode = debug.mode, level = level - 1)
      return(result)
    }else return(NULL) } 
  if( length(which(df.subset$active > 10))==0 | length(unique(df.subset$cluster)) == 1 ){
    if( level == 1 ){
      result <- search.cyclic.box.meta(df.original, bin.threshold = bin.threshold, width.threshold = width.threshold, 
                                       debug.mode = debug.mode, level = level - 1)
      return(result)
    }else return(NULL)
  }
  
  df.melt <- melt( df.subset, id.vars='cluster', measure.vars=c('active', 'reactive'))
  cluster.info <- cast( df.melt, cluster ~ variable, c(mean,sd) )
  cluster.info <- merge( cluster.info, ddply( df.melt, .(cluster), summarize, n.pt = length(variable)/2) )
  
  if( all(is.na(cluster.info$active_sd)) ){
    if( level == 1 ){
      result <- search.cyclic.box.meta(df, bin.threshold = bin.threshold, width.threshold = width.threshold, 
                                       debug.mode = debug.mode, level = level - 1)
      return(result)
    }else return(NULL)
  }
  cluster.info <- subset(cluster.info, !is.na(active_sd) & abs(active_mean) > 35)
  if( length(unique(sign(cluster.info$active_mean))) < 2 ){
    if( level == 1 ){
      result <- search.cyclic.box.meta(df.original, bin.threshold = bin.threshold, width.threshold = width.threshold, 
                                       debug.mode = debug.mode, level = level - 1)
      return(result)
    }else return(NULL)
  }
  
  
  rising  <- cluster.info[ cluster.info$active_mean > 0, ][ which.max(cluster.info$n.pt[ cluster.info$active_mean > 0 ]) ,] 
  
  if( use.periodicity ){
    event   <- df.list[[ rising$cluster ]]$idx
    t.event <- total.event(event)
    if( any(!t.event$is.detected) ){
      candidate <- lapply( which(!t.event$is.detected),
                           function(i) subset( df, ifelse( t.event$is.detected[i-1], 
                                                           t.event$eventTime[i-1], 
                                                           round(t.event$eventTime[(i-1):i]) )< idx & 
                                                 idx < ifelse( t.event$is.detected[i+1], 
                                                               t.event$eventTime[i+1], 
                                                               round(t.event$eventTime[i:(i+1)])) &
                                                 (active != 0 | reactive != 0) ))
      candidate <- ldply( candidate, function(x) x[which.min(abs(x$active-rising$active_mean)),] ) 
      buffer    <- sapply( c(.01,.99), function(x) qnorm(x,rising$active_mean,rising$active_sd) ) + c(-10,10)
      candidate <- candidate[ which( buffer[1] < candidate$active & candidate$active < buffer[2]),]
      names(candidate)[1] <- 'day'
      df.subset <- rbind(  df.subset, candidate )
      if( names(df.list[[ rising$cluster ]])[1] != 'day' ){
        names(df.list[[ rising$cluster ]])[1] <- 'day'
      }
      df.list[[ rising$cluster ]] <- rbind( df.list[[ rising$cluster ]], candidate )
      rising$active_mean   <- mean(df.list[[rising$cluster]]$active)
      rising$active_sd     <- sd(df.list[[rising$cluster]]$active)
      rising$reactive_mean <- mean(df.list[[rising$cluster]]$reactive)
      rising$reactive_sd   <- sd(df.list[[rising$cluster]]$reactive)
      rising$active_min    <- min(df.list[[rising$cluster]]$active)
      rising$active_max    <- max(df.list[[rising$cluster]]$active)
      rising$reactive_min  <- min(df.list[[rising$cluster]]$reactive)
      rising$reactive_max  <- max(df.list[[rising$cluster]]$reactive)
    }
  }
  
  edge.diff <- ddply( cluster.info[ cluster.info$active_mean < 0, ], .(cluster), 
                      function(x) (x$active_mean + rising$active_mean)^2 + (x$reactive_mean + rising$reactive_mean)^2)
  falling <- cluster.info[ cluster.info$cluster == edge.diff$cluster[which.min(edge.diff$V1)], ]
  
  if( use.periodicity ){
    event   <- df.list[[ falling$cluster ]]$idx
    t.event <- total.event(event)
    if( any(!t.event$is.detected) ){
      candidate <- lapply( which(!t.event$is.detected),
                           function(i) subset( df,  ifelse( t.event$is.detected[i-1], 
                                                            t.event$eventTime[i-1], 
                                                            round(t.event$eventTime[(i-1):i]) )< idx & 
                                                 idx < ifelse( t.event$is.detected[i+1], 
                                                               t.event$eventTime[i+1], 
                                                               round(t.event$eventTime[i:(i+1)])) &
                                                 (active != 0 | reactive != 0) ))
      candidate <- ldply( candidate, function(x) x[which.min(abs(x$active-falling$active_mean)),] ) 
      buffer    <- sapply( c(.01,.99), function(x) qnorm(x,falling$active_mean,falling$active_sd) ) + c(-10,10)
      candidate <- candidate[ which( buffer[1] < candidate$active & candidate$active < buffer[2]),]
      buffer    <- sapply( c(.01,.99), function(x) qnorm(x,falling$reactive_mean,falling$reactive_sd) ) + c(-10,10)
      candidate <- candidate[ which( buffer[1] < candidate$reactive & candidate$reactive < buffer[2]),]
      names(candidate)[1] <- 'day'
      df.subset <- rbind(  df.subset, candidate )
      if( '.id' %in% names(df.list[[ falling$cluster ]]) ) 
        names(df.list[[ falling$cluster ]])[names(df.list[[ falling$cluster ]])=='.id'] <- 'day'
      df.list[[ falling$cluster ]] <- rbind( df.list[[ falling$cluster ]], candidate )
      falling$active_mean   <- mean(df.list[[falling$cluster]]$active)
      falling$active_sd     <- sd(df.list[[falling$cluster]]$active)
      falling$reactive_mean <- mean(df.list[[falling$cluster]]$reactive)
      falling$reactive_sd   <- sd(df.list[[falling$cluster]]$reactive)
      falling$active_min    <- min(df.list[[falling$cluster]]$active)
      falling$active_max    <- max(df.list[[falling$cluster]]$active)
      falling$reactive_min  <- min(df.list[[falling$cluster]]$reactive)
      falling$reactive_max  <- max(df.list[[falling$cluster]]$reactive)
    }
  }

  if( nrow(rising) * nrow(falling) == 0 ) return(NULL)
    
  if( sign(rising$active_mean) == sign(falling$active_mean) ) falling <- rising
  if(debug.mode){
    x1 = qnorm( 0.01, rising$active_mean,   rising$active_sd   )
    x2 = qnorm( 0.99, rising$active_mean,   rising$active_sd   )
    y1 = qnorm( 0.01, rising$reactive_mean, rising$reactive_sd )
    y2 = qnorm( 0.99, rising$reactive_mean, rising$reactive_sd )
    rect( x1, y1, x2, y2, border = "blue")
    points( df.subset$active[  df.subset$cluster == rising$cluster ], 
            df.subset$reactive[  df.subset$cluster == rising$cluster ], col='red' )
    
    x1 = qnorm( 0.01, falling$active_mean,   falling$active_sd   )
    x2 = qnorm( 0.99, falling$active_mean,   falling$active_sd   )
    y1 = qnorm( 0.01, falling$reactive_mean, falling$reactive_sd )
    y2 = qnorm( 0.99, falling$reactive_mean, falling$reactive_sd )
    rect( x1, y1, x2, y2, border = "blue")
    points( df.subset$active[  df.subset$cluster == falling$cluster ], 
            df.subset$reactive[  df.subset$cluster == falling$cluster ], col='red' )
    if( sign(x1) != sign(x2) ) falling <- rising
  } 
  
  which.on  <- which( (df$active %in% df.subset$active[ df.subset$cluster == rising$cluster ]) & 
                        (df$reactive %in% df.subset$reactive[ df.subset$cluster == rising$cluster ]) )
  which.off <- which( (df$active %in% df.subset$active[ df.subset$cluster == falling$cluster ]) & 
                        (df$reactive %in% df.subset$reactive[ df.subset$cluster == falling$cluster ]) )
  #duration <- sapply( which.on[ which.on < max(which.off)], 
  #                    function(x) min(which.off[ which.off > x ]) - x )
  duration <- sapply( seq_along((which.on[ which.on < max(which.off)])), 
          function(i) ifelse( any( (which.off > which.on[i]) & (which.on[i+1] > which.off)),
                              min( which.off[ which.off > which.on[i] ] ) - which.on[i], NA ))
  duration <- duration[!duration %in% boxplot.stats(duration)$out] # outlier 제거
  duration <- duration[ !is.na(duration) ]
  
  r.edge <- c( 'ap_height' = rising$active_mean,    'ap_sigma' = rising$active_sd,
               'ap_min'    = rising$active_min,     'ap_max'   = rising$active_max,
               'rp_height' = rising$reactive_mean,  'rp_sigma' = rising$reactive_sd,
               'rp_min'    = rising$reactive_min,   'rp_max'   = rising$reactive_max )
  f.edge <- c( 'ap_height' = falling$active_mean,   'ap_sigma' = falling$active_sd,
               'ap_min'    = falling$active_min,    'ap_max'   = falling$active_max,
               'rp_height' = falling$reactive_mean, 'rp_sigma' = falling$reactive_sd,
               'rp_min'    = falling$reactive_min,  'rp_max'   = falling$reactive_max)
  on.off <- c('working_time' = mean(duration), 'working_time_sigma' = sd(duration))
  if( is.na(on.off['working_time_sigma']) ) on.off['working_time_sigma'] <- 10
  if( abs(r.edge['rp_height']) < 10 | abs(f.edge['rp_height']) < 10){
    r.edge[c('rp_height','rp_sigma')] <- 0
    f.edge[c('rp_height','rp_sigma')] <- 0
  }
  if(!is.finite(on.off['duty_cycle'])) on.off['duty_cycle'] <- 0
  return( list( r.edge, f.edge, on.off ))
}

remove.Peacky <- function(df){
  df$active_power   <- envelopeDetector(df, use.active.power=T)$LowerEnvelope
  df$reactive_power <- envelopeDetector(df, use.active.power=F)$UpperEnvelope
  return(df)  
}

search.interval <- function( prior.interval, break.pts, edge ){
  
  converged = T
  while( converged ){
    if( max(break.pts) == prior.interval[2] ) break
    prior.interval.old <- prior.interval
    prior.interval[2] <- min(break.pts[break.pts > max(prior.interval.old)])
    if( all(prior.interval == prior.interval.old) ) break
    converged <- !all(findInterval( edge$active, prior.interval ) == findInterval( edge$active, prior.interval.old ))
  }
  
  converged = T
  while( converged ){
    if( min(break.pts) == prior.interval[1] ) break
    prior.interval.old <- prior.interval
    prior.interval[1] <- max(break.pts[break.pts < min(prior.interval.old)])
    if( all(prior.interval == prior.interval.old) ) break
    converged <- !all(findInterval( edge$active, prior.interval ) == findInterval( edge$active, prior.interval.old ))
  }
  return( prior.interval )
}

