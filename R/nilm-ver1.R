meta.to.box.Shape.first <- function(data,meta.list, Debug.mode = F){
  
  data.list <- list(data)
  box.list  <- list()
  for( i in 1:length(meta.list) ){
    temp <- find.box.shape( meta.list[[i]], data.list[[i]] )
    data.list[[length(data.list)+1]] <- data.frame( timestamp = data$timestamp, 
                                                    active_power = temp$ap.box.reverse, 
                                                    reactive_power = temp$rp.box.reverse )
    box.list[[length(box.list)+1]] <- data.frame( timestamp = data$timestamp, 
                                                  p = temp$ap.box, 
                                                  q = temp$rp.box )
    if( Debug.mode ) plot( temp$ap.box, type='l'); title(main = "active power"  )
    if( any(temp$rp.box != 0) ){plot( temp$rp.box, type='l'); title(main = "reactive power")}  
  }
  if( Debug.mode ) plot(temp$ap.box.reverse, type='l'); title('AP : residual')
  if( Debug.mode ) plot(temp$rp.box.reverse, type='l'); title('RP : residual')  
  return(box.list)
}

generate.meta.info <- function( data, max.iter = 3, debug.mode = F, split.algorithm = 1){
  generate.CyclicBox.meta( data=data, max.iter=max.iter, debug.mode=debug.mode, split.algorithm=split.algorithm)}

generate.StandbyPower.meta <- function( data ){
  val <- min(envelopeDetector(data)$UpperEnvelope)
  return(list(from.para.to.JSON3( data=data, StandbyPowerValue=val )))
}
  
generate.CyclicBox.meta <- function( data, max.iter = 3, debug.mode = F, split.algorithm = 1, apply.envelopeDetector=F, bin.threshold = 5, use.edgeMatching = T ){
  
  require(lubridate)
  require(plyr)
  
  if( split.algorithm == 1 ){
    n.hour <- 4
    # 4시간 연속 사용량이 가장 낮은 구간 찾기
    tmp <- mutate(data, date = floor_date(timestamp,'hour') , day = floor_date(timestamp,'day'))
    tmp <- ddply( tmp, .(day), transform, sufficient = length(unique(date)) > 4 )
    tmp <- subset( tmp, sufficient )
    hourly.sum <- ddply( tmp, .(day, date), summarize, usage = mean(active_power))
    hourly.sum <- ddply( hourly.sum, .(day), summarize,
                         str.t = date[which.min(filter(usage, rep(1,4), sides=2))-1] )
    hourly.sum <- dlply(hourly.sum, .(day), function(x) interval(x$str.t, x$str.t + dhours(4)))
    data.list  <- lapply( hourly.sum, function(p) subset( data, timestamp %within% p ))  
  }else if( split.algorithm == 2 ){
    data.list <- split.data1( data=data,  debug.mode=debug.mode, apply.envelopeDetector = apply.envelopeDetector)
  }else if( split.algorithm == 3 ){
    data.list <- list(data)
  }
  
  result     <- list()
  para.lists <- list()
  json.lists <- list()
  active.list   <- lapply( data.list, function(dt) list(dt$active_power))
  reactive.list <- lapply( data.list, function(dt) list(dt$reactive_power))
  
  basal.threshold <- .2 
  for( iter in 1:max.iter ){
    
    preprocessed.data.list <- 
      mapply( function( df, ap.list, rp.list ){
      data.frame( timestamp=df$timestamp, 
                  active_power   = getLast(ap.list), 
                  reactive_power = getLast(rp.list) )}, 
      data.list, active.list, reactive.list, SIMPLIFY = F)
    
    
    if( split.algorithm == 1 ){
      remove.data <- sapply( preprocessed.data.list, function(x) 
        length(which(x$active_power > quantile( data$active_power, .05) + 30)) / 
          length(x$active_power) < .1 )  
    }else{
      remove.data <- sapply( preprocessed.data.list, function(x) 
        all( x$active_power < quantile( data$active_power, .05) + 30))
    }
    if( any(remove.data) ){
      preprocessed.data.list[ which(remove.data) ] <- NULL
      data.list[ which(remove.data) ]     <- NULL
      active.list[ which(remove.data) ]   <- NULL
      reactive.list[ which(remove.data) ] <- NULL
    } 
    if( length(preprocessed.data.list) == 0 ) return(json.lists)
    
    # edge detect
    edge <- ldply( preprocessed.data.list, edge.detect )
    
    condition = T
    if( iter > 1 ) condition = sum(rbind.fill(lapply( result, function(df) getLast(df)))$active_power) != 0
    
    #para <- search.meta.info.old(edge, condition, debug.mode = debug.mode)
    if( length(para.lists) > 0 ) para.old <- para
    para <- search.cyclic.box.meta(edge, debug.mode = debug.mode, condition=condition )
    if( is.null(para) ) break
    if( length(para.lists) > 0 && identical(para.old,para) ) para <- search.cyclic.box.meta(edge, debug.mode = debug.mode, condition=condition, level = 0 )
    if( sign(para[[1]]['ap_height']) == sign(para[[2]]['ap_height']) ){
      para <- search.cyclic.box.meta(edge, debug.mode = debug.mode, width.threshold=10 )
      if( sign(para[[1]]['ap_height']) == sign(para[[2]]['ap_height']) ) break
    }
      
#     if( ( sign(qnorm(.01,para[[1]]['ap_height'],para[[1]]['ap_sigma'])) !=  
#             sign(qnorm(.99,para[[1]]['ap_height'],para[[1]]['ap_sigma'])) ) |
#           ( sign(qnorm(.01,para[[2]]['ap_height'],para[[2]]['ap_sigma'])) !=  
#               sign(qnorm(.99,para[[2]]['ap_height'],para[[2]]['ap_sigma'])) ) |
#           ( sign(qnorm(.01,para[[1]]['rp_height'],para[[1]]['rp_sigma'])) !=  
#               sign(qnorm(.99,para[[1]]['rp_height'],para[[1]]['rp_sigma'])) ) |
#           ( sign(qnorm(.01,para[[2]]['rp_height'],para[[2]]['rp_sigma'])) !=  
#               sign(qnorm(.99,para[[2]]['rp_height'],para[[2]]['rp_sigma'])) ) )
#       para <- search.cyclic.box.meta(edge, debug.mode = debug.mode, width.threshold=10 )
#     if( is.null(para) ) break
#     if( para[[3]]['working_time'] < 30 ){
#       preprocessed.data.list <- lapply( preprocessed.data.list, remove.Peacky)
#       edge <- ldply( preprocessed.data.list, edge.detect )
#       para <- search.cyclic.box.meta(edge, debug.mode = debug.mode )
#     }
    r.edge <- para[[1]]
    f.edge <- para[[2]]
    on.off <- para[[3]]    
    
print( r.edge )
print( f.edge )
print( on.off )
    meta <- from.para.to.JSON( data, r.edge, f.edge, on.off, 1 )
    temp.list <- lapply( preprocessed.data.list, function(df) find.box.shape(meta,df) )
    active.list <- mapply( function(ap.list,box.removed) 
      append( ap.list, list(box.removed$ap.box.reverse)), active.list, temp.list, SIMPLIFY = F )
    reactive.list <- mapply( function(rp.list,box.removed) 
      append( rp.list, list(box.removed$rp.box.reverse)), reactive.list, temp.list, SIMPLIFY = F )
    result <- append( result, list( mapply( function(box.removed,df) 
      data.frame( timestamp = df$timestamp, active_power = box.removed$ap.box, reactive_power = box.removed$rp.box ),
      temp.list, data.list, SIMPLIFY = F )))
    
    temp <- rbind.fill( temp.list )
    if( max(abs(temp$ap.box)) + max(abs(temp$rp.box)) != 0 ){
      
      #plot( temp$ap.box, type='l'); title(main = "active power")
      #plot( temp$rp.box, type='l'); title(main = "reactive power")
      
      # update parameters
      #box.result   <- series.to.box.lists( data$timestamp, temp$ap.box, temp$rp.box )
      #rising.edge  <- ddply( box.result[[1]], .(str), function(x){find.jump.value(data,x$str)} )
      #falling.edge <- ddply( box.result[[1]], .(end), function(x){find.jump.value(data,x$end)} )
      #duration.on  <- box.result[[1]]$duration
      #duration.off <- box.result[[1]]$str[-1] - box.result[[1]]$end[-nrow(box.result[[1]])]
      #duration.off <- duration.off[!duration.off %in% boxplot.stats(duration.off)$out] # outlier 제거
      box.result.list <- mapply( function(df,box.list) 
        series.to.box.lists( df$timestamp, box.list$ap.box, box.list$ap.box )[[1]],
        data.list, temp.list, SIMPLIFY = F )
      box.result   <- rbind.fill( box.result.list )
      box.result   <- box.result[ box.result$max.p > 10, ]
      duration.on  <- box.result$duration
      duration.on  <- duration.on[!duration.on %in% boxplot.stats(duration.on)$out] # outlier 제거
      duration.off <- box.result$str[-1] - box.result$end[-nrow(box.result)]
      duration.off <- duration.off[ duration.off > 2]
      duration.off <- duration.off[!duration.off %in% boxplot.stats(duration.off)$out] # outlier 제거
      
      #       r.edge <- c('ap_height' = mean(rising.edge$ap), 'ap_sigma' = sd(rising.edge$ap),
      #                   'rp_height' = mean(rising.edge$rp), 'rp_sigma' = sd(rising.edge$rp))
      #       f.edge <- c('ap_height' = mean(falling.edge$ap), 'ap_sigma' = sd(falling.edge$ap),
      #                   'rp_height' = mean(falling.edge$rp), 'rp_sigma' = sd(falling.edge$rp))
      on.off <- c('working_time' = mean(duration.on), 
                  'working_time_sigma' = sd(duration.on),
                  'duty_cycle' = mean(duration.on) / (mean(duration.on)+mean(duration.off)))
      if( !is.finite(on.off['duty_cycle']) | is.na(on.off['duty_cycle']) ) on.off['duty_cycle'] <- 0
      if( is.na(on.off['working_time_sigma']) ) on.off['working_time_sigma'] <- 10
      para.lists[[length(para.lists)+1]] <- list( 'rising_edge' = r.edge, 'falling_edge' = f.edge, 'cycle' = on.off)
      json.lists[[length(json.lists)+1]] <- from.para.to.JSON( data, r.edge, f.edge, on.off, nrow(box.result) )
    }
  }
  
  #preprocessed.data <- data.frame( timestamp      = data$timestamp, 
  #                                 active_power   = active.list[[length(active.list)]],
  #                                 reactive_power = reactive.list[[length(reactive.list)]] )
  #edgeMatching.result <- extract.M.info.using.edgeMatching(preprocessed.data)
if( use.edgeMatching ){
  preprocessed.data.list <- mapply( function(df, ap.list, rp.list){
    data.frame( timestamp=df$timestamp, active_power=getLast(ap.list), reactive_power=getLast(rp.list))}, 
    data.list, active.list, reactive.list, SIMPLIFY = F)
  matching.fn <- function(x) extract.M.info.using.edgeMatching(x, debug.mode = debug.mode)
  edgeMatching.result <- lapply( preprocessed.data.list, matching.fn )
  edgeMatching.result <- unlist( edgeMatching.result, recursive = F )
  if( is.null(edgeMatching.result) ) return(json.lists)
  
  #json.lists <- merge.metaInfo( json.lists, edgeMatching.result )
  json.lists <- merge.metaInfo.new( json.lists, edgeMatching.result )
  
}
  
  json.lists <- json.lists[laply( json.lists, function(x) abs(x$rising_edge + x$falling_edge)[1] < 100 )]
  json.lists <- json.lists[laply( json.lists, function(x) abs(x$rising_edge + x$falling_edge)[3] < 100 )]
  json.lists <- json.lists[laply( json.lists, function(x) !is.na(x$cycle[1]) && x$cycle[1] > 0 )]
  #json.lists <- json.lists[laply( json.lists, function(x) x$rising[1] < 500 )]
  json.lists <- json.lists[laply( json.lists, function(x) x$rising[1] > 35 )]
  json.lists <- json.lists[laply( json.lists, function(x) x$box.no >= bin.threshold )]
  names(json.lists) <- seq_along(json.lists)
  return(json.lists)
}

getLast <- function( list.data ) list.data[[length(list.data)]]




extract.M.info.using.edgeMatching <- function( data, debug.mode = F ){
  
  library(clue)
  library(ggplot2)
  
  THRESHOLD_RISING  <-  30
  THRESHOLD_FALLING <- -30
  
  data$ts <- data$timestamp
  data$ap <- data$active_power
  data$rp <- data$reactive_power
  
  data$ap <- remove.consecutive.jump.pt2(data$ap)
  data$rp <- remove.consecutive.jump.pt2(data$rp)
  edgeList <- edgeDetector( data, DEBUG_MODE=F) #little bit slow
  if( any( sapply( edgeList, nrow ) == 0 ) ) return(NULL)
  risingEdges  <- edgeList[[1]]
  fallingEdges <- edgeList[[2]]
  
  # Name columns
  # This can be moved to edgeDetector 
  colnames(risingEdges) <- c( 'idx', paste0('risingSteadyStateHeight', c('Ap','Rp')), 
                              outer( c('rising','falling','flat'), 
                                     c('Length','HeightAp','HeightRp'), paste0 ))
  
  colnames(fallingEdges) <- c( 'idx', paste0('fallingSteadyStateHeight', c('Ap','Rp')), 
                               outer( c('falling','rising','flat'), 
                                      c('Length','HeightAp','HeightRp'), paste0 ))
  
  
  # Subset edges having SteadyStateHeightAp > Threshold
  sr <- subset( risingEdges,  risingSteadyStateHeightAp  >= THRESHOLD_RISING  )
  sf <- subset( fallingEdges, fallingSteadyStateHeightAp <= THRESHOLD_FALLING )
  
  ap <- data$ap; base <- min(ap)
  solution <- edgeMatching(sr, sf, ap, base) 
  if( nrow(solution[[1]]) == 0 ) return(NULL)
  mr <- solution[[1]] # matched rising edges
  mf <- solution[[2]] # matched falling edges
  mr <- mr[!is.na(mf$idx),]
  mf <- mf[!is.na(mf$idx),]
  data$idx <- 1:nrow(data)
  if( nrow(mr) == 0 ) return(NULL)
  if( debug.mode ) print(plotMatching(1, max(data$idx), data, mr, mf))
  
  meta <- extract.meta.info( data, mr$idx, mf$idx)
  
  meta <- meta[,names(meta) %in% c('active.ON','active.OFF','reactive.ON','reactive.OFF','duration.ON')]
  Fun    <- function(i,j) norm( meta[i,] - meta[j,],type='2'); VecFun <- Vectorize( Fun )
  tmp <- outer( 1:nrow(meta), 1:nrow(meta), VecFun )
  tmp[ lower.tri(tmp) ] <- 0
  
  similarity.threshold <- 100
  ind <- which((0 < tmp) & (tmp < similarity.threshold), arr.ind=T); ind <- as.data.frame(ind)
  ind.cluster <- dlply( ind, .(row), function(x) c(x[1,1],x[,2]) )
  
  if( length(ind.cluster) > 1 ){
    for( i in 1:(length(ind.cluster)-1) ){
      for( j in (i+1):length(ind.cluster) ){
        if( length(intersect( ind.cluster[[i]], ind.cluster[[j]])) > 0 ){
          ind.cluster[[i]] <- union( ind.cluster[[i]], ind.cluster[[j]] )
          ind.cluster[[j]] <- NA
        }
      }
    }
    ind.cluster <- ind.cluster[!is.na(ind.cluster)]
  }else{
    ind.cluster <- as.list(1:nrow(meta))
  }
  
  cluster.mu    <- ldply( ind.cluster, function(x) apply( meta[x,], 2, mean) )
  cluster.sigma <- ldply( ind.cluster, function(x) apply( meta[x,], 2, sd  ) )
  cluster.mu$duration.OFF <- 0    
  cluster.sigma[is.na(cluster.sigma)] <- 1
  json.meta.lists <- list()
  for( i in 1:nrow(cluster.mu) ){
    r.edge <- c( 'ap_height' = cluster.mu$active.ON[i], 'ap_sigma' = cluster.sigma$active.ON[i], 
                 'rp_height' = cluster.mu$reactive.ON[i], 'rp_sigma' = cluster.sigma$reactive.ON[i] )
    f.edge <- c( 'ap_height' = cluster.mu$active.OFF[i], 'ap_sigma' = cluster.sigma$active.OFF[i], 
                 'rp_height' = cluster.mu$reactive.OFF[i], 'rp_sigma' = cluster.sigma$reactive.OFF[i] )
    on.off <- c('working_time' = cluster.mu$duration.ON[i], 'working_time_sigma' = cluster.sigma$duration.ON[i],
                'duty_cycle' = 0 )
    
    json.meta.lists[[length(json.meta.lists)+1]] <- from.para.to.JSON( data, r.edge, f.edge, on.off, 1 )
  }
  json.meta.lists <- json.meta.lists[laply( json.meta.lists, function(x) x$cycle[1] > 0 )]
  
  return(json.meta.lists)
}

extract.meta.info <- function( data, str.ind, end.ind ){
  ind <- data.frame( str = str.ind, end = end.ind )
  ind$duration.ON <- ind$end - ind$str + 1
  ind <- cbind(ind,ldply( str.ind, function(x){ find.jump.value(data,x)}))
  ind <- cbind(ind,ldply( end.ind, function(x){ find.jump.value(data,x)}))
  names(ind)[4:7] <- c('active.ON','reactive.ON','active.OFF','reactive.OFF')
  return(ind)
}

find.jump.value <- function( data, ind ){
  
  if( is.vector(data) ){
    data        <- remove.consecutive.jump.pt2(data)
    first.state <- find.stationary.state( data[ c(1:ind)], forward=F)[[2]]
    secon.state <- find.stationary.state( data[-c(1:ind)], forward=T)[[2]]
    jump.value  <- secon.state - first.state
    return(c('power' = jump.value))
  }
  if( !('ap' %in% names(data)) ) data$ap <- data$active_power
  if( !('rp' %in% names(data)) ) data$rp <- data$reactive_power
  data$ap <- remove.consecutive.jump.pt2(data$ap)
  data$rp <- remove.consecutive.jump.pt2(data$rp)
  
  first.state.ap <- find.stationary.state( data$ap[ c(1:ind)], forward=F)[[2]]
  secon.state.ap <- find.stationary.state( data$ap[-c(1:ind)], forward=T)[[2]]
  first.state.rp <- find.stationary.state( data$rp[ c(1:ind)], forward=F)[[2]]
  secon.state.rp <- find.stationary.state( data$rp[-c(1:ind)], forward=T)[[2]]
  jump.value.ap <- secon.state.ap - first.state.ap
  jump.value.rp <- secon.state.rp - first.state.rp
  return(c('ap' = jump.value.ap, 'rp' = jump.value.rp))
}

edge.detect <- function(data, threshold.active = 10, threshold.reactive = 10, remove.peak = F){
  
  if( remove.peak ){
    tmp  <- envelopeDetector( data, use.active.power=T ); data$active_power <- tmp$LowerEnvelope
    tmp  <- envelopeDetector( data, use.active.power=F ); data$active_power <- tmp$UpperEnvelope
  } 
  edge <- data.frame( active   = diff(remove.consecutive.jump.pt2(data$active_power)), 
                      reactive = diff(remove.consecutive.jump.pt2(data$reactive_power)) )
  edge$active[ abs(edge$active) < threshold.active] <- 0
  edge$reactive[ abs(edge$reactive) < threshold.reactive] <- 0
  return(edge)  
}


find.heavy.load <- function( meta.json, data, debug.mode ){
  
  box.reverse <- find.box.shape.reverse5( data=data, json.para=meta.json, threshold=10, level=2, recursive=T )
  box.shape   <- data$active_power - box.reverse[[1]]
  result <- data.frame( 'timestamp' = data$timestamp, 
                        'ap.box' = box.shape, 'rp.box' = 0, 
                        'ap.box.reverse' = box.reverse[[1]], 
                        'rp.box.reverse' = data$reactive_power )
  
  if( debug.mode ){
    
    data   <- merge( data[,c('timestamp','active_power','reactive_power')], result, by='timestamp' )
    
    r.edge <- box.reverse[[2]]
    f.edge <- box.reverse[[3]]
    r.edge <- data[r.edge,c('timestamp','active_power')] 
    f.edge <- data[f.edge,c('timestamp','active_power')] 
    
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
    #points(x=r.edge$timestamp, y=r.edge$active_power, col='red',  pch=17)
    #points(x=f.edge$timestamp, y=f.edge$active_power, col='blue', pch=25)
    abline( v=r.edge$timestamp, col='red' ) 
    abline( v=f.edge$timestamp, col='blue')
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
    
    par(mfrow=c(3,2), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power,   ann=FALSE, xaxt="n", type='l')
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l', yaxt='n'); axis(side=4)
    plot(x=data$timestamp, y=data$ap.box, ann=FALSE, xaxt="n", type='l')
    plot(x=data$timestamp, y=data$rp.box, ann=FALSE, xaxt="n", type='l', yaxt='n'); axis(side=4)
    plot(x=data$timestamp, y=data$ap.box.reverse, ann=FALSE, xaxt="n", type='l')
    plot(x=data$timestamp, y=data$rp.box.reverse, ann=FALSE, xaxt="n", type='l', yaxt='n'); axis(side=4)
    title(paste('Step2 : Box removing 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Active / Reactive power", 1, 1, outer=TRUE)
    mtext("Residual / Box / Original",     2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
    
    #     draw.fig <- function(d){
    #       df <- subset( data.melt, day == d & variable %in% c('active_power','reactive_power') )
    #       re <- subset( r.edge,    day == d )
    #       fe <- subset( f.edge,    day == d )
    #       
    #       fig1 <- ggplot( df, aes(x=timestamp, y=value)) + geom_line(aes(colour=variable, group=variable)) +
    #         geom_point( data=re, aes(y=active_power), colour='red',   bg='red',   shape=24, size=3 ) +
    #         geom_point( data=fe, aes(y=active_power), colour='black', bg='white', shape=25, size=3 ) +
    #         theme(legend.position='bottom') +
    #         labs(title= paste('Step1 : Edge detection 결과\n', min(df$timestamp), '--', max(df$timestamp),'\n'))
    #       
    #       df <- subset( data.melt, day == d )
    #       
    #       fig2 <- ggplot( df, aes(x=timestamp, y=value)) + 
    #         geom_line(aes(colour=variable, group=variable)) + 
    #         facet_grid(variable~., scales='free') +
    #         theme(legend.position='NONE') +
    #         labs(title= paste('Step2 : Box removing 결과\n', min(df$timestamp), '--', max(df$timestamp),'\n'))
    #       
    #       return( list( fig1, fig2 ))      
    #     }
    #     print( lapply( unique(data.melt$day), draw.fig ) )
    
  }
  result$p <- result$ap.box
  result$q <- result$rp.box
  return( result )  
}

find.box.shape <- function( meta.json, data, show.fig = F ){
  
  if( abs(meta.json$rising_edge['rp_height']) < 35 ) meta.json$rising_edge['rp_height'] <- 0
  if( !is.na(meta.json$rising_edge['rp_height']) & !is.na(meta.json$falling_edge['rp_height']) & 
        (meta.json$rising_edge['rp_height'] != 0) & (meta.json$falling_edge['rp_height'] != 0) ){
    
    box.reverse.ap <- find.box.shape.reverse4( data$active_power, data$reactive_power, 10
                                               , meta.json$rising_edge[c('rp_height','rp_sigma')]
                                               , meta.json$falling_edge[c('rp_height','rp_sigma')]
                                               , meta.json$cycle[c('working_time','working_time_sigma')] )
    
    if( is.null(box.reverse.ap) | !is.list(box.reverse.ap) ){
      result <- data.frame( 'timestamp' = data$timestamp, 
                            'ap.box' = 0, 'rp.box' = 0, 
                            'ap.box.reverse' = data$active_power, 
                            'rp.box.reverse' = data$reactive_power )
      return(result)
    }
    box.reverse.rp <- find.box.shape.reverse4( data$reactive_power, data$reactive_power, 10
                                               , meta.json$rising_edge[c('rp_height','rp_sigma')]
                                               , meta.json$falling_edge[c('rp_height','rp_sigma')]
                                               , meta.json$cycle[c('working_time','working_time_sigma')] )
    if( is.null(box.reverse.rp) | !is.list(box.reverse.rp)  ){
      result <- data.frame( 'timestamp' = data$timestamp, 
                            'ap.box' = 0, 'rp.box' = 0, 
                            'ap.box.reverse' = data$active_power, 
                            'rp.box.reverse' = data$reactive_power )
      return(result)
    }
    
    box.reverse.rp <- box.reverse.rp[['residual']] 
    
  }else{
    
    meta.json$rising_edge[c('rp_height','rp_sigma')] <- meta.json$rising_edge[c('ap_height','ap_sigma')] 
    meta.json$falling_edge[c('rp_height','rp_sigma')] <- meta.json$falling_edge[c('ap_height','ap_sigma')] 
    box.reverse.ap <- find.box.shape.reverse4( data$active_power, data$reactive_power, 10 
                                               , meta.json$rising_edge[c('rp_height','rp_sigma')]
                                               , meta.json$falling_edge[c('rp_height','rp_sigma')]
                                               , meta.json$cycle[c('working_time','working_time_sigma')]
                                               , is.rp.zero = T  )
    if( is.null(box.reverse.ap) | !is.list(box.reverse.ap) ){
      result <- data.frame( 'timestamp' = data$timestamp, 
                            'ap.box' = 0, 'rp.box' = 0, 
                            'ap.box.reverse' = data$active_power, 
                            'rp.box.reverse' = data$reactive_power )
      return(result)
    }
    box.reverse.rp <- data$reactive_power
  }
  
  box.ap <- data$active_power   - box.reverse.ap[['residual']]
  box.rp <- data$reactive_power - box.reverse.rp
  
  result <- data.frame( 'timestamp' = data$timestamp, 
                        'ap.box' = box.ap, 'rp.box' = box.rp, 
                        'ap.box.reverse' = box.reverse.ap[['residual']], 
                        'rp.box.reverse' = box.reverse.rp )
  
  under.peak <- which( apply( cbind( shift( result$ap.box != 0, -1 ), 
                                     result$ap.box == 0, 
                                     shift( result$ap.box != 0, 1 )), 1, all ))
  if( length(under.peak) > 0 ){
    under.peak <- under.peak[(result$ap.box[ under.peak - 1 ] - result$ap.box[ under.peak + 1 ]) < 10]
    if( length(under.peak) > 0 ){
      peak.value <- (result$ap.box[under.peak-1]+result$ap.box[under.peak+1])/2
    result$ap.box[ under.peak ] <- result$ap.box[ under.peak ] + peak.value
    result$ap.box.reverse[ under.peak ] <- result$ap.box.reverse[ under.peak ] - peak.value}
  }
 
  
  if( show.fig ){
    
    data <- merge( data[,c('timestamp','active_power','reactive_power')], result, by='timestamp' )
    data.melt <- reshape::melt( data, id.vars='timestamp' )

    r.edge <- data[box.reverse.ap[['edge']]$re,c('timestamp','active_power')] 
    f.edge <- data[box.reverse.ap[['edge']]$fe,c('timestamp','active_power')] 
    
    data.melt$day <- day(data.melt$timestamp)
    r.edge$day    <- day(r.edge$timestamp)
    f.edge$day    <- day(f.edge$timestamp)
    
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power, ann=FALSE,   xaxt="n", type='l')
    #points(x=r.edge$timestamp, y=r.edge$active_power, col='red',  pch=17)
    #points(x=f.edge$timestamp, y=f.edge$active_power, col='blue', pch=25)
    abline( v=r.edge$timestamp, col='red' ) 
    abline( v=f.edge$timestamp, col='blue')
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))
        
    par(mfrow=c(3,2), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(x=data$timestamp, y=data$active_power,   ann=FALSE, xaxt="n", type='l')
    plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l', yaxt='n'); axis(side=4)
    plot(x=data$timestamp, y=data$ap.box, ann=FALSE, xaxt="n", type='l')
    plot(x=data$timestamp, y=data$rp.box, ann=FALSE, xaxt="n", type='l', yaxt='n'); axis(side=4)
    plot(x=data$timestamp, y=data$ap.box.reverse, ann=FALSE, xaxt="n", type='l')
    plot(x=data$timestamp, y=data$rp.box.reverse, ann=FALSE, xaxt="n", type='l', yaxt='n'); axis(side=4)
    title(paste('Step2 : Box removing 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Active / Reactive power", 1, 1, outer=TRUE)
    mtext("Residual / Box / Original",     2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))

#     draw.fig <- function(d){
#       df <- subset( data.melt, day == d & variable %in% c('active_power','reactive_power') )
#       re <- subset( r.edge,    day == d )
#       fe <- subset( f.edge,    day == d )
#       
#       fig1 <- ggplot( df, aes(x=timestamp, y=value)) + geom_line(aes(colour=variable, group=variable)) +
#         geom_point( data=re, aes(y=active_power), colour='red',   bg='red',   shape=24, size=3 ) +
#         geom_point( data=fe, aes(y=active_power), colour='black', bg='white', shape=25, size=3 ) +
#         theme(legend.position='bottom') +
#         labs(title= paste('Step1 : Edge detection 결과\n', min(df$timestamp), '--', max(df$timestamp),'\n'))
#       
#       df <- subset( data.melt, day == d )
#       
#       fig2 <- ggplot( df, aes(x=timestamp, y=value)) + 
#         geom_line(aes(colour=variable, group=variable)) + 
#         facet_grid(variable~., scales='free') +
#         theme(legend.position='NONE') +
#         labs(title= paste('Step2 : Box removing 결과\n', min(df$timestamp), '--', max(df$timestamp),'\n'))
#       
#       return( list( fig1, fig2 ))      
#     }
#     print( lapply( unique(data.melt$day), draw.fig ) )
    
  }
  result$p <- result$ap.box
  result$q <- result$rp.box
  return(result)
}

merge.metaInfo <- function( meta.main, meta.sub, similarity.threshold = 100 ){
  
  Fun    <- function(i,j){
    para1 <- c( meta[[i]]$rising_edge[c(1,3)],  meta[[i]]$falling_edge[c(1,3)], meta[[i]]$cycle[1] )
    para2 <- c( meta[[j]]$rising_edge[c(1,3)],  meta[[j]]$falling_edge[c(1,3)], meta[[j]]$cycle[1] )
    return( norm( para1 - para2, type = '2'))
  }
  VecFun <- Vectorize( Fun )
  
  add.two.mean  <- function( mu1, n1, mu2, n2 ) (mu1*n1 + mu2*n2) / (n1+n2)
  add.two.sigma <- function( mu1, sigma1, n1, mu2, sigma2, n2 ){
    sqrt( ((sigma1^2 + mu1^2)*n1 + (sigma2^2 + mu2^2)*n2) / (n1+n2) - add.two.mean( mu1, n1, mu2, n2 )^2 )
  } 
  
  merge.two.meta <- function( meta1, meta2 ){
    meta  <- meta1
    n1 <- meta1$box.no
    n2 <- meta2$box.no
    meta$rising_edge['ap_height'] <- 
      add.two.mean( meta1$rising_edge['ap_height'], n1, meta2$rising_edge['ap_height'], n2)
    meta$rising_edge['rp_height'] <- 
      add.two.mean( meta1$rising_edge['rp_height'], n1, meta2$rising_edge['rp_height'], n2)
    meta$falling_edge['ap_height'] <- 
      add.two.mean( meta1$falling_edge['ap_height'], n1, meta2$falling_edge['ap_height'], n2)
    meta$falling_edge['rp_height'] <- 
      add.two.mean( meta1$falling_edge['rp_height'], n1, meta2$falling_edge['rp_height'], n2)
    meta$cycle['working_time'] <- 
      add.two.mean( meta1$cycle['working_time'], n1, meta2$cycle['working_time'], n2)
    
    meta$rising_edge['ap_sigma'] <- 
      add.two.sigma( meta1$rising_edge['ap_height'], meta1$rising_edge['ap_sigma'], n1,
                     meta2$rising_edge['ap_height'], meta2$rising_edge['ap_sigma'], n2 )
    meta$rising_edge['rp_sigma'] <- 
      add.two.sigma( meta1$rising_edge['rp_height'], meta1$rising_edge['rp_sigma'], n1,
                     meta2$rising_edge['rp_height'], meta2$rising_edge['rp_sigma'], n2 )
    meta$falling_edge['ap_sigma'] <- 
      add.two.sigma( meta1$falling_edge['ap_height'], meta1$falling_edge['ap_sigma'], n1,
                     meta2$falling_edge['ap_height'], meta2$falling_edge['ap_sigma'], n2 )
    meta$falling_edge['rp_sigma'] <- 
      add.two.sigma( meta1$falling_edge['rp_height'], meta1$falling_edge['rp_sigma'], n1,
                     meta2$falling_edge['rp_height'], meta2$falling_edge['rp_sigma'], n2 )
    meta$cycle['working_time_sigma'] <- 
      add.two.sigma( meta1$cycle['working_time'], meta1$cycle['working_time_sigma'], n1,
                     meta2$cycle['working_time'], meta2$cycle['working_time_sigma'], n2 )
    meta$box.no <- n1 + n2
    return(meta)
  }
  
  meta <- c( meta.main, meta.sub )
  tmp <- outer( 1:length(meta), 1:length(meta), VecFun )
  tmp[ lower.tri(tmp) ] <- 0
  
  ind <- which((0 < tmp) & (tmp < similarity.threshold), arr.ind=T); ind <- as.data.frame(ind)
  ind.cluster <- dlply( ind, .(row), function(x) c(x[1,1],x[,2]) )
  if( length(ind.cluster) > 1 ){
    for( ii in 1:2 ){
      for( i in 1:(length(ind.cluster)-1) ){
        for( j in (i+1):length(ind.cluster) ){
          if( length(intersect( ind.cluster[[i]], ind.cluster[[j]])) > 0 ){
            ind.cluster[[i]] <- union( ind.cluster[[i]], ind.cluster[[j]] )
            ind.cluster[[j]] <- NA
          }
        }
      }
    }
    ind.cluster <- ind.cluster[!is.na(ind.cluster)]
  }else{
    ind.cluster <- as.list(1:length(meta))
  }
  
  if( length(ind.cluster) > 0){
    for( i in 1:length(ind.cluster) ){
      j <- min(ind.cluster[[i]])
      for( k in ind.cluster[[i]][ind.cluster[[i]]>j] ){
        meta[[j]]  <- merge.two.meta(meta[[j]], meta[[k]])
        meta[[k]] <- NA
      }      
    }
  } 
  meta <- meta[!is.na(meta)]
  return( meta )
}

merge.metaInfo.new <- function( meta.main, meta.sub, alpha.threshold = .01 ){
  
  Fun <- function(i,j){
    
    para1 <- c( meta[[i]]$rising_edge,  meta[[i]]$falling_edge, meta[[i]]$cycle )
    para2 <- c( meta[[j]]$rising_edge,  meta[[j]]$falling_edge, meta[[j]]$cycle )
    
    left.threshold <- alpha.threshold
    rite.threshold <- 1 - alpha.threshold
    
    ap.left <- qnorm( left.threshold, para1['ap_height'],    para1['ap_sigma'])
    ap.rite <- qnorm( rite.threshold, para1['ap_height'],    para1['ap_sigma'])
    rp.left <- qnorm( left.threshold, para1['rp_height'],    para1['rp_sigma'])
    rp.rite <- qnorm( rite.threshold, para1['rp_height'],    para1['rp_sigma'])
    wt.left <- qnorm( left.threshold, para1['working_time'], para1['working_time_sigma'])
    wt.rite <- qnorm( rite.threshold, para1['working_time'], para1['working_time_sigma'])
    
    ap.included <- ( ap.left <= para2['ap_height'] ) & ( para2['ap_height'] <= ap.rite )
    rp.included <- ( rp.left <= para2['ap_height'] ) & ( para2['ap_height'] <= rp.rite )
    wt.included <- ( wt.left <= para2['working_time'] ) & ( para2['working_time'] <= wt.rite )

    if( (abs(para1['rp_height']) < 10) | (abs(para2['rp_height']) < 10) ){
      is.included <- ap.included & wt.included 
    }else{ 
      is.included <- ap.included & wt.included & rp.included
    }
    
    return( is.included )
  }
  VecFun <- Vectorize( Fun )
    
  add.two.mean  <- function( mu1, n1, mu2, n2 ) (mu1*n1 + mu2*n2) / (n1+n2)
  add.two.sigma <- function( mu1, sigma1, n1, mu2, sigma2, n2 ){
    sqrt( ((sigma1^2 + mu1^2)*n1 + (sigma2^2 + mu2^2)*n2) / (n1+n2) - add.two.mean( mu1, n1, mu2, n2 )^2 )
  } 
  
  merge.two.meta <- function( meta1, meta2 ){
    
    meta  <- meta1
    n1 <- meta1$box.no
    n2 <- meta2$box.no
    
    meta$rising_edge['ap_height'] <- 
      add.two.mean( meta1$rising_edge['ap_height'], n1, meta2$rising_edge['ap_height'], n2)
    meta$rising_edge['rp_height'] <- 
      add.two.mean( meta1$rising_edge['rp_height'], n1, meta2$rising_edge['rp_height'], n2)
    meta$falling_edge['ap_height'] <- 
      add.two.mean( meta1$falling_edge['ap_height'], n1, meta2$falling_edge['ap_height'], n2)
    meta$falling_edge['rp_height'] <- 
      add.two.mean( meta1$falling_edge['rp_height'], n1, meta2$falling_edge['rp_height'], n2)
    meta$cycle['working_time'] <- 
      add.two.mean( meta1$cycle['working_time'], n1, meta2$cycle['working_time'], n2)
    
    meta$rising_edge['ap_sigma'] <- 
      add.two.sigma( meta1$rising_edge['ap_height'], meta1$rising_edge['ap_sigma'], n1,
                     meta2$rising_edge['ap_height'], meta2$rising_edge['ap_sigma'], n2 )
    meta$rising_edge['rp_sigma'] <- 
      add.two.sigma( meta1$rising_edge['rp_height'], meta1$rising_edge['rp_sigma'], n1,
                     meta2$rising_edge['rp_height'], meta2$rising_edge['rp_sigma'], n2 )
    meta$falling_edge['ap_sigma'] <- 
      add.two.sigma( meta1$falling_edge['ap_height'], meta1$falling_edge['ap_sigma'], n1,
                     meta2$falling_edge['ap_height'], meta2$falling_edge['ap_sigma'], n2 )
    meta$falling_edge['rp_sigma'] <- 
      add.two.sigma( meta1$falling_edge['rp_height'], meta1$falling_edge['rp_sigma'], n1,
                     meta2$falling_edge['rp_height'], meta2$falling_edge['rp_sigma'], n2 )
    meta$cycle['working_time_sigma'] <- 
      add.two.sigma( meta1$cycle['working_time'], meta1$cycle['working_time_sigma'], n1,
                     meta2$cycle['working_time'], meta2$cycle['working_time_sigma'], n2 )
    meta$box.no <- n1 + n2
    return(meta)
  }
  
  meta <- c( meta.main, meta.sub )
  tmp <- outer( 1:length(meta), 1:length(meta), VecFun )
  tmp[ lower.tri(tmp) ] <- F
  
  ind <- which(tmp, arr.ind=T); ind <- as.data.frame(ind)
  ind.cluster <- dlply( ind, .(row), function(x) c(x[1,1],x[,2]) )
  if( length(ind.cluster) > 1 ){
    for( ii in 1:2 ){
      for( i in 1:(length(ind.cluster)-1) ){
        for( j in (i+1):length(ind.cluster) ){
          if( length(intersect( ind.cluster[[i]], ind.cluster[[j]])) > 0 ){
            ind.cluster[[i]] <- union( ind.cluster[[i]], ind.cluster[[j]] )
            ind.cluster[[j]] <- NA
          }
        }
      }
    }
    ind.cluster <- ind.cluster[!is.na(ind.cluster)]
  }else{
    ind.cluster <- as.list(1:length(meta))
  }
  
  if( length(ind.cluster) > 0){
    for( i in 1:length(ind.cluster) ){
      j <- min(ind.cluster[[i]])
      for( k in ind.cluster[[i]][ind.cluster[[i]]>j] ){
        meta[[j]]  <- merge.two.meta(meta[[j]], meta[[k]])
        meta[[k]] <- NA
      }      
    }
  } 
  meta <- meta[!is.na(meta)]
  return( meta )
}



