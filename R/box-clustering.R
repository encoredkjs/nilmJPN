find.stationary.state <- function( power, forward = T, threshold = 10 ){
  
  if( length(power) == 1 ) return(list(0,power))
  if( !forward ) power <- rev(power)
  if( length(power) > 4 ){
    for( i in 1:(length(power)-4) ){
      subinterval <- power[i:(i+4)]
      #if( max(abs(diff(subinterval))) < threshold ) return(list(i-1,median(subinterval)))
      if( diff(range(subinterval)) < threshold ) return(list(i-1,median(subinterval)))
    }  
  }
  i <- which.min(abs(diff(power)))
  return(list(i, power[i+1]))
}

remove.consecutive.jump.pt2 <- function( power, threshold = 10, level = 1 ){
  
  if( level == 0 ) return(power)
  consecutive.jump.pt <- which( abs(diff(power)) > threshold )
  consecutive.jump.pt <- union( consecutive.jump.pt, consecutive.jump.pt + 1 )
  
  if( length(consecutive.jump.pt) == 0 ) return(power)
  for( i in consecutive.jump.pt ){
    power.new <-  power[c(i-1,i+1)[which(c(i-1,i+1)>1)]][ which.min(abs(power[i] - power[c(i-1,i+1)[which(c(i-1,i+1)>1)]])) ]
    if( abs(power[i]-power.new) > threshold) power[i] <- power.new
    if( i > 2 & i < (length(power)-1))
      if( (abs(power[i] - power[i-2]) > threshold) & (abs(power[i] - power[i+2]) > threshold) ){
        if( i > 2 ) power[i] <- power[c(i-2,i+2)][ which.min(abs(power[i] - power[c(i-2,i+2)])) ]
        else if( i > 1 ) power[i] <- power[c(i-1,i+2)][ which.min(abs(power[i] - power[c(i-1,i+2)])) ]
      }
  }
  power <- remove.consecutive.jump.pt2( power, threshold, level = level - 1 )
  return(power)
}

remove.consecutive.jump.pt.ind <- function( power, threshold = 10, level = 1 ){
  
  if( level == 0 ) return(power)
  consecutive.jump.pt <- which( abs(diff(power)) > threshold )
  consecutive.jump.pt <- union( consecutive.jump.pt, consecutive.jump.pt + 1 )
  if( length(consecutive.jump.pt) == 0 ) return(power)
  
  ind <- 1:length(power)
  for( i in consecutive.jump.pt ){
    if( (1 < i)&(i < length(power)) ){
      double.diff <- diff(diff(power[c(i-1,i,i+1)]))
      if( double.diff >= 0 ){ind[i] <- ind[i-1]; power[i] <- power[i-1]} 
      else{ind[i]<-ind[i+1]; power[i] <- power[i+1]}
    } 
    if( (2 < i)&(i < length(power)-1) ){
      double.diff <- diff(diff(power[c(i-2,i,i+2)]))
      if( double.diff >= 0 ){
        if( abs(diff(power[c(i-2,i)])) > threshold ){ind[(i-1):i]<-ind[i-2]; power[(i-1):i] <- power[i-2]} 
      } 
      else{
        if( abs(diff(power[c(i+2,i)])) > threshold ){ind[i:(i+1)]<-ind[i+2]; power[i:(i+1)] <- power[i+2]} 
      } 
    }  
  }
  return(ind)
}

find.box.shape.reverse <- function( power, outer.threshold = 10, inner.threshold = 10, level = 2, recursive = T ){
  
  if( (length(power) < 5) | (level==0) ) return(power)
  box.reverse <- power; power.original <- power
  power       <- remove.consecutive.jump.pt2( power, inner.threshold, level = 2 )
  #power[which( power-box.reverse > outer.threshold )] <- box.reverse[which( power-box.reverse > outer.threshold )]
  power.diff  <- round(diff(power))
  jump.position <- which(abs(power.diff) > outer.threshold)
  if( length(jump.position) == 0 ) return(power.original)
  
  increase <- integer()
  decrease <- integer()
  interval <- integer()
  jump.str <- integer()
  jump.end <- integer()
  for( i in 1:length(jump.position) ){
    if( sign(power.diff[jump.position[i]]) == 1){
      if( i == length(jump.position) ){
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.left <- jump.position[i]   - stationary.left
        box.reverse[ (stationary.left+1):length(box.reverse) ] <- box.reverse[stationary.left]
      }else if( sign(power.diff[jump.position[i]]) != sign(power.diff[jump.position[i+1]]) ){
        increase <- c( increase,power.diff[jump.position[i]] )
        decrease <- c( decrease,power.diff[jump.position[i+1]] )
        interval <- c( interval,jump.position[i+1] - jump.position[i] )
        jump.str <- c( jump.str, jump.position[i]+1)
        jump.end <- c( jump.end, jump.position[i+1])
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        rite.ind <- jump.position[i+1]:(jump.position[i+1]+14);  rite.ind <- pmin(rite.ind,length(power))
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F); left.pt <- stationary.left[[2]]
        stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T); rite.pt <- stationary.rite[[2]]
        stationary.left <- jump.position[i]   - stationary.left[[1]]
        stationary.rite <- jump.position[i+1] + stationary.rite[[1]]
        
        subinterval.ind <- c((stationary.left+1):(stationary.rite-1))
        subinterval <- power[subinterval.ind]
        
        test <- find.stationary.state( subinterval,forward=T); remove.left.pt <- test[[1]]
        if( length(remove.left.pt) > 0 && remove.left.pt > 0 ) subinterval[1:remove.left.pt] <- test[[2]]
        test <- find.stationary.state( subinterval,forward=F); remove.rite.pt <- test[[1]]
        if( length(remove.rite.pt) > 0 && remove.rite.pt > 0 ) subinterval[(length(subinterval)+1) - (1:remove.rite.pt)] <- test[[2]]

        box.reverse[ subinterval.ind ] <- min(c(left.pt,rite.pt))
        box.reverse[ box.reverse > power.original ] <- power.original[ box.reverse > power.original ]
#         if( remove.left.pt > 0 ){
#           subinterval.ind <- subinterval.ind[-(1:remove.left.pt)]
#           subinterval     <- subinterval[-(1:remove.left.pt)]
#         } 
#         if( remove.rite.pt > 0 ){
#           subinterval.ind <- subinterval.ind[-(length(subinterval.ind+1)-(1:remove.rite.pt))]
#           subinterval     <- subinterval[-(length(subinterval+1)-(1:remove.rite.pt))]
#         } 
        
        if( recursive ){
          if(level==1){
            box.reverse[ which(box.reverse > power.original) ] <- power.original[ which(box.reverse > power.original) ]
            return( box.reverse )
          } 
          if(length(subinterval)>5)
          if( max(abs(diff(subinterval))) > inner.threshold ){
            while( length(which(abs(diff(subinterval)) > inner.threshold)) > 0 ){
              subinterval.old <- subinterval
              subbox <- ( subinterval - find.box.shape.reverse( subinterval, level = (level-1) ))
              if(subbox[1]>0){
              #  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox - subbox[1], 0 )
                box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox , 0 )
              }else if(subbox[length(subbox)]>0){
              #  box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox - subbox[length(subbox)], 0 )
                box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox , 0 )
              }else{
                box.reverse[subinterval.ind] <- pmax( box.reverse[subinterval.ind] + subbox, 0 )
              }
              subinterval <- subinterval - subbox
              if( all(subinterval-subinterval.old == 0) ) break
            }
          } 
        }
      } 
    }else if( i == 1 & jump.position[1] < length(power)) 
      box.reverse[ 1:(jump.position[1])] <- power[jump.position[1]+1]
  }
  box.reverse[ box.reverse > power.original ] <- power.original[ box.reverse > power.original ]
  return( box.reverse )
}

find.box.shape.reverse2 <- function( active.power, reactive.power, threshold = 10 ){
  
  box.reverse     <- active.power; active.power.original <- active.power
  reactive.power  <- remove.consecutive.jump.pt2( reactive.power, threshold )
  active.power    <- remove.consecutive.jump.pt2(   active.power, threshold )
  power.diff      <- round(diff(reactive.power))
  jump.position   <- which(abs(power.diff) > threshold)
  if( length(jump.position) == 0 ) return(active.power)
  
  increase <- integer()
  decrease <- integer()
  interval <- integer()
  jump.str <- integer()
  jump.end <- integer()
  
  for( i in 1:(length(jump.position)) ){
    if( sign(power.diff[jump.position[i]]) == 1){
      if( i == length(jump.position) ) {
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.left <- jump.position[i]   - stationary.left
        box.reverse[ (stationary.left+1):length(box.reverse) ] <- box.reverse[stationary.left]
      }
      else if( sign(power.diff[jump.position[i]]) != sign(power.diff[jump.position[i+1]]) ){
        
        increase <- c( increase,power.diff[jump.position[i]] )
        decrease <- c( decrease,power.diff[jump.position[i+1]] )
        interval <- c( interval,jump.position[i+1] - jump.position[i] )
        jump.str <- c( jump.str, jump.position[i]+1)
        jump.end <- c( jump.end, jump.position[i+1])
        
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        rite.ind <- jump.position[i+1]:(jump.position[i+1]+14);  rite.ind <- pmin(rite.ind,length(active.power))
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T)[[1]]
        stationary.left <- jump.position[i]   - stationary.left
        stationary.rite <- jump.position[i+1] + stationary.rite + 1
        
        #subinterval <- box.reverse[(jump.position[i]+1):(jump.position[i+1])]
        subinterval <- box.reverse[(stationary.left+1):(stationary.rite-1)]
        remove.pt1  <- find.stationary.state( subinterval, forward=T )[[1]]
        remove.pt2  <- find.stationary.state( subinterval, forward=F )[[1]]
        subinterval.ind <- c((stationary.left+1):(stationary.rite-1))
        if( remove.pt1 > 0 ) subinterval.ind <- subinterval.ind[-c(1:remove.pt1)]
        if( remove.pt2 > 0 ) subinterval.ind <- subinterval.ind[-c((length(subinterval.ind)-remove.pt2+1):length(subinterval.ind))]
        subinterval <- box.reverse[subinterval.ind]
        
        base.line <- min(box.reverse[c(stationary.left,stationary.rite)])
        if( !all(base.line >= subinterval) ){
          box.reverse[(stationary.left+1):(stationary.rite-1)] <- base.line
          while( length(which(abs(diff(subinterval)) > threshold)) > 0 ){
            subbox <- pmax( subinterval - pmax(find.box.shape.reverse( subinterval ), 0), 0 )
            
             if( subbox[1] > 0 )
               subbox[subbox>0] <- subbox[subbox>0] + seq( box.reverse[stationary.left] - base.line - median(head(subbox)), 
                                                           0, length.out=length(which(subbox>0)))
             
             if( subbox[length(subbox)] > 0 )
               subbox[subbox>0] <- subbox[subbox>0] + seq( 0, box.reverse[stationary.rite] - base.line - median(tail(subbox)), 
                                                           length.out=length(which(subbox>0)))
             
            box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox
            subinterval <- subinterval - subbox
            #           if( stationary.left > 0 ) 
            #             box.reverse[ (stationary.left+1):(jump.position[i]) ] <- 
            #             active.power.original[ (stationary.left+1):(jump.position[i]) ] 
            #           if( stationary.rite > 0 )
            #             box.reverse[ (jump.position[i+1]+1):(stationary.rite-1) ] <- 
            #             active.power.original[ (jump.position[i+1]+1):(stationary.rite-1) ]
            if( all(subbox==0) ) break
          }
          if( remove.pt1 > 0 )
            box.reverse[ stationary.left+c(1:remove.pt1) ] <- seq( box.reverse[stationary.left], 
                                                                   box.reverse[stationary.left+remove.pt1+1], 
                                                                   length.out=length(remove.pt1))
          
          if( remove.pt2 > 0 )
            box.reverse[ stationary.rite-c(remove.pt2:1) ] <- seq( box.reverse[stationary.rite-remove.pt2-1], 
                                                                   box.reverse[stationary.rite], 
                                                                   length.out=length(remove.pt2))
          
        }else{
          increase <- increase[-length(increase)]
          decrease <- decrease[-length(decrease)]
          interval <- interval[-length(interval)]
          jump.str <- jump.str[-length(jump.str)]
          jump.end <- jump.end[-length(jump.end)]
        }
      }
    }else if( i == 1 & jump.position[1] < length(power)){
      rite.ind <- jump.position[i]:(jump.position[i]+14);  rite.ind <- pmin(rite.ind,length(power))
      stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T)[[1]]
      stationary.rite <- jump.position[i] + stationary.rite + 1
      box.reverse[ 1:(stationary.rite-1)] <- box.reverse[stationary.rite]
    } 
  }

box.reverse[ box.reverse > active.power.original ] <- active.power.original[ box.reverse > active.power.original ]
  return( list( box.reverse, data.frame( increase, decrease, interval, jump.str, jump.end) ))
}

find.box.shape.reverse4 <- function( active.power, reactive.power, threshold = 10, 
                                     reactive.r.edge.parameters,
                                     reactive.f.edge.parameters,
                                     reactive.duration.parameters, 
                                     debug.mode = F,
                                     is.rp.zero = F 
                                     ){

  if( is.rp.zero ){
    reactive.power.original <- reactive.power
    reactive.power          <- active.power
  }
  
  # reactive.r.edge.parameters['rp_height']가 음수인 경우에는 -1 곱하기
  upside.down  <- F
  if( all( active.power == reactive.power )){
    if( sign(reactive.r.edge.parameters['rp_height']) == -1 ){
      active.power <- active.power * (-1)
      upside.down  <- T
    }
  }
  
  box.reverse     <- active.power; active.power.original <- active.power
  reactive.power  <- remove.consecutive.jump.pt2( reactive.power, threshold )
  active.power    <- remove.consecutive.jump.pt2(   active.power, threshold )
  power.diff      <- diff(reactive.power)
  jump.position   <- which(abs(power.diff) > threshold)
  jump.value      <- power.diff[jump.position]
  
#   remove.pt <- apply( outer( which(diff(jump.position) == 1), c(0,1), '+'), 1,
#                       function(x) x[which.min(abs(jump.value[x]))] ) 
#   if( length(remove.pt) > 0 ){
#     remaining.pt <- apply( outer( which(diff(jump.position) == 1), c(0,1), '+'), 1, 
#                            function(x) x[-which.min(abs(jump.value[x]))] ) 
#     
#     jump.value[remaining.pt] <- jump.value[remaining.pt] + jump.value[remove.pt]
#     jump.value    <- jump.value[-remove.pt]
#     jump.position <- jump.position[-remove.pt]
#   }
   
  if( length(jump.position) <= 1 ) return(NULL) 
  
  if( is.rp.zero ){
    remove.pt     <- sapply( jump.position, function(i) find.jump.value( reactive.power.original, i))
    remove.pt     <- abs(unlist(remove.pt[sapply( remove.pt, length ) > 0])) < 35
    jump.position <- jump.position[ remove.pt ]
    jump.value    <- jump.value[ remove.pt ]
  }
  
  if( length(jump.position) <= 1 ) return(NULL) 
  
  r.lower.bound <- qnorm( .01, reactive.r.edge.parameters[1], reactive.r.edge.parameters[2] ) 
  r.upper.bound <- qnorm( .99, reactive.r.edge.parameters[1], reactive.r.edge.parameters[2] )
  f.lower.bound <- qnorm( .01, reactive.f.edge.parameters[1], reactive.f.edge.parameters[2] )
  f.upper.bound <- qnorm( .99, reactive.f.edge.parameters[1], reactive.f.edge.parameters[2] )
r.lower.bound <- r.lower.bound * (1-sign(r.lower.bound) * .2)
r.upper.bound <- r.upper.bound * (1+sign(r.upper.bound) * .2)
f.lower.bound <- f.lower.bound * (1-sign(f.lower.bound) * .2)
f.upper.bound <- f.upper.bound * (1+sign(f.upper.bound) * .2)

  r.edge <- which( (r.lower.bound < jump.value) & (jump.value < r.upper.bound) )
  if( length(r.edge) == 0 ) return( NULL)
  
  fn <- function(i,j){
    edges <- (sign(jump.value[i])*jump.value < 0) & 
      (jump.position[i] < jump.position) & (jump.position<jump.position[j])
    if( !any(edges) ) return(NULL)
    prob1 <-  dnorm( jump.value[edges], reactive.f.edge.parameters[1], reactive.f.edge.parameters[2])
    if( length(reactive.duration.parameters) > 1 ){
      if(reactive.duration.parameters[2] > 1000){
        prob2 <-  1. / ( jump.position[ edges ] - (jump.position[i] + reactive.duration.parameters[1]))
        prob2 <- abs(prob2)
      }else{
        prob2 <-  dnorm( jump.position[ edges ], 
                         jump.position[i] + reactive.duration.parameters[1], reactive.duration.parameters[2])
      }
    }else prob2 <- 1 
    score <- prob1 * prob2
    if( all(score==0) ) return(NULL)
    score.max <- which.max(score)
    return( which( jump.position ==( jump.position[ edges ][score.max])))
  }
  
  tmp <- matrix( c(r.edge,c(tail(r.edge,-1),NA)), ncol=2 ); tmp[is.na(tmp)] <- length(jump.position)
  f.edge <- apply( tmp, 1, function(x) fn(x[1],x[2]))
  r.edge <- r.edge[!sapply( f.edge, is.null)]
  f.edge <- unlist(f.edge)
  if( length(f.edge) == 0 ) return( NULL)
  
  
  if( debug.mode ){
    
    par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
    plot(   active.power, type='l', ann=FALSE,   xaxt="n") 
    abline( v=jump.position[r.edge], col='red' ) 
    abline( v=jump.position[f.edge], col='blue')
    plot( reactive.power, type='l', ann=FALSE,   xaxt="n") 
    abline( v=jump.position[r.edge], col='red' ) 
    abline( v=jump.position[f.edge], col='blue')
    title(paste('Step1 : Edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)), 
          outer=TRUE)
    mtext("Timestamp", 1, 1, outer=TRUE)
    mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
    par(mfrow=c(1,1))

  }
  edge <- data.frame( re = jump.position[r.edge], fe = jump.position[f.edge] )
  if(any( edge$re == edge$fe )) edge <- edge[-which( edge$re == edge$fe ),]
  if( nrow(edge) == 0 ) return(NULL)
  
  for( i in 1:(nrow(edge)) ){
        left.ind <- (edge$re[i]-14):edge$re[i];  left.ind <- pmax(left.ind,1)
        rite.ind <- edge$fe[i]:(edge$fe[i]+14);  rite.ind <- pmin(rite.ind,length(active.power))
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T)[[1]]
        stationary.left <- edge$re[i] - stationary.left
        stationary.rite <- min( edge$fe[i] + stationary.rite + 1, length(box.reverse))
        
        subinterval <- box.reverse[(stationary.left+1):(stationary.rite-1)]
        remove.pt1  <- find.stationary.state( subinterval, forward=T )[[1]]
        remove.pt2  <- find.stationary.state( subinterval, forward=F )[[1]]
        subinterval.ind <- c((stationary.left+1):(stationary.rite-1))
        if( remove.pt1 > 0 ) subinterval.ind <- subinterval.ind[-c(1:remove.pt1)]
        if( remove.pt2 > 0 ) subinterval.ind <- subinterval.ind[-c((length(subinterval.ind)-remove.pt2+1):length(subinterval.ind))]
        subinterval <- box.reverse[subinterval.ind]
         
        base.line <- min(box.reverse[c(stationary.left,stationary.rite)])
        if( !all(base.line >= subinterval) ){
          box.reverse[(stationary.left+1):(stationary.rite-1)] <- base.line
          while( length(which(abs(diff(subinterval)) > threshold)) > 0 ){
            subbox <- pmax( subinterval - find.box.shape.reverse( subinterval ), 0 )
            
#             if( subbox[1] > 0 ){
#               box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox + 
#                 box.reverse[stationary.left] - base.line - subbox[1]
#               
#               #subbox[subbox>0] <- subbox[subbox>0] +  box.reverse[stationary.left] - base.line - subbox[1]
#               #seq( box.reverse[stationary.left] - base.line - median(head(subbox)), 0, length.out=length(which(subbox>0)))
#               
#             }else if( subbox[length(subbox)] > 0 ){
#               box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox + 
#                 box.reverse[stationary.rite] - base.line - subbox[length(subbox)]
#               
#               #subbox[subbox>0] <- subbox[subbox>0] + box.reverse[stationary.rite] - base.line - subbox[1]
#               #seq( 0, box.reverse[stationary.rite] - base.line - median(tail(subbox)), length.out=length(which(subbox>0)))
#               
#             }else{
              box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox  
#             }
            
            subinterval <- subinterval - subbox
            
            subinterval.lower <- envelopeDetector( data.frame(active_power = subinterval,
                                                              reactive_power=0), use.active.power=T )$LowerEnvelope
            if(!all((subinterval-subinterval.lower)==0)){
              box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subinterval-subinterval.lower
              subinterval <- subinterval.lower
            }
              
            #           if( stationary.left > 0 ) 
            #             box.reverse[ (stationary.left+1):(jump.position[i]) ] <- 
            #             active.power.original[ (stationary.left+1):(jump.position[i]) ] 
            #           if( stationary.rite > 0 )
            #             box.reverse[ (jump.position[i+1]+1):(stationary.rite-1) ] <- 
            #             active.power.original[ (jump.position[i+1]+1):(stationary.rite-1) ]
            if( all(subbox==0) ) break
          }
          if( remove.pt1 > 0 ){
            if( length(remove.pt1) == 1){
              box.reverse[ stationary.left+1 ] <- box.reverse[stationary.left]
            }else{
              box.reverse[ stationary.left+c(1:remove.pt1) ] <- seq( box.reverse[stationary.left], 
                                                                     box.reverse[stationary.left+remove.pt1+1], 
                                                                     length.out=length(remove.pt1))
            }
          }
            
          
          if( remove.pt2 > 0 ){
            if( length(remove.pt2) == 1 ){
              box.reverse[ stationary.rite-1 ] <- box.reverse[stationary.rite]
            }else{  
              box.reverse[ stationary.rite-c(remove.pt2:1) ] <- seq( box.reverse[stationary.rite-remove.pt2-1], 
                                                                     box.reverse[stationary.rite], 
                                                                     length.out=length(remove.pt2))
            }
          }
box.reverse[(stationary.left):(stationary.rite)] <- box.reverse[(stationary.left):(stationary.rite)] + 
  pmin( (median(subinterval) - base.line) * 2, 0 )

          box.reverse[ box.reverse > active.power.original ] <- active.power.original[ box.reverse > active.power.original ]
          
        }
      }
  
   box.reverse[ box.reverse > active.power.original ] <- active.power.original[ box.reverse > active.power.original ]

  box.shape <- active.power - box.reverse
  box.lists <- series.to.box.lists(1:length(box.shape), box.shape, 0 )[[1]]
if( is.null(box.lists) ) return(active.power)

  box.lists <- box.lists[ which(box.lists$end - box.lists$str > 5), ]
  if( is.null(box.lists) | is.null(nrow(box.lists)) | nrow(box.lists) == 0 ) return( active.power )

  box.reverse[ 1:(box.lists$str[1]-1) ] <- active.power.original[ 1:(box.lists$str[1]-1) ]
  if( nrow(box.lists) > 1){
    for( i in 1:(nrow(box.lists)-1)){
      box.reverse[ (box.lists$end[i]+1):(box.lists$str[i+1]-1) ] <- 
        active.power.original[ (box.lists$end[i]+1):(box.lists$str[i+1]-1) ]
    }
  }
  box.reverse[ (box.lists$end[nrow(box.lists)]+1):length(box.reverse) ] <- 
  active.power.original[ (box.lists$end[nrow(box.lists)]+1):length(box.reverse) ]

  #box.lists <- subset(box.lists, (r.lower.bound > mode.p) | (mode.p > r.upper.bound))
  #for( i in 1:nrow(box.lists) ){
  #  box.reverse[ box.lists$str[i]:box.lists$end[i] ] <- ap[ box.lists$str[i]:box.lists$end[i] ]
  #}

if( upside.down ) box.reverse <- box.reverse * (-1)  
  
   return( list( 'residual' = box.reverse[!is.na(box.reverse)], 'edge' = edge) )
}

find.box.shape.reverse3 <- function( active.power, reactive.power, threshold = 10, 
                                     reactive.r.edge.threshold,
                                     reactive.f.edge.threshold,
                                     reactive.duration.threshold 
){
  
  box.reverse     <- active.power; active.power.original <- active.power
  reactive.power  <- remove.consecutive.jump.pt2( reactive.power, threshold )
  active.power    <- remove.consecutive.jump.pt2(   active.power, threshold )
  power.diff      <- round(diff(reactive.power))
  jump.position   <- which(abs(power.diff) > threshold)
  r.edge <- which( ( power.diff > reactive.r.edge.threshold[1]) & ( power.diff < reactive.r.edge.threshold[2]) )
  f.edge <- which( (-power.diff > reactive.f.edge.threshold[1]) & (-power.diff < reactive.f.edge.threshold[2]) )
  print(r.edge)
  print(f.edge)
  edge <- list()
  for( i in r.edge ){
    print(i)
    left <- (i + reactive.duration.threshold[1])
    rite <- (i + reactive.duration.threshold[2])
    if( any((left<f.edge)&(f.edge<rite)) ){
      edge[[length(edge)+1]] <- c( i, median(f.edge[which((left<f.edge)&(f.edge<rite))]) )
    }
  }
  if( length(edge) == 0 ) return(active.power)
  
  for( i in 1:(length(edge)) ){
    
    left.ind <- (edge[[i]][1]-14):edge[[i]][1];  left.ind <- pmax(left.ind,1)
    rite.ind <- edge[[i]][2]:(edge[[i]][2]+14);  rite.ind <- pmin(rite.ind,length(active.power))
    stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
    stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T)[[1]]
    stationary.left <- edge[[i]][1] - stationary.left
    stationary.rite <- edge[[i]][2] + stationary.rite + 1
    
    subinterval <- box.reverse[(stationary.left+1):(stationary.rite-1)]
    remove.pt1  <- find.stationary.state( subinterval, forward=T )[[1]]
    remove.pt2  <- find.stationary.state( subinterval, forward=F )[[1]]
    subinterval.ind <- c((stationary.left+1):(stationary.rite-1))
    if( remove.pt1 > 0 ) subinterval.ind <- subinterval.ind[-c(1:remove.pt1)]
    if( remove.pt2 > 0 ) subinterval.ind <- subinterval.ind[-c((length(subinterval.ind)-remove.pt2+1):length(subinterval.ind))]
    subinterval <- box.reverse[subinterval.ind]
    
    base.line <- min(box.reverse[c(stationary.left,stationary.rite)])
    if( !all(base.line >= subinterval) ){
      box.reverse[(stationary.left+1):(stationary.rite-1)] <- base.line
      while( length(which(abs(diff(subinterval)) > threshold)) > 0 ){
        subbox <- pmax( subinterval - pmax(find.box.shape.reverse( subinterval ), 0), 0 )
        
        #             if( subbox[1] > 0 ){
        #               box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox + 
        #                 box.reverse[stationary.left] - base.line - subbox[1]
        #               
        #               #subbox[subbox>0] <- subbox[subbox>0] +  box.reverse[stationary.left] - base.line - subbox[1]
        #               #seq( box.reverse[stationary.left] - base.line - median(head(subbox)), 0, length.out=length(which(subbox>0)))
        #               
        #             }else if( subbox[length(subbox)] > 0 ){
        #               box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox + 
        #                 box.reverse[stationary.rite] - base.line - subbox[length(subbox)]
        #               
        #               #subbox[subbox>0] <- subbox[subbox>0] + box.reverse[stationary.rite] - base.line - subbox[1]
        #               #seq( 0, box.reverse[stationary.rite] - base.line - median(tail(subbox)), length.out=length(which(subbox>0)))
        #               
        #             }else{
        box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subbox  
        #             }
        
        subinterval <- subinterval - subbox
        
        subinterval.lower <- envelopeDetector( data.frame(active_power = subinterval,
                                                          reactive_power=0), use.active.power=T )$LowerEnvelope
        if(!all((subinterval-subinterval.lower)==0)){
          box.reverse[subinterval.ind] <- box.reverse[subinterval.ind] + subinterval-subinterval.lower
          subinterval <- subinterval.lower
        }
        
        #           if( stationary.left > 0 ) 
        #             box.reverse[ (stationary.left+1):(jump.position[i]) ] <- 
        #             active.power.original[ (stationary.left+1):(jump.position[i]) ] 
        #           if( stationary.rite > 0 )
        #             box.reverse[ (jump.position[i+1]+1):(stationary.rite-1) ] <- 
        #             active.power.original[ (jump.position[i+1]+1):(stationary.rite-1) ]
        if( all(subbox==0) ) break
      }
      if( remove.pt1 > 0 )
        box.reverse[ stationary.left+c(1:remove.pt1) ] <- seq( box.reverse[stationary.left], 
                                                               box.reverse[stationary.left+remove.pt1+1], 
                                                               length.out=length(remove.pt1))
      
      if( remove.pt2 > 0 )
        box.reverse[ stationary.rite-c(remove.pt2:1) ] <- seq( box.reverse[stationary.rite-remove.pt2-1], 
                                                               box.reverse[stationary.rite], 
                                                               length.out=length(remove.pt2))
      box.reverse[ box.reverse > active.power.original ] <- active.power.original[ box.reverse > active.power.original ]
      
    }
  }
  
  box.reverse[ box.reverse > active.power.original ] <- active.power.original[ box.reverse > active.power.original ]
  return( box.reverse )
}

box.detect <- function( data, feederID, box.order = c(rep('reactive',4),rep('active',3)) ){
  
  original.p    <- data$active_power[ data$feeder_id == feederID]
  original.q    <- data$reactive_power[ data$feeder_id == feederID]
  active.base   <- max( min(data$active_power[ data$feeder_id == feederID]), 0 )
  data$active_power[ data$feeder_id == feederID] <- 
    data$active_power[ data$feeder_id == feederID] - active.base
  
  print('high.power')
  data$high.power[ data$feeder_id == feederID] <- 
    pmax( data$active_power[ data$feeder_id == feederID] - 
            find.high.power.reverse( data$active_power[ data$feeder_id == feederID], outer.threshold=550), 0 )
  data$active_power[ data$feeder_id == feederID] <- 
    data$active_power[ data$feeder_id == feederID] - data$high.power[ data$feeder_id == feederID]
  
  print('peak')
  if( 'envelope.active' %in% names(data) ){
    envelope.active <- data$envelope.active[ data$feeder_id == feederID,]
  }else envelope.active   <- envelopeDetector(data[ data$feeder_id == feederID,])
  if( 'envelope.reactive' %in% names(data) ){
    envelope.reactive <- data$envelope.reactive[ data$feeder_id == feederID,]
  }else envelope.reactive <- envelopeDetector(data[ data$feeder_id == feederID,], use.active.power=F)
    
  print('box')
  active.power  <- list()
  active.box    <- list()
  active.power[[1]]    <- envelope.active$active_power
  active.power[[2]]    <- envelope.active$LowerEnvelope
  active.box[[1]]      <- active.power[[1]] - active.power[[2]]
  
  if( any('reactive' %in% box.order) ){
    reactive.power  <- list()
    reactive.box    <- list()
    reactive.power[[1]]  <- envelope.reactive$reactive_power
    reactive.power[[2]]  <- envelope.reactive$UpperEnvelope
    reactive.box[[1]]    <- reactive.power[[1]] - reactive.power[[2]]
  } 
  
  # main algorithm
  for( i in 1:length(box.order) ){
    if( box.order[i] == 'reactive' ){
      nap <- length(  active.power) 
      nrp <- length(reactive.power)
      if( length(diff(reactive.power[[nrp]])[ abs(diff(reactive.power[[nrp]])) > 10]) > 0 ) {
        active.power[[nap+1]]   <- find.box.shape.reverse2(active.power[[nap]], reactive.power[[nrp]])[[1]]
        reactive.power[[nrp+1]] <- find.box.shape.reverse(reactive.power[[nrp]])
        active.box[[nap]]       <- pmax(  active.power[[nap]] -   active.power[[nap+1]], 0)
        reactive.box[[nrp]]     <- pmax(reactive.power[[nrp]] - reactive.power[[nrp+1]], 0)
        active.power[[nap+1]]   <- active.power[[nap]] - active.box[[nap]]
        reactive.power[[nrp+1]] <- reactive.power[[nrp]] - reactive.box[[nrp]]
        if( all( active.box[[nap]] == 0 ) ){
          active.power[[nap+1]] <- NULL
          active.box[[nap]] <- NULL
        }
      }        
    }else{
      nap <- length(  active.power )
      if( length(diff(active.power[[nap]])[ abs(diff(active.power[[nap]])) > 10]) > 0 ) {
        active.power[[nap+1]] <- find.box.shape.reverse(active.power[[nap]])
        active.box[[nap]]     <- pmax(active.power[[nap]] -   active.power[[nap+1]],0)
        active.power[[nap+1]] <- active.power[[nap]] - active.box[[nap]]
      }
    }
  }
    
  baseline.usage <- active.base * length(active.power[[1]])
  box.usage <- unlist(lapply( active.box, sum ))
  percent <- c( baseline.usage, box.usage, sum(active.power[[length(active.power)]]) ) * 
    length(data$active_power[ data$feeder_id == feederID]) / length(active.power[[1]])
  percent <- percent / sum(original.p) * 100
  
  if( TRUE ){
    df <- data.frame( timestamp = envelope.active$timestamp )
    df$original  <- original.p
    df$base      <- active.base
    df$high.power <- data$high.power[ data$feeder_id == feederID]
    df <- cbind( df, as.data.frame( active.box ) )
    names(df)[-c(1:4)] <- paste0('Box',1:length(active.box))
    df$residual <- active.power[[length(active.power)]]
    
    if(length(which(colSums(df == 0) == nrow(df)))>0) df <- df[,-which(colSums(df == 0) == nrow(df))]
    box.index <- rbind.fill( lapply( as.numeric(gsub('Box','',names(df)[grep('Box',names(df))])), 
                                     function(x) data.frame( series.to.box.lists( df$timestamp, 
                                                                                  df[,paste0('Box',x)], 
                                                                                  0)[[1]],
                                                             nBox = x) ))
    # Box 별로 time-series t, p, q
    box.lists <- lapply( as.numeric(gsub('Box','',names(df)[grep('Box',names(df))])), 
                         function(x) series.to.box.lists( df$timestamp, 
                                                          df[,paste0('Box',x)], 0)[[2]] )
    
    # Box의 최대값이 threshold(500) 초과하는 값은 high power로 분류
    if( any(box.index$max.p > 500) ){
      high.power <- adply( box.index[ box.index$max.p > 500, c('nBox','no')], 1, 
                           function(x) box.lists[[x$nBox]][[x$no]] )
      high.power <- ddply( high.power, .(timestamp), summarize, p = sum(p) )
      high.power <- merge( high.power, data.frame( timestamp = df$timestamp), all.y = T )
      high.power$p[ is.na(high.power$p) ] <- 0  
      high.power.logical <- T
      
      box.index <- box.index[ which(box.index$max.p > 500), ]
      for( i in unique(box.index$nBox) ){
        tmp <- subset( box.index, nBox == i )
        df[unlist(sapply( 1:nrow(tmp), function(x) tmp$str[x]:tmp$end[x] )),paste0('Box',i)] <- 0
      }
      df$high.power <- df$high.power + high.power$p
    }else high.power.logical <- F
    df$original.q <- original.q
    
    df.melt <- melt( df, id.vars='timestamp')
    
    print('figure')
    fig <- NULL
    if( FALSE ){
    fig <- ggplot(df.melt,aes(x=timestamp,y=value,colour=variable)) + 
      geom_line() + 
      facet_grid(variable~., scales='free') +
      theme( legend.position='NONE') +
      ggtitle( paste('siteID :',unique(data$site_id),
                     '- Day :',unique(as.Date(data$timestamp)),
                     '- feeder :',feederID,
                     '\n',paste0( round(ddply( df.melt, .(variable), summarize, value=sum(value))$value / 3600), collapse='/')) )
    
    #     ggsave( plot=fig, path='fig', 
    #             filename=paste0('siteID',unique(data$site_id),
    #                             'Day',unique(as.Date(data$timestamp)),
    #                             'Feeder',feederID,'.png' ),
    #             width = 50, height = 20, limitsize=F)
    }
    if( any('reactive' %in% box.order) ){
      df.reactive <- as.data.frame( reactive.box )
      names(df.reactive) <- paste0('reactive.box',1:length(reactive.box))      
      df <- cbind( df, df.reactive )
      nrp <- length(reactive.box)
      nap <- length(  active.box)
      if( nrp < nap ) df[,paste0('reactive.box',(nrp+1):nap)] <- 0
      df$reactive.residual <- reactive.power[[length(reactive.power)]]
    }    
  }
  return(list(df,fig))
}


# Box 별 대표값
representative.value <- function(x){
  table.p <- table(x$p)
  table.q <- table(x$q)
  return( data.frame('mode.p' = median(as.numeric(names(table.p)[table.p == max(table.p)])),
                     'mode.q' = median(as.numeric(names(table.q)[table.q == max(table.q)])),
                     'min.p'  = min(x$p),
                     'min.q'  = min(x$q),
                     'max.p'  = max(x$p),
                     'max.q'  = max(x$q),
                     'median.p'  = median(x$p),
                     'median.q'  = median(x$q),
                     'sum.p'  = sum(x$p), 
                     'sum.q'  = sum(x$q)))
}

shift <- function(vec,n) 
  if (n == 0){
    vec
  }else if (n > 0){
    c(rep(NA, n), head(vec, -n))
  }else if (n < 0){
    c(tail(vec, n), rep(NA, -n))
  }
    


# 유효/무효 전력 등을 이용해서 나눈 timeseries p,q 값을 개별 Box 단위로 나누기
series.to.box.lists <- function( series.time, series.p, series.q  = 0 ){
  
  series <- data.frame( timestamp = series.time, 
                        p = pmax( series.p, 0 ), 
                        q = series.q ) 
  n <- nrow(series)
  if( sum(series$p) == 0 ) return(NULL)
  
  str <- which( diff(series$p) > 0 & (series$p[-n]==0) ) + 1
  end <- which( diff(series$p) < 0 & (series$p[-1]==0) )  
  if( series$p[1] != 0 ) str <- c(1,str)
  if( series$p[n] != 0 ) end <- c(end,n)
  
  duration <- end - str + 1
  
  box.lists <- apply( data.frame(str,end), 1, function(x){ series[x['str']:x['end'],] } ) 
  result1 <- cbind( data.frame(str,end,duration), rbind.fill(lapply( box.lists, representative.value )) )
  result1$no <- 1:nrow(result1)
  result2 <- box.lists
  return( list( result1, result2 ) ) 
}

box.clustering <- function( data, feederID, box.order = c(rep('reactive',4),rep('active',3)), nCluster = 6 ){
  
  box.detect.results <- box.detect( data, feederID, box.order )
  
  # 겯과 그림으로 그리기
  print( box.detect.results[[2]] )

  df <- box.detect.results[[1]]
  # Box 별로 시작지점, 끝지점, 대표값, duration 등 찾기 
  box.index <- rbind.fill( lapply( as.numeric(gsub('Box','',names(df)[grep('Box',names(df))][-1])), 
                                   function(x) data.frame( series.to.box.lists( df$timestamp, 
                                                                                df[,paste0('Box',x)], 
                                                                                df[,paste0('reactive.box',x)])[[1]],
                                                           nBox = (x-1)) ))
  box.index$duration <- box.index$end - box.index$str + 1
  
  # Box 별로 time-series t, p, q
  box.lists <- lapply( as.numeric(gsub('Box','',names(df)[grep('Box',names(df))][-1])), 
                       function(x) series.to.box.lists( df$timestamp, 
                                                        df[,paste0('Box',x)], 
                                                        df[,paste0('reactive.box',x)])[[2]] )
  
  # Box의 최대값이 threshold(500) 초과하는 값은 high power로 분류
  if( any(box.index$max.p > 500) ){
    high.power <- adply( box.index[ box.index$max.p > 500, c('nBox','no')], 1, 
                         function(x) box.lists[[x$nBox]][[x$no]] )
    high.power <- ddply( high.power, .(timestamp), summarize, p = sum(p) )
    high.power <- merge( high.power, data.frame( timestamp = df$timestamp), all.y = T )
    high.power$p[ is.na(high.power$p) ] <- 0  
    high.power.logical <- T
  }else high.power.logical <- F
  box.index <- box.index[ which(box.index$max.p <= 500), ]
  
  # duration이 5초 미만인 경우 short peak로 분류
  peak.power <- adply( box.index[ box.index$duration < 5, c('nBox','no')], 1, 
                       function(x) box.lists[[x$nBox]][[x$no]] )
  peak.power <- ddply( peak.power, .(timestamp), summarize, p = sum(p) )
  peak.power <- merge( peak.power, data.frame( timestamp = df$timestamp), all.y = T )
  if( nrow(peak.power) > 0 ){
    if( length(is.na(peak.power$p)) > 0 ) peak.power$p[ is.na(peak.power$p) ] <- 0
    box.index <- box.index[ box.index$duration >= 5, ]
    peak.power.logical = T  
  }else peak.power.logical = F
  
  
  box.cluster.data.frame <- function( cluster.no ){
    box.cluster <-  adply( box.index[ box.index$cluster == cluster.no, c('nBox','no')], 1, 
                           function(x) box.lists[[x$nBox]][[x$no]] )
    box.cluster <- ddply( box.cluster, .(timestamp), summarize, p = sum(p))
    box.cluster <- merge( box.cluster, data.frame( timestamp = df$timestamp), all.y = T )
    box.cluster$p[ is.na(box.cluster$p) ] <- 0
    names(box.cluster)[ names(box.cluster) == 'p' ] <- paste0('Cluster',cluster.no)
    return(box.cluster)
  }
  
  # high power와 short peak를 제외한 Box에 대해 대표값 등을 이용해서 clustering
  box.index$cluster <- kmeans( box.index[,c('mode.p','max.p','min.p','duration')], centers=nCluster)$cluster
  seperatedBox.lists <- lapply( 1:nCluster, box.cluster.data.frame )
  merge.all <- function(x, y) merge(x, y, all=TRUE, by="timestamp")
  seperatedBox <- Reduce(merge.all, seperatedBox.lists)
  tmp          <- data.frame( timestamp = df$timestamp, 
                              original = df$original, 
                              base = df$base )
  if( high.power.logical ) tmp$High.power <- high.power$p
  if( peak.power.logical ) tmp$short.peak <- peak.power$p
  seperatedBox <- cbind( tmp, data.frame( peak = df$Box1,
                                                   seperatedBox[,names(seperatedBox) !='timestamp'], 
                                                   residual = df$residual ))
  
  seperatedBox.melt <- melt( seperatedBox, id.vars = 'timestamp' )
  results <- ddply( seperatedBox.melt, .(variable), summarize, value=sum(value))
  results <- results$value / sum(seperatedBox$original) * 100
  results <- paste0(round(results,2),collapse='/')
  
  # 그래프 출력
  fig <- ggplot( seperatedBox.melt, aes(x=timestamp, y=value, colour=variable)) + 
    geom_line() +
    facet_grid(variable~., scales='free') +
    ggtitle( paste( results,'\n')  ) +
    theme(legend.position='NONE',strip.text.y = element_text(angle =   0),legend.position='NONE')

  print(fig)

  return(append(box.detect.results,list(seperatedBox,fig)))
  
}

find.high.power.reverse <- function( power, outer.threshold = 500, inner.threshold = 10 ){
  
  if( length(power) < 5 ) return(power)
  box.reverse <- power; power.original <- power
  power       <- remove.consecutive.jump.pt2( power, inner.threshold )
  power.diff  <- round(diff(power))
  jump.position <- which(abs(power.diff) > outer.threshold)
  if( length(jump.position) == 0 ) return(power)
  
  increase <- integer()
  decrease <- integer()
  interval <- integer()
  jump.str <- integer()
  jump.end <- integer()
  for( i in 1:length(jump.position) ){
    if( sign(power.diff[jump.position[i]]) == 1){ # up
      if( i == length(jump.position) ){
        left.ind <- (jump.position[i]-14):jump.position[i]; left.ind <- pmax(left.ind,1)
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.left <- jump.position[i]   - stationary.left
        box.reverse[(stationary.left+1):length(box.reverse) ] <- box.reverse[c(stationary.left)] 
      }
      else if( sign(power.diff[jump.position[i]]) != sign(power.diff[jump.position[i+1]]) ){
        increase <- c( increase,power.diff[jump.position[i]] )
        decrease <- c( decrease,power.diff[jump.position[i+1]] )
        interval <- c( interval,jump.position[i+1] - jump.position[i] )
        jump.str <- c( jump.str, jump.position[i]+1)
        jump.end <- c( jump.end, jump.position[i+1])
        left.ind <- (jump.position[i]-14):jump.position[i];      left.ind <- pmax(left.ind,1)
        rite.ind <- jump.position[i+1]:(jump.position[i+1]+14);  rite.ind <- pmin(rite.ind,length(power))
        stationary.left <- find.stationary.state( box.reverse[left.ind],forward=F)[[1]]
        stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T)[[1]]
        stationary.left <- jump.position[i]   - stationary.left
        stationary.rite <- jump.position[i+1] + stationary.rite
        
        subinterval <- box.reverse[(jump.position[i]+1):(jump.position[i+1])]
        box.reverse[(stationary.left+1):(stationary.rite-1)] <- min(box.reverse[c(stationary.left,stationary.rite)])
        if( length(abs(diff(subinterval)) > inner.threshold) > 0 ){
          subbox <- subinterval - pmax(find.box.shape.reverse( subinterval ), 0)
          box.reverse[(jump.position[i]+1):(jump.position[i+1])] <- 
            box.reverse[(jump.position[i]+1):(jump.position[i+1])] + subbox
        }
      } 
    }else if( i == 1 & jump.position[1] < length(power)){
      rite.ind <- jump.position[i+1]:(jump.position[i+1]+14);  rite.ind <- pmin(rite.ind,length(power))
      stationary.rite <- find.stationary.state( box.reverse[rite.ind],forward=T)[[1]]
      stationary.rite <- jump.position[i+1] + stationary.rite
      box.reverse[1:(stationary.rite-1)] <- box.reverse[c(stationary.rite)] 
    }
  }
  box.reverse[ box.reverse > power.original ] <- power.original[ box.reverse > power.original ]
  return( box.reverse )
}




