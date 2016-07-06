showPatternScan.ricecooker_15Hz <- function (data) {
  
  # exam
  if ( as.numeric(data$timestamp[2]) - as.numeric(data$timestamp[1]) > 0.5 ) {
    stop('Make sure that the input data is from 15Hz signal!')
  }
  
  '%nin%'<- Negate('%in%')
  standby.pwr <- 1
  
  # warming_info
  ap.h1.min <- 75 # 90
  ap.h1.max <- 115 # 102.362
  esti_usage.sec <- 4
  esti_usage.watt <- 84.2
  sec.period <- 30
  sec.period.margin <- 5
  time.slot.sec <- 30
  time.slot_margin <- 6
  
  
  # cooking_info
  c.ap.h1.min <- 710 # 756.215
  c.ap.h1.max <- 900 # 862.313
  c.cooking.time <- 3000 # 40 min. approx. cooking 2 (2400) ~ (cooking3) 50 min. approx (3000)
  c.cooking.watt <- 431.604 # cooking3 :  # cooking 2 (402.06) ,, (cooking3) 431.604 ,, 
  c.min.fluc.num <- 7
  c.t.boundary_minute <- 40
  
  # parameters
  w.eff_size <- 7
  c.eff_size <- 7
  
  thres_rp.delta <- 30
  
  w.min_minute <- 4.5 # to split into lumps
  
  c.min_t.minute <- 1
  
  # parameters
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"), 
                                            floor_date(max(data$timestamp), unit = "second"), 'secs'))
  answer.log <- data.frame( answer.log, p= rep(standby.pwr,nrow(answer.log)), q= rep(0,nrow(answer.log)))
  
  # examine raw data
  # data <- data[order(data$timestamp),] 
  
  s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
  s.pattern <- s.pattern[order(s.pattern$start.idx),]
  
  if (nrow(s.pattern) == 0) return(answer.log)
  
  # step) warming signal computation
  
  # medium-power signals for warming
  s.pattern.warming1 <- subset(s.pattern, (h1 >= ap.h1.min) & (h1 <= ap.h1.max) & (abs(sub.delta) <= thres_rp.delta) & !is.na(min.slope) ) # original set
  s.pattern.warming2 <- subset(s.pattern, (h1 >= ap.h1.min) & (h1 <= ap.h1.max) & is.na(min.slope) ) # added set
  
  s.pattern.warming <- rbind(s.pattern.warming1, s.pattern.warming2)
  
  if (nrow(s.pattern.warming) < w.eff_size){
    
    print("There's no SIGNAL for warming mode.")
  } else {
    
    chosen.s.info <- s.pattern.warming[order(s.pattern.warming$start.idx),] 
    lump.idx <- c(0, which(diff(as.numeric(chosen.s.info$start.timestamp)) >= w.min_minute*60 ), nrow(chosen.s.info))
    chosen.lump <- which(diff(lump.idx) >= w.eff_size)
    
    if (length(chosen.lump)  == 0) {
      
      print("There's no LUMP for warming mode.")
    } else {
      
      lump_log <- list()
      
      # examine each lump
      for (l_idx in 1:length(chosen.lump)) {
        
        tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
        tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
        
        lump_log[[length(lump_log)+1]] <- data.frame(chosen.s.info[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx )
      }
      
      lump_log.total <- rbind.fill(lump_log) 
      
      lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))) )
      
      lump_summary <- subset(lump_summary, med.t >= (time.slot.sec - time.slot_margin) & 
                               med.t <= (time.slot.sec + time.slot_margin) )
      
      if (nrow(lump_summary) == 0) {
        
        print("No lump with proper size of time slot.")
      } else {
        
        lump_log.total <- rbind.fill(lump_log[lump_summary$lump_idx])
        
        # insert warming signals (for results)
        new_timestamp <- floor_date(lump_log.total$start.timestamp, unit = "second")
        
        warming.start.idx <- which(answer.log$timestamp %in% new_timestamp)
        warming.on.idx <- unique(unlist(lapply( warming.start.idx, function(x) x + 0:(round(esti_usage.sec)-1) )))
        
        ### avoid overflow of the idx
        warming.on.idx <- warming.on.idx[warming.on.idx <= nrow(answer.log)] 
        answer.log$p[warming.on.idx] <- esti_usage.watt
        
        # print summarized information
        lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
                              min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1), 
                              med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
                              lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
        
        print(lump_summary)
      }
    }
  }
  
  # step) cooking signal computation
  cooking_s.pattern <- subset(s.pattern, sub.delta >= -thres_rp.delta & 
                                sub.delta <= thres_rp.delta )
  
  s.pattern.cooking <- cooking.candidate.search_15Hz(data, start.pattern = cooking_s.pattern, min_watt = c.ap.h1.min, max_watt = c.ap.h1.max, 
                                                     min_t.minute = 0, 
                                                     max_t.minute = c.t.boundary_minute,
                                                     min_fluc.num = c.min.fluc.num)
  
  if (nrow(s.pattern.cooking) == 0){
    
    print("There's no SIGNAL for cooking mode.")
  } else {  
    
    # cooking signal compensation
    cooking.lump_start.idx <- c(0, which(diff(as.numeric(s.pattern.cooking$start.timestamp)) > c.t.boundary_minute*60 )) + 1
    cooking.lump.start.log <- s.pattern.cooking[cooking.lump_start.idx, ]
    
    print("detected cooking mode:")
    cooking_compensated.start.log <- cooking.lump.start.log
    
    print(cooking_compensated.start.log)
    
    if (nrow(cooking_compensated.start.log) != 0) {
      
      # insert cooking signals (for results)
      cooking.start.timestamp <- floor_date(cooking_compensated.start.log$start.timestamp, unit = "second")
      cooking.start.idx <- which(answer.log$timestamp %in% cooking.start.timestamp)

      if( length(cooking.start.idx) > 0 ){
        data('riceCookerRealData')
        cooking.end.idx <- pmin( cooking.start.idx + nrow(riceCookerRealData) - 1, nrow(answer.log) )
        
        for( i in 1:length(cooking.start.idx) ){
          len <- cooking.end.idx[i] - cooking.start.idx[i] + 1
          answer.log$p[cooking.start.idx[i]:cooking.end.idx[i]] <- riceCookerRealData$active_power[1:len]
        }
      }
# 
#       tmp <- mapply( function(x) {tmp_c.on.idx <- x + 0:(-1)
#       tmp_c.on.idx <- tmp_c.on.idx[tmp_c.on.idx <= nrow(answer.log)] 
#       answer.log$p[cooking.on.idx] <- riceCookerRealData$active_power[1:length(cooking.on.idx)]
#       }, cooking.start.idx )
#       
#       cooking.on.idx <- unique(unlist(tmp))
#       
#       ### avoid overflow of the idx
#       cooking.on.idx <- cooking.on.idx[cooking.on.idx <= nrow(answer.log)]
#       answer.log$p[cooking.on.idx] <- c.cooking.watt
      # end
      
    }
  }
  
  return(answer.log)
}

showPatternScan.TV_15Hz <- function (data) {
  
  # exam
  if ( as.numeric(data$timestamp[2]) - as.numeric(data$timestamp[1]) > 0.5 ) {
    stop('Make sure that the input data is from 15Hz signal!')
  }
  
  '%nin%'<- Negate('%in%')
  standby.pwr <- 0.35
  
  # oscillation info
  ap.abs.h1.L <- 17.5
  ap.abs.h1.H <- 32.5
  
  # grouping info
  lump.slot.sec_max <- 180 # 90
  lump.size_min <- 30
  inter.lump.sec_min <- 6
  
  # consumption
  TV.watt <- 240 # 2 options: 1) 256.3519 2) 325.8776 r) 239.7943
  
  # parameters
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"), 
                                            floor_date(max(data$timestamp), unit = "second"), 'secs'))
  answer.log <- data.frame( answer.log, p= rep(standby.pwr,nrow(answer.log)), q= rep(0,nrow(answer.log)))
  
  # examine raw data
  # data <- data[order(data$timestamp),] 
  s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
  
  # examine reverse-magnitude data
  data.rev <- data
  data.rev$active_power <- -data.rev$active_power
  data.rev$reactive_power <- -data.rev$reactive_power
  s.pattern.rev <- DetectPattern_1Hz_new_modifiedforRicecooker(data.rev, position = "start", main_type = "active", sub_type = "reactive")
  s.pattern.rev$h1 <- -s.pattern.rev$h1
  s.pattern.rev$h2 <- -s.pattern.rev$h2
  s.pattern.rev$sub.delta <- -s.pattern.rev$sub.delta
  s.pattern.rev$delta <- -s.pattern.rev$delta
  
  s.pattern.TV <- subset( rbind(s.pattern, s.pattern.rev) , abs(h1) >= ap.abs.h1.L & abs(h1) <= ap.abs.h1.H )
  
  if (nrow(s.pattern.TV) < lump.size_min){
    
    print("There's no SIGNAL for TV.")
  } else {
    
    chosen.s.info <- s.pattern.TV[order(s.pattern.TV$start.idx),] 
    lump.idx <- c(0, which(diff(as.numeric(chosen.s.info$start.timestamp)) >= lump.slot.sec_max ), nrow(chosen.s.info))
    chosen.lump <- which(diff(lump.idx) >= lump.size_min)
    
    if (length(chosen.lump)  == 0) {
      
      print("There's no LUMP for TV.")
    } else {
      
      lump_log <- list()
      
      # examine each lump
      for (l_idx in 1:length(chosen.lump)) {
        
        tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
        tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
        
        lump_log[[length(lump_log)+1]] <- data.frame(chosen.s.info[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx )
      }
      
      lump_log.total <- rbind.fill(lump_log) 
      
      # print summarized information
      lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
                            min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1), 
                            med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
                            lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
      
      
      # insert signals (for results)
      TV.start.timestamp <- floor_date(lump_summary$lump.start, unit = "second")
      TV.start.idx <- which(answer.log$timestamp %in% TV.start.timestamp)
      TV.end.timestamp <- floor_date(lump_summary$lump.end, unit = "second")
      TV.end.idx <- which(answer.log$timestamp %in% TV.end.timestamp)
      TV.on.idx <- unique(unlist(mapply( function(x,y) x:y, TV.start.idx, TV.end.idx )))
      
      ### avoid overflow of the idx
      TV.on.idx <- TV.on.idx[TV.on.idx <= nrow(answer.log)]
      answer.log$p[TV.on.idx] <- TV.watt
      # end
      
      # error compensation
      if (nrow(lump_summary) != 1) {
        
        lump_time.gap <- as.numeric(c( lump_summary$lump.start[-1], max(data$timestamp) ) ) - as.numeric(lump_summary$lump.end) 
        lump_error.compensation <- inter.lump.sec_min* 60 >= lump_time.gap
        lump_summary <- data.frame(lump_summary, compensation =lump_error.compensation, comp.time =  c( lump_summary$lump.start[-1], max(data$timestamp) ) )
      } else {
        
        lump_time.gap <- as.numeric( max(data$timestamp)) - as.numeric(lump_summary$lump.end)  
        lump_error.compensation <- inter.lump.sec_min* 60 >= lump_time.gap
        lump_summary <- data.frame(lump_summary, compensation =lump_error.compensation, comp.time =  max(data$timestamp) )
      }
      
      print(lump_summary)
      compensation.list <- subset(lump_summary, compensation == T)
      
      # insert signals (for results)
      TV.start.timestamp <- floor_date(compensation.list$lump.end, unit = "second")
      TV.start.idx <- which(answer.log$timestamp %in% TV.start.timestamp)
      TV.end.timestamp <- floor_date(compensation.list$comp.time, unit = "second")
      TV.end.idx <- which(answer.log$timestamp %in% TV.end.timestamp)
      TV.on.idx <- unique(unlist(mapply( function(x,y) x:y, TV.start.idx, TV.end.idx )))
      
      ### avoid overflow of the idx
      TV.on.idx <- TV.on.idx[TV.on.idx <= nrow(answer.log)]
      answer.log$p[TV.on.idx] <- TV.watt
      # end
      
    }
  }
  
  return(answer.log)
}

showPatternScan.general <- function (data, appliance) {
  
  if (appliance == '에어컨') {
    
    print('Make sure that the input data is from 15Hz signal!')
    
    # exam
    if ( as.numeric(data$timestamp[2]) - as.numeric(data$timestamp[1]) > 0.5 ) {
      stop('Make sure that the input data is from 15Hz signal!')
    }
    pwr_type <- 'rp'
    pwr_shape <- 'h'
    
    #rising edge
    s.h.low <- 1200 
    s.h.high <- 2000
    s.h2.low <- -Inf
    s.h2.high <- Inf 
    s.sub.low <- -Inf
    s.sub.high <- Inf 
    
    # falling edge
    e.low <- 35 # 61.4955
    e.high <- 118 # 88.69425
    e.sub.low <- 350 # 437
    e.sub.high <- 786 # 686.6827
    
    # time
    time.on_sec <- 450 # 327.843 : from 's.pattern' to 'e.pattern'
    time.flex <- 1 # flexibility w/ prob 0 to 1
    
    # power
    ap.off <- 1
    ap.on <- 498
    
    compensation <- T
    C_factor <- 1
    
    # additional function
    
    
  } else if (appliance == '전자레인지') {
    
    print('Make sure that the input data is from 1Hz signal!')
    
    # exam
    if ( as.numeric(data$timestamp[2]) - as.numeric(data$timestamp[1]) < 0.5 ) {
      stop('Make sure that the input data is from 1Hz signal!')
    }
    
    pwr_type <- 'rp'
    pwr_shape <- 'h'
    
    #rising edge
    s.h.low <- 500 # 538.873
    s.h.high <- 730 # 652.517
    s.h2.low <- -Inf
    s.h2.high <- Inf
    s.sub.low <- -Inf
    s.sub.high <- Inf 
    
    # falling edge
    e.low <- 200 # 266.8974
    e.high <- 460 # 404.9399
    e.sub.low <- 930 # 1039.322
    e.sub.high <- 1260 # 1152.798
    
    # time
    time.on_sec <- 120 # from 's.pattern' to 'e.pattern'
    time.flex <- 1 # flexibility w/ prob 0 to 1
    
    # power
    ap.off <- 0.7
    ap.on <- 1081 # 1) 4 min -- 1113.005, 2) 1min -- 1047.612 3) 2min -- 1118.317 4) 2 min -- 1048.807 
    
    compensation <- T
    C_factor <- 0
    
  } else if (appliance == '청소기') {
    
    print('Make sure that the input data is from 15Hz signal!')
    
    # exam
    if ( as.numeric(data$timestamp[2]) - as.numeric(data$timestamp[1]) > 0.5 ) {
      stop('Make sure that the input data is from 15Hz signal!')
    }
    
    pwr_type <- 'rp'
    pwr_shape <- 'h'
    
    #rising edge
    s.h.low <- 85 # 119.683
    s.h.high <- 185 # 145.417
    s.h2.low <- -165 # -128.919
    s.h2.high <- -55  # -102.793
    s.sub.low <- 130 # 178.445
    s.sub.high <- 260 # 199.797
    
    # falling edge
    e.low <- 5 # around 15
    e.high <- 30 # around 15 
    e.sub.low <- 150 # 176.0965
    e.sub.high <- 210 # 184.8756
    
    # time
    time.on_sec <- 120 # from 's.pattern' to 'e.pattern'
    time.flex <- 1 # flexibility w/ prob 0 to 1
    
    # power
    ap.off <- 0.25
    ap.on <- 179 # 5 min. 30 sec.: 176.6576 
    
    compensation <- T
    C_factor <- 1
    
  } else if (appliance == '전기주전자') {
    
    print('Make sure that the input data is from 15Hz signal!')
    
    # exam
    if ( as.numeric(data$timestamp[2]) - as.numeric(data$timestamp[1]) > 0.5 ) {
      stop('Make sure that the input data is from 15Hz signal!')
    }
    
    pwr_type <- 'ap'
    pwr_shape <- 'h'
    
    #rising edge
    s.h.low <- 450 # 503.665
    s.h.high <- 550 # 503.665
    s.h2.low <- -20 # -0.105
    s.h2.high <- Inf # 0
    s.sub.low <- -15 # 2.737
    s.sub.high <- 18 # 2.737
    
    # falling edge
    e.low <- 450 # 504.7744
    e.high <- 550 # 504.7744
    e.sub.low <- -15 # 2.552750
    e.sub.high <- 20 # 2.552750
    
    # time
    time.on_sec <- 175 # from 's.pattern' to 'e.pattern'
    time.flex <- 1 # set this parameter so 'time.on_sec' is less than 350 secs.
    
    # power
    ap.off <- 0.2
    ap.on <- 502 # 274 sec. : 501.7843
    
    compensation <- F
    C_factor <- 1
    
  } else {
    
    stop('Appliance name is not valid.')
  }
  
  return(showPatternScan.general_main(data = data, pwr_type = pwr_type, pwr_shape = pwr_shape, s.h.low = s.h.low, s.h.high = s.h.high, s.h2.low = s.h2.low, s.h2.high = s.h2.high,
                                      s.sub.low = s.sub.low, s.sub.high = s.sub.high, e.low = e.low, e.high = e.high, e.sub.low = e.sub.low, e.sub.high = e.sub.high, 
                                      time.on_sec = time.on_sec, time.flex = time.flex, ap.off = ap.off, ap.on = ap.on, compensation = compensation, C_factor = C_factor)  )
  
}

showPatternScan.general_main <- function (data, pwr_type, pwr_shape, s.h.low, s.h.high, s.h2.low, s.h2.high, s.sub.low, s.sub.high,
                                          e.low, e.high, e.sub.low, e.sub.high, 
                                          time.on_sec, time.flex, ap.off, ap.on, compensation, C_factor) {
  
  '%nin%'<- Negate('%in%')
  
  # parameters
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"), 
                                            floor_date(max(data$timestamp), unit = "second"), 'secs'))
  answer.log <- data.frame( answer.log, p= rep(ap.off,nrow(answer.log)), q= rep(0,nrow(answer.log)))
  
  # examine raw data
  # data <- data[order(data$timestamp),] 
  
  if (pwr_type == 'ap'){
    
    s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
  } else if (pwr_type == 'rp') {
    
    s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "reactive", sub_type = "active")
  } else {
    
    stop('Power type is not valid.')
  }
  
  s.pattern <- s.pattern[order(s.pattern$start.idx), ]
  
  # step 1) rising edge
  
  if (pwr_shape == 'h'){
    
    chosen.s.pattern <- subset( s.pattern, h1 >= s.h.low & h1 <= s.h.high & h2 >= s.h2.low & 
                                  h2 <= s.h2.high & sub.delta >= s.sub.low & sub.delta <= s.sub.high )
  } else if (pwr_shape == 'd') {
    
  } else {
    
    stop('Power shape is not valid.')
  }
  
  if (nrow(chosen.s.pattern) == 0){
    
    print("There's no SIGNAL for the appliance.")
  } else {
    
    if (pwr_type == 'ap'){
      
      e.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "end", main_type = "active", sub_type = "reactive", c_factor = C_factor)
    } else if (pwr_type == 'rp') {
      
      e.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "end", main_type = "reactive", sub_type = "active", c_factor = C_factor)
    }
    
    
    
    chosen.e.pattern <- subset(e.pattern, h2 >= e.low & h2 <= e.high & sub.delta >= e.sub.low & sub.delta <= e.sub.high )
    
    # edge matching
    #============================================================================================  
    # < compute watt at power-on >
    # matching_info <- mapply(function(s,e) {z <- subset(chosen.e.pattern, start.timestamp %within% interval(s,e) )
    #                                        if (nrow(z) == 1) {
    #                                           z1 <- subset(data, timestamp %within% interval(s,z$end.timestamp) )
    #                                           return(sum(z1$active_power) / nrow(z1))
    #                                        } }, chosen.s.pattern$start.timestamp,
    #                                        c(chosen.s.pattern$start.timestamp[-1], max(data$timestamp) ), SIMPLIFY = F)
    #============================================================================================  
    
    matching_info <- mapply(function(s,e,idx) {   z <- subset(chosen.e.pattern, start.timestamp %within% interval(s,e) )
    if (nrow(z) >= 1) {
      
      tmp_time.on <- as.numeric(z$end.timestamp) - as.numeric(s)
      tmp_valid.f.edge <- rep(F, length(tmp_time.on))
      chosen.f.edge.idx <- which.min(abs(time.on_sec - tmp_time.on) )
      if ( tmp_time.on[chosen.f.edge.idx] >= time.on_sec*(1-time.flex) && 
           tmp_time.on[chosen.f.edge.idx] <= time.on_sec*(1+time.flex) ) {
        
        tmp_valid.f.edge[chosen.f.edge.idx] <- T
      } 
      return(data.frame(z, r.edge.timestamp = s, time.on = tmp_time.on, valid.f.edge = tmp_valid.f.edge, r.edge.idx = idx) )
      
    } },
    chosen.s.pattern$start.timestamp, c(chosen.s.pattern$start.timestamp[-1], max(data$timestamp) ), seq(1,nrow(chosen.s.pattern)),  
    SIMPLIFY = F)
    
    matching_info <- rbind.fill(matching_info)
    print(matching_info)
    
    
    
    if (!is.null(matching_info)){
      
      matched.set <- subset(matching_info, valid.f.edge == T)
      redundant.r.edge.idx <- which(seq(1,nrow(chosen.s.pattern)) %nin% matched.set$r.edge.idx)
      redundant.f.edge <- subset(matching_info, valid.f.edge == F)
      
    } else {
      matched.set <- data.frame()
      if (nrow(chosen.s.pattern) != 0) {
        redundant.r.edge.idx <- seq(1,nrow(chosen.s.pattern)) 
      } else {
        redundant.r.edge.idx <- c()
      }
      redundant.f.edge <- data.frame()
    }
    
    
    # result 1) matched set
    
    # insert signals (for results)
    if (nrow(matched.set) != 0) {
      
      r.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp[matched.set$r.edge.idx], unit = "second")
      r.start.idx <- which(answer.log$timestamp %in% r.start.timestamp)
      r.end.timestamp <- floor_date(matched.set$end.timestamp, unit = "second")
      r.end.idx <- which(answer.log$timestamp %in% r.end.timestamp)
      r.on.idx <- unique(unlist(mapply( function(x,y) x:y, r.start.idx, r.end.idx )))
      
      ### avoid overflow of the idx
      r.on.idx <- r.on.idx[r.on.idx <= nrow(answer.log)]
      answer.log$p[r.on.idx] <- ap.on
      # end
      
    }
    
    
    # result 2) redundant r.edge set
    
    if (compensation == T) {
      
      if (length(redundant.r.edge.idx) != 0) {
        
        # insert signals (for results)
        rr.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp[redundant.r.edge.idx], unit = "second")
        rr.start.idx <- which(answer.log$timestamp %in% rr.start.timestamp)
        
        # all rising edges are considered to set appropriate time slots for consideration
        all.r.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp, unit = "second")
        all.r.start.idx <- which(answer.log$timestamp %in% all.r.start.timestamp)
        
        rr.end.limit <- unique( c(all.r.start.idx , nrow(answer.log)) ) # compute each minimum time slot
        rr.end.idx <- mapply(function(x) min(rr.end.limit[rr.end.limit > x]), rr.start.idx)
        r.on.idx <- unique(unlist(mapply( function(x,y) x:min((x+time.on_sec-1), y), rr.start.idx, rr.end.idx )))
        
        ### avoid overflow of the idx
        r.on.idx <- r.on.idx[r.on.idx <= nrow(answer.log)]
        answer.log$p[r.on.idx] <- ap.on
        # end
        
      }
      
    }
    
    
    
    # result 3) redundant f.edge set
    
  }
  
  return(answer.log)
  
}


# showPatternScan.general_15Hz <- function (data, appliance) {
# 
#   if (appliance == '에어컨') {
#     
#     pwr_type <- 'rp'
#     pwr_shape <- 'h'
#     
#     #rising edge
#     s.h.low <- 1200 
#     s.h.high <- 2000
#     
#     # falling edge
#     e.low <- 50 # 80.4955
#     e.high <- 118 # 88.69425
#     e.sub.low <- 546 # 646.7661
#     e.sub.high <- 786 # 686.6827
#     
#     # time
#     time.on_sec <- 328 # 327.843 : from 's.pattern' to 'e.pattern'
#     time.flex <- 0.5 # flexibility w/ prob 0 to 1
#     
#     # power
#     ap.off <- 1
#     ap.on <- 657.5
#     
#     # additional function
#     
#     
#   } else if (appliance == '전자레인지') {
#     
#     pwr_type <- 'rp'
#     pwr_shape <- 'h'
#     
#     # power
#     ap.off <- 1
#     
#   } else {
#     
#     stop('Appliance name is not valid.')
#   }
# 
#   '%nin%'<- Negate('%in%')
#   
#   # parameters
#   str.t <- as.character(min(data$timestamp))
#   end.t <- as.character(max(data$timestamp))
#   answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"), 
#                                             floor_date(max(data$timestamp), unit = "second"), 'secs'))
#   answer.log <- data.frame( answer.log, p= rep(ap.off,nrow(answer.log)), q= rep(0,nrow(answer.log)))
#   
#   # examine raw data
#   data <- data[order(data$timestamp),] 
#   
#   if (pwr_type == 'ap'){
#     
#     s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
#   } else if (pwr_type == 'rp') {
#     
#     s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "reactive", sub_type = "active")
#   } else {
#     
#     stop('Power type is not valid.')
#   }
#   
#   s.pattern <- s.pattern[order(s.pattern$start.idx), ]
#   
#   # step 1) rising edge
#   
#   if (pwr_shape == 'h'){
#     
#     chosen.s.pattern <- subset( s.pattern, h1 >= s.h.low & h1 <= s.h.high)
#   } else if (pwr_shape == 'd') {
#     
#   } else {
#     stop('Power shape is not valid.')
#   }
#   
#   if (nrow(chosen.s.pattern) == 0){
#     
#     stop("There's no SIGNAL for the appliance.")
#   } else {
#     
#     e.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "end", main_type = "reactive", sub_type = "active", c_factor = 1)
#     
#     chosen.e.pattern <- subset(e.pattern, h2 >= e.low & h2 <= e.high & sub.delta >= e.sub.low & sub.delta <= e.sub.high )
#     
#     # edge matching
#     #============================================================================================  
#     # < compute watt at power-on >
#     # matching_info <- mapply(function(s,e) {z <- subset(chosen.e.pattern, start.timestamp %within% interval(s,e) )
#     #                                        if (nrow(z) == 1) {
#     #                                           z1 <- subset(data, timestamp %within% interval(s,z$end.timestamp) )
#     #                                           return(sum(z1$active_power) / nrow(z1))
#     #                                        } }, chosen.s.pattern$start.timestamp,
#     #                                        c(chosen.s.pattern$start.timestamp[-1], max(data$timestamp) ), SIMPLIFY = F)
#     #============================================================================================  
#     
#     matching_info <- mapply(function(s,e,idx) {   z <- subset(chosen.e.pattern, start.timestamp %within% interval(s,e) )
#                                               if (nrow(z) >= 1) {
#                                              
#                                                 tmp_time.on <- as.numeric(z$end.timestamp) - as.numeric(s)
#                                                 tmp_valid.f.edge <- rep(F, length(tmp_time.on))
#                                                 chosen.f.edge.idx <- which.min(abs(time.on_sec - tmp_time.on) )
#                                                 if ( tmp_time.on[chosen.f.edge.idx] >= time.on_sec*(1-time.flex) && 
#                                                      tmp_time.on[chosen.f.edge.idx] <= time.on_sec*(1+time.flex) ) {
#                                                   
#                                                   tmp_valid.f.edge[chosen.f.edge.idx] <- T
#                                                 } 
#                                                 return(data.frame(z, time.on = tmp_time.on, valid.f.edge = tmp_valid.f.edge, r.edge.idx = idx) )
#                                                                    
#                                               } },
#                                           chosen.s.pattern$start.timestamp, c(chosen.s.pattern$start.timestamp[-1], max(data$timestamp) ), seq(1,nrow(chosen.s.pattern)),  
#                                           SIMPLIFY = F)
#     
#     matching_info <- rbind.fill(matching_info)
#     print(matching_info)
#     
#     matched.set <- subset(matching_info, valid.f.edge == T)
#     redundant.r.edge.idx <- which(seq(1,nrow(chosen.s.pattern)) %nin% matched.set$r.edge.idx)
#     redundant.f.edge <- subset(matching_info, valid.f.edge == F)
#     
#     
#     # result 1) matched set
# 
#     if (nrow(matched.set) != 0) {
#       
#       # insert signals (for results)
#       r.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp[matched.set$r.edge.idx], unit = "second")
#       r.start.idx <- which(answer.log$timestamp %in% r.start.timestamp)
#       r.end.timestamp <- floor_date(matched.set$end.timestamp, unit = "second")
#       r.end.idx <- which(answer.log$timestamp %in% r.end.timestamp)
#       r.on.idx <- unique(unlist(mapply( function(x,y) x:y, r.start.idx, r.end.idx )))
#       
#       ### avoid overflow of the idx
#       r.on.idx <- r.on.idx[r.on.idx <= nrow(answer.log)]
#       answer.log$p[r.on.idx] <- ap.on
#       # end
#     }
# 
# 
#     
#     # result 2) redundant r.edge set
#     
#     if (length(redundant.r.edge.idx) != 0) {
#       
#       # insert signals (for results)
#       rr.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp[redundant.r.edge.idx], unit = "second")
#       rr.start.idx <- which(answer.log$timestamp %in% rr.start.timestamp)
#       
#       # all rising edges are considered to set appropriate time slots for consideration
#       all.r.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp, unit = "second")
#       all.r.start.idx <- which(answer.log$timestamp %in% all.r.start.timestamp)
#       
#       rr.end.limit <- unique( c(all.r.start.idx , nrow(answer.log)) ) # compute each minimum time slot
#       rr.end.idx <- mapply(function(x) min(rr.end.limit[rr.end.limit > x]), rr.start.idx)
#       r.on.idx <- unique(unlist(mapply( function(x,y) x:min((x+time.on_sec-1), y), rr.start.idx, rr.end.idx )))
#       
#       ### avoid overflow of the idx
#       r.on.idx <- r.on.idx[r.on.idx <= nrow(answer.log)]
#       answer.log$p[r.on.idx] <- ap.on
#       # end
#       
#     }
# 
#     
#     # result 3) redundant f.edge set
#     
#   }
#   
#   return(answer.log)
# 
# }
# 
# showPatternScan.general_1Hz <- function (data, appliance) {
#   
#   if (appliance == '전자레인지') {
#     
#     print('Make sure that the input data is from 1Hz signal !!!')
#     pwr_type <- 'rp'
#     pwr_shape <- 'h'
#     
#     #rising edge
#     s.h.low <- 500 # 538.873
#     s.h.high <- 730 # 652.517
#     
#     # falling edge
#     e.low <- 200 # 266.8974
#     e.high <- 460 # 404.9399
#     e.sub.low <- 930 # 1039.322
#     e.sub.high <- 1260 # 1152.798
#     
#     # time
#     time.on_sec <- 300 # from 's.pattern' to 'e.pattern'
#     time.flex <- 1 # flexibility w/ prob 0 to 1
#     
#     # power
#     ap.off <- 0.7
#     ap.on <- 1081 # 1) 4 min -- 1113.005, 2) 1min -- 1047.612 3) 2min -- 1118.317 4) 2 min -- 1048.807 
#     
#     
#   } else {
#     
#     stop('Appliance name is not valid.')
#   }
#   
#   '%nin%'<- Negate('%in%')
#   
#   # parameters
#   str.t <- as.character(min(data$timestamp))
#   end.t <- as.character(max(data$timestamp))
#   answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"), 
#                                             floor_date(max(data$timestamp), unit = "second"), 'secs'))
#   answer.log <- data.frame( answer.log, p= rep(ap.off,nrow(answer.log)), q= rep(0,nrow(answer.log)))
#   
#   # examine raw data
#   data <- data[order(data$timestamp),] 
#   
#   if (pwr_type == 'ap'){
#     
#     s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
#   } else if (pwr_type == 'rp') {
#     
#     s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "reactive", sub_type = "active")
#   } else {
#     
#     stop('Power type is not valid.')
#   }
#   
#   s.pattern <- s.pattern[order(s.pattern$start.idx), ]
#   
#   # step 1) rising edge
#   
#   if (pwr_shape == 'h'){
#     
#     chosen.s.pattern <- subset( s.pattern, h1 >= s.h.low & h1 <= s.h.high)
#   } else if (pwr_shape == 'd') {
#     
#   } else {
#     stop('Power shape is not valid.')
#   }
#   
#   if (nrow(chosen.s.pattern) == 0){
#     
#     stop("There's no SIGNAL for the appliance.")
#   } else {
#     
#     e.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "end", main_type = "reactive", sub_type = "active", c_factor = 1)
#     
#     chosen.e.pattern <- subset(e.pattern, h2 >= e.low & h2 <= e.high & sub.delta >= e.sub.low & sub.delta <= e.sub.high )
#     
#     # edge matching
#     #============================================================================================  
#     # < compute watt at power-on >
#     # matching_info <- mapply(function(s,e) {z <- subset(chosen.e.pattern, start.timestamp %within% interval(s,e) )
#     #                                        if (nrow(z) == 1) {
#     #                                           z1 <- subset(data, timestamp %within% interval(s,z$end.timestamp) )
#     #                                           return(sum(z1$active_power) / nrow(z1))
#     #                                        } }, chosen.s.pattern$start.timestamp,
#     #                                        c(chosen.s.pattern$start.timestamp[-1], max(data$timestamp) ), SIMPLIFY = F)
#     #============================================================================================  
#     
#     matching_info <- mapply(function(s,e,idx) {   z <- subset(chosen.e.pattern, start.timestamp %within% interval(s,e) )
#     if (nrow(z) >= 1) {
#       
#       tmp_time.on <- as.numeric(z$end.timestamp) - as.numeric(s)
#       tmp_valid.f.edge <- rep(F, length(tmp_time.on))
#       chosen.f.edge.idx <- which.min(abs(time.on_sec - tmp_time.on) )
#       if ( tmp_time.on[chosen.f.edge.idx] >= time.on_sec*(1-time.flex) && 
#            tmp_time.on[chosen.f.edge.idx] <= time.on_sec*(1+time.flex) ) {
#         
#         tmp_valid.f.edge[chosen.f.edge.idx] <- T
#       } 
#       return(data.frame(z, r.edge.timestamp = s, time.on = tmp_time.on, valid.f.edge = tmp_valid.f.edge, r.edge.idx = idx) )
#       
#     } },
#     chosen.s.pattern$start.timestamp, c(chosen.s.pattern$start.timestamp[-1], max(data$timestamp) ), seq(1,nrow(chosen.s.pattern)),  
#     SIMPLIFY = F)
#     
#     matching_info <- rbind.fill(matching_info)
#     print(matching_info)
#     
#     matched.set <- subset(matching_info, valid.f.edge == T)
#     redundant.r.edge.idx <- which(seq(1,nrow(chosen.s.pattern)) %nin% matched.set$r.edge.idx)
#     redundant.f.edge <- subset(matching_info, valid.f.edge == F)
#     
#     
#     # result 1) matched set
#     
#     # insert signals (for results)
#     if (nrow(matched.set) != 0) {
#       
#       r.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp[matched.set$r.edge.idx], unit = "second")
#       r.start.idx <- which(answer.log$timestamp %in% r.start.timestamp)
#       r.end.timestamp <- floor_date(matched.set$end.timestamp, unit = "second")
#       r.end.idx <- which(answer.log$timestamp %in% r.end.timestamp)
#       r.on.idx <- unique(unlist(mapply( function(x,y) x:y, r.start.idx, r.end.idx )))
#       
#       ### avoid overflow of the idx
#       r.on.idx <- r.on.idx[r.on.idx <= nrow(answer.log)]
#       answer.log$p[r.on.idx] <- ap.on
#       # end
#       
#     }
# 
#     # result 2) redundant r.edge set
#     
#     if (length(redundant.r.edge.idx) != 0) {
#       
#       # insert signals (for results)
#       rr.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp[redundant.r.edge.idx], unit = "second")
#       rr.start.idx <- which(answer.log$timestamp %in% rr.start.timestamp)
#       
#       # all rising edges are considered to set appropriate time slots for consideration
#       all.r.start.timestamp <- floor_date(chosen.s.pattern$start.timestamp, unit = "second")
#       all.r.start.idx <- which(answer.log$timestamp %in% all.r.start.timestamp)
#       
#       rr.end.limit <- unique( c(all.r.start.idx , nrow(answer.log)) ) # compute each minimum time slot
#       rr.end.idx <- mapply(function(x) min(rr.end.limit[rr.end.limit > x]), rr.start.idx)
#       r.on.idx <- unique(unlist(mapply( function(x,y) x:min((x+time.on_sec-1), y), rr.start.idx, rr.end.idx )))
#       
#       ### avoid overflow of the idx
#       r.on.idx <- r.on.idx[r.on.idx <= nrow(answer.log)]
#       answer.log$p[r.on.idx] <- ap.on
#       # end
#       
#     }
# 
#     # result 3) redundant f.edge set
#     
#   }
#   
#   return(answer.log)
#   
# }

#============================================================================================  
# < supplementary functions >
#============================================================================================

