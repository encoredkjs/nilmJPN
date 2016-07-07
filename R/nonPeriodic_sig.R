###############################################################################################
# Functions for JPN
###############################################################################################
generate.PatternScan.meta.ricecooker_JPN <- function (data, Hz, eff_size = 5, periodicity = 0.2){
  
  # options(scipen = 15)
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"),
                                            floor_date(max(data$timestamp), unit = "second"), 'secs'))
  answer.log <- data.frame( answer.log, p= rep(0,nrow(answer.log)), q= rep(0,nrow(answer.log)))
  
  # grouping parameters
  search_iter <- 650 # increased by 50 due to the JPN site '10012085'
  quantileProb <- 0.25
  
  e_pattern <- DetectPattern_1Hz_new(data, position = "end", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
  
  # major signal search
  major_min_watt <- 250
  major_max_watt <- 1600
  major_medSec_min <- 10
  major_medSec_max <- 600
  major_lumpGap <- 3600
  
     # DB flag - 0: cooking || 1: warming || 2: both
  DB_major <- data.frame(time = c(15, 16, 20, 32, 45, 50, 60, 75, 128), 
                         flag = c( 0,  2,  2,  1,  1,  1,  2,  1,   1))
  DB_toleranceSec <- 0.4

  majorSig_info <- parallel_majorSig_search( patterns = e_pattern, min_watt = major_min_watt, max_watt = major_max_watt, group_periodicity = periodicity, group_quantProb = quantileProb, 
                                             group_searchIter = search_iter, group_medSec_min = major_medSec_min, group_medSec_max = major_medSec_max, effective_size = eff_size, 
                                             lump_timeGapThres = major_lumpGap, timeTolerance = DB_toleranceSec, DB_majorSig = DB_major)
  
  if( length(majorSig_info) == 0) {
    print('Detection has failed - no info for the major signal')
    # result frame

    
    return(answer.log) # return(list())
  }
  
  thres_timeDuration <- 3600 * 2
  proper_cookTime <- 40 * 60
  min_cookOn <- 5 * 60
  ap_ignoreThres <- 5 
  
  validInfo <- sigOrchestration_riceCooker(data, Hz, answer.log, major_SigInfo = majorSig_info, major_DBInfo = DB_major, thresTime = thres_timeDuration, 
                                           range_energyEsti = proper_cookTime, min_usageTime = min_cookOn, min_apThres = ap_ignoreThres)
  
  return(validInfo)


  

  
  
  
    cooking.max_t.sec <- 40*60
    cooking.min_fluc.num <- 3
    
    CDB_tolerance <- .02 # 2%
    
    c.lump_eff_size <- 4

  # tmp_param
  eff_group_size = 40
  thres_ap.h = 15
  thres_rp.delta = 20
  
  # warming mode search
  w.thres_ap.max <- 170
  
  # time parameter
  if (nation == 'KOR') w.med_sec.min <- 10 else if (nation == 'JPN') w.med_sec.min <- 2 
  w.med_sec.max <- 240
  
  # w.lump (basic parameters)
  w.lump.gap_sec.max <- 270
  w.lump.med_sec.margin <- 1  
  
  # w.lump (candidate DB)
  if (nation == 'KOR') {
    w.cand.DB <- data.frame(med.t = c(16, 30, 32, 48, 64, 80, 96, 112, 128, 144, 160), #30 is from the brand 'Coochan'
                            min.med.rate2 = c(1, 1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10 ) )
    # magnitude DB (off the record): 1(106, 70, (91.5, 103.5), 42)
    # magnitude DB (off the record): 2(40.5), 3(100), 4(60)
    
  } else if (nation == 'JPN') {
    w.cand.DB <- data.frame(med.t = c(4, 16, 32, 48), 
                            min.med.rate2 = c(1, 1, 1/2, 1/3) )
    w.cand.extDB <- data.frame(med.t = c(25, 50, 50, 75, 75, 75, 100, 100, 100), 
                               min.med.rate2 = c(1, 1, 1/2, 2/3, 1/3, 1/5, 1, 3/4, 1/2) )
  }
  
  
  DB_tolerance <- .02 # 2%
  DB_rate_quantile <- .08 # 8%
  
  DB_timespan <- .1 # compared to maximum duration
  
  
  # cooking mode search
  if (nation == 'KOR') cooking.min_watt <- 750 else if (nation == 'JPN') cooking.min_watt <- 250
  cooking.max_watt <- 1600
  cooking.min_t.sec <- 7*60
  cooking.max_t.sec <- 40*60
  cooking.min_fluc.num <- 0
  
  # c.lump (basic parameters)
  c.lump.gap_sec.max <- 3600
  c.lump.med_sec.margin <- 0.0625 # 1/16 : In case of 16 sec, the margin is 1 sec.
  
  # cooking repetition DB
  if (nation == 'KOR') {
    c.cand.DB <- data.frame(med.t = c(16, 32))
  } else if (nation == 'JPN') {
    c.cand.DB <- data.frame(med.t = c(15, 16, 60))
  }
  
  CDB_tolerance <- .02 # 2%
  
  # connection between warming & cooking
  wnc.forward.time_diff.max_hr <- 4
  wnc.backward.time_diff.max_hr <- 0.5
  lump_merge.max_hr <- 1/3
  
  
  # cooking mode refine
  margin_sub.delta <- 20
  margin_delta <- 50 # watt
  margin_sec <- 120*60 #second
  
  # cooking.med_watt <- 1200
  
  # warming & cooking signal search (in parallel)
  w.pattern.s <- DetectPattern_1Hz_new(data, position = "start", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
  
  warm.info_low <- parallel_warmSig_search_15Hz(data, w.pattern.s, ap_h1_min = thres_ap.h, ap_h1_max = w.thres_ap.max, ap_delta_min = thres_ap.h, ap_delta_max = w.thres_ap.max,
                                                rp_delta_thres = thres_rp.delta, main.g.name = "h1", sub.g.name = "delta", periodicity, search_iter, w.med_sec.min, w.med_sec.max,
                                                eff_size, eff_group_size, w.lump.gap_sec.max, w.lump.med_sec.margin, w.cand.DB, DB_tolerance, DB_rate_quantile, DB_timespan)
  # CAUTION!!!
  # warm.info_low <- list()
  
  warm.info_high <- list()
  if (!any(unlist(warm.info_low$type) > 0)){
    
    if (nation == 'KOR'){
      warm.info_high <- parallel_warmSig_search_15Hz(data, w.pattern.s, ap_h1_min = w.thres_ap.max, ap_h1_max = cooking.min_watt, ap_delta_min = -Inf, ap_delta_max = Inf,
                                                     rp_delta_thres = Inf, main.g.name = "h1", sub.g.name = "sub.delta", periodicity, search_iter, w.med_sec.min, w.med_sec.max, 
                                                     eff_size, eff_group_size, w.lump.gap_sec.max, w.lump.med_sec.margin, w.cand.DB, DB_tolerance, DB_rate_quantile, DB_timespan)
      
    } else if (nation == 'JPN'){
      JPN_w.thres_ap.extMax <- 1000
      w.cand.DB <- w.cand.extDB
      warm.info_high <- parallel_warmSig_search_15Hz(data, w.pattern.s, ap_h1_min = w.thres_ap.max, ap_h1_max = JPN_w.thres_ap.extMax, ap_delta_min = -Inf, ap_delta_max = Inf,
                                                     rp_delta_thres = Inf, main.g.name = "h1", sub.g.name = "sub.delta", periodicity, search_iter, w.med_sec.min, w.med_sec.max, 
                                                     eff_size, eff_group_size, w.lump.gap_sec.max, w.lump.med_sec.margin, w.cand.DB, DB_tolerance, DB_rate_quantile, DB_timespan)
    }
    
  }
  
  cook.info <- parallel_cookSig_search_15Hz(data, cooking.min_watt, cooking.max_watt, cooking.max_t.sec, cooking.min_fluc.num, periodicity,
                                            search_iter, c.med_sec.min = w.med_sec.min, c.med_sec.max = w.med_sec.max, eff_size, c.lump.gap_sec.max,
                                            c.lump.med_sec.margin, c.cand.DB, CDB_tolerance)
  
  valid.info <- sig.Orchestration_RC_15Hz(data, warm.info_low, warm.info_high, cook.info, w.cand.DB, c.cand.DB, DB_timespan,
                                          cooking.forward_search = wnc.forward.time_diff.max_hr, cooking.backward_search = wnc.backward.time_diff.max_hr,
                                          lump_merge = lump_merge.max_hr, cooking.min_sec = cooking.min_t.sec,
                                          margin_sub.delta, margin_delta, margin_sec)
  
  if (length(valid.info) == 0) {
    print("Valid signal search for a rice cooker has failed.")
    return(list())
  }
  
  # make parameter group to return
  parameters <- c('eff_size' = eff_size, 'thres_ap.h' = thres_ap.h, 'thres_rp.delta' = thres_rp.delta, 'w.lump.gap_sec.max' = w.lump.gap_sec.max,
                  'w.lump.med_sec.margin' = w.lump.med_sec.margin, 'c.max_t.sec' = cooking.max_t.sec, 'c.lump.gap_sec.max' = c.lump.gap_sec.max,
                  'c.lump.med_sec.margin' = c.lump.med_sec.margin, 'c.lump_merge.max_hr' = lump_merge.max_hr, 'r.margin_sub.delta' = margin_sub.delta,
                  'r.margin_delta' = margin_delta, 'r.margin_sec' = margin_sec)
  
  # JSON result
  return(JSON_result_15Hz(data, information = valid.info, ap_thres = thres_ap.h, min_cook.time = cooking.min_t.sec, parameters, 
                          warm_DB_time = w.cand.DB$med.t))
  
}

parallel_majorSig_search <- function(patterns, min_watt, max_watt, group_periodicity, group_quantProb, group_searchIter, group_medSec_min,
                                     group_medSec_max, effective_size, lump_timeGapThres, timeTolerance, DB_majorSig){
  
  # search periodicity in signature level
  major_pattern <- patterns %>% filter(h2 >= -max_watt) %>% filter(h2 <= -min_watt)
  if( nrow(major_pattern) == 0 ){
    print("There is no candidate for the major signals")
    return(list())
  } 

  major_pattern$delta <- -major_pattern$delta

  major_group_info <- detectGroup_riceCooker(data = major_pattern, resolution = group_searchIter, main.g.name = "delta", sub.g.name = "sub.delta", p_factor = group_periodicity,
                                             q_factor = group_quantProb, med.time_min = group_medSec_min, med.time_max = group_medSec_max, group_size = effective_size)

  if (length(major_group_info) == 0) {
    print("Group search to retrieve the major DB has failed.")
    return(list())
  }

  print("major groups:")
  group_list <- cbind(major_group_info[[length(major_group_info)]], data.frame(orderNum = seq( length(major_group_info)-1 ) )) %>% arrange(desc(sum))
  print(group_list)
  
  # examine candidates 1) lapply for each group; 2) lapply for each DB time
  # for ( candIdx in 1:nrow(group_list) ) { print(candIdx)
  candResult <- lapply(seq(nrow(group_list)), function(candIdx){
    chosen_group <- major_group_info[[ group_list$orderNum[candIdx] ]] %>% arrange(start.timestamp)
    time_diff <- difftime( tail(chosen_group$start.timestamp,-1), head(chosen_group$start.timestamp,-1), units='secs')
    group_detail <- lapply(DB_majorSig$time, function(DB_time){ 
      
                                             chosen_rows <- which( ( (DB_time-timeTolerance) <= time_diff) & ((DB_time+timeTolerance) >= time_diff))
                                             chosen_rows <- unique( c(chosen_rows, chosen_rows+1) ) %>% sort() %>% chosen_group[.,]
                                             if(nrow(chosen_rows) != 0) {
                                               
                                               # split the chosen rows into several lumps depending on time gap
                                               splitInterval <- cumsum( c(TRUE, difftime( tail(chosen_rows$start.timestamp,-1),
                                                                                          head(chosen_rows$start.timestamp,-1), units='secs') >= lump_timeGapThres))
                                               splitLumps <- bind_cols(chosen_rows, data.frame(lumpIdx = splitInterval)) 
                                               # splitLumps <- split.data.frame( chosen_rows, splitInterval) %>% bind_rows(.id = 'lumpIdx')
                                               
                                               # summarize information
                                               lumpSummary <- ddply(splitLumps, .(lumpIdx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
                                                                    min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), max.t = max(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
                                                                    min.d = min(delta), med.d = median(delta), max.d = max(delta), sd.d = sd(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
                                                                    lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
                                               
                                               lumpSummary <- lumpSummary %>% filter(sum >= effective_size)
                                               if(nrow(lumpSummary) != 0){
                                                 lumpLog <- splitLumps %>% filter( lumpIdx %in% lumpSummary$lumpIdx)
                                                 return(list(lumpSummary, lumpLog, DB_time, candIdx))
                                               } else return(NULL)
                                             }
                                             return(NULL) })
    return(group_detail[!sapply(group_detail, is.null)])
  })
  
  
      
  # return results
  ID_type <- list()
  ID_total.log <- list()
  ID_log.summary <- list()
  ID <- list()
  for(f_idx in seq(length(candResult))) {
    
    if(length(candResult[[f_idx]]) != 0){
      for(s_idx in seq(length(candResult[[f_idx]]))){
        
        ID_log.summary[[length(ID_log.summary)+1]] <- candResult[[f_idx]][[s_idx]][[1]]
        ID_total.log[[length(ID_total.log)+1]] <- candResult[[f_idx]][[s_idx]][[2]]
        ID_type[[length(ID_type)+1]] <- candResult[[f_idx]][[s_idx]][[3]]
        ID[[length(ID)+1]] <- candResult[[f_idx]][[s_idx]][[4]]
        
      }
    }
  }
  
  # return entire pattern information
  entirePatternLog <- list()
  for(f_idx in seq(nrow(group_list))) entirePatternLog[[length(entirePatternLog)+1]] <- major_group_info[[ group_list$orderNum[f_idx] ]] %>% arrange(start.timestamp)

  return(list( wholePattern = entirePatternLog, candGroup = ID, candType = ID_type, candPattern = ID_total.log, candSummary = ID_log.summary ))
}

detectGroup_riceCooker <- function(data, resolution, main.g.name, sub.g.name, p_factor, q_factor = 0.05, med.time_min, med.time_max, group_size){
  
  findGaps <- function(x,n){
    x <- sort(x)
    x.diff <- data.frame( val = diff(x), idx = 1:(length(x)-1) )
    wall.idx <- x.diff$idx[ order( x.diff$val, decreasing=T ) ][1:n] # choose first N gaps
    result <- sapply( wall.idx, function(i) mean(x[i+c(0,1)]))
    return( sort(result) )
  }
  
  summarize_timestamp <- function( timestamp ){
    tsDiff <- diff(as.numeric(sort(timestamp)))
    data.frame( 
      'sum' = length(timestamp), 
      'min.t'  = min(tsDiff),
      'min.t2' = quantile(tsDiff, q_factor),
      'med.t'  = median(tsDiff), 
      'lost.sig.num' = sum( pmax( round(tsDiff / median(tsDiff)) -1 ,0)))
  }
  
  logSummary <- list()
  logDetail <- list()
  o.list <- data
  
  r_idx <- 1
  iter_idx <- 1
  while( iter_idx <= resolution ){
    
    if ( nrow(o.list) +1 <= r_idx ){
      print("stop resolution increment for searcing groups")
      break
    }
    
    xValue    <- o.list[,main.g.name] %>% unlist(use.names = FALSE)
    yValue    <- o.list[, sub.g.name] %>% unlist(use.names = FALSE)
    o.list$xDivision <- findInterval( xValue, findGaps( xValue, r_idx ) )
    o.list$yDivision <- findInterval( yValue, findGaps( yValue, r_idx ) )
    
    removeData <- o.list %>% 
      group_by( xDivision, yDivision ) %>%
      filter( n() < group_size )
    
    o.list <- o.list %>% 
      group_by( xDivision, yDivision ) %>%
      filter( n() >= group_size )
    
    r_idx <- pmax( r_idx - length(group_size(removeData)), 1 )
    
    if (nrow(o.list) > 0){ # in case 'o.list' is empty
      
      group.info <- o.list %>%
        do( summarize_timestamp(.$start.timestamp) ) %>%
        mutate( min.med.rate  = min.t/med.t, 
                min.med.rate2 = min.t2/med.t )
      
      group.info <- merge( group.info, o.list %>% 
                             summarise_each_( funs(median), c('h1','delta','sub.delta') )) %>%
        dplyr::rename( med.h1 = h1, med.d = delta, med.sub.d = sub.delta ) %>%
        mutate( lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
    }
    
    ### effective group (conservative search by default)
    if (nrow(group.info) > 0){
      
      eff_group.info <- group.info %>%
        filter( med.t >= med.time_min ) %>%
        filter( med.t <= med.time_max ) %>%
        filter( min.med.rate2 >= p_factor ) %>%
        filter( xDivision != 0 )
      
      if( nrow(eff_group.info) > 0 ){
        
        logSummary[[length(logSummary)+1]] <- eff_group.info
        
        for(g_idx in 1:nrow(eff_group.info) ){
          logDetail[[length(logDetail) +1]] <- o.list %>%
            filter( xDivision == eff_group.info$xDivision[g_idx] & yDivision == eff_group.info$yDivision[g_idx] ) 
          
          o.list <- o.list %>%
            filter( !(xDivision == eff_group.info$xDivision[g_idx] & yDivision == eff_group.info$yDivision[g_idx]) ) 
        }
      }
      
      # index increment for while loop
      r_idx <- pmax( r_idx - nrow(eff_group.info) + 1, 1 )
    } else {
      r_idx <- r_idx + 1
    }
    iter_idx <- iter_idx + 1
    
  }
  
  # return information
  if (length(logSummary) != 0) {
    logDetail[[length(logDetail) +1]] <- rbind.fill(logSummary)
    return(logDetail)
  } else {
    return(list())
  }
  
}

sigOrchestration_riceCooker <- function(data, Hz, answer.log, major_SigInfo, major_DBInfo, thresTime, range_energyEsti, min_usageTime, min_apThres){
  
  # possible properties of signals for the orchestration
  # 1) flag in DB 1) time duration # 1) active power # 1) sum of signatures # 1) number of lumps
  
  # step 1] classify the signals: short(0) or long(1)
  sigClass <- lapply(major_SigInfo$candSummary, function(summaryInfo){
    if( max(summaryInfo$lump.duration) >= thresTime) return(1) else return(0)
  })
  
  # step 2] process short(0) signals to identify the cooking mode
  ID_cookSig <- 0
  # cookEnergy <- data.frame(startTime = 0, endTime = 0,  SecON = 0, whON = 0, apON = 0, minSigNum = 0) 
  if(!any(sigClass == 0)) {
    print('no cooking candidates in major signals')

  } else {
    print('using proper short signals for cooking candidates')
    shortCand <- which(sigClass == 0)
    
    # details of the energy consumption
    energyEsti <- lapply(shortCand, function(sigIdx){
      lumpSummary <- major_SigInfo$candSummary[[sigIdx]]
      lumpLog <- major_SigInfo$candPattern[[sigIdx]] %>% mutate( ap.thres = (data$active_power[start.idx] + data$active_power[end.idx]) / 2 )
      
      estiOutput <- ldply(seq(nrow(lumpSummary)), function(rowIdx){
        tmp_lumpIdx <- lumpSummary$lumpIdx[rowIdx]
        tmp_LumpPattern <- lumpLog %>% filter(lumpIdx == tmp_lumpIdx)
        # relatively old function
        return(energyEstimation_caseCook(data, Hz, tmp_LumpPattern, range_energyEsti, min_usageTime))
      }) # %>% ldply(., .id= NULL)
      
      return(cbind(estiOutput, sigNum = rep(lumpSummary$sum, each = 2)))
    })
    
    # examine proper metric
    validMetric <- ldply(energyEsti, function(energyData){
      
      properUsage <- energyData %>% filter(SecON >= min_usageTime)
      if(nrow(properUsage) == 0) return(data.frame(sigNum_org = 0,sigNum_proper = 0,nRow_org = 0,nRow_proper = 0) ) else{
        return(data.frame(sigNum_org = min(energyData$sigNum), sigNum_proper = min(properUsage$sigNum), nRow_org = nrow(energyData)/2, nRow_proper = nrow(properUsage)/2) )
      }
    })

    # candMetric <- validMetric %>% mutate(chosenMetric = sigNum_proper * nRow_proper) %>% .$chosenMetric 
    candMetric <- validMetric %>% mutate(chosenMetric = sigNum_proper) %>% .$chosenMetric 
    
    # print("ricecooker sigNum:")    
    # print(max(candMetric)) # 10 to 51
    
    if(all(candMetric == 0)){
      print('no cooking candidates with proper energy usage')

    } else {
      tmp_ID <- which.max(candMetric)
      proper_lumpNum <- validMetric$nRow_proper[tmp_ID]
      proper_minSigNum <- validMetric$sigNum_proper[tmp_ID]   
      ID_cookSig <- shortCand[tmp_ID]
      cookEnergy <- energyEsti[[tmp_ID]] %>% filter(SecON >= min_usageTime) %>% 
                    summarise(startTime = median(startTime), endTime = median(endTime), SecON = median(SecON), whON = median(whON), apON = median(apON), minSigNum = min(sigNum))
      
      # to split the cook signals into two categories (later in consumption output)
      cookSummary <- major_SigInfo$candSummary[[ID_cookSig]]
      log_cookEnergy <-energyEsti[[tmp_ID]]
    }
  }

  # step 3] process long(1) signals + residual short(0) signals to detect warming mode
  ### 3-1] choose candidates
  ID_warmSig <- 0
  if(!any(sigClass == 1)) {
    considerCand <- which(sigClass == 0)
    considerCand <- considerCand[considerCand != ID_cookSig]
    
    if(length(considerCand) == 0){
      print('no warming candidates in major signals')
      longCand <- NULL
    } else { # conservative approach
      # print('using residual signals for warming candidates')
      # longCand <- considerCand
      print('no warming candidates in major signals')
      longCand <- NULL
    }
    
  } else { # There exist long signals
    print('using proper long signals for warming candidates')
    longCand <- which(sigClass == 1)
  }
  
  ### 3-2] refine signal
  if(!is.null(longCand)){
    
    refinedLog <- lapply(longCand, function(sigIdx){
      lumpSummary <- major_SigInfo$candSummary[[sigIdx]]
      
      lumpLog <- major_SigInfo$candPattern[[sigIdx]]
      h2Min <- min(lumpLog$h2)
      h2Max <- max(lumpLog$h2)
      lump_timeIdx <- lumpLog %>% group_by(lumpIdx) %>% summarise(min_startIdx = min(start.idx), max_endIdx = max(end.idx))
      
      reduced_entireLog <- major_SigInfo$wholePattern[[ major_SigInfo$candGroup[[sigIdx]] ]] %>% filter(h2 >= h2Min) %>% filter(h2 <= h2Max)
      
      newOutput <- apply(lump_timeIdx, 1, function(x) return(reduced_entireLog %>% filter(start.idx >= as.numeric(x['min_startIdx']), end.idx <= as.numeric(x['max_endIdx'])) ))
      return( bind_rows(newOutput, .id = NULL) )  
    })
    
    idx_chosenCand <- which.max(sapply(refinedLog, nrow))
    ID_warmSig <- longCand[idx_chosenCand]
    refined_warmLog <- refinedLog[[idx_chosenCand]]
    
    # energy estimation (w/ original pattern log)
    warmEnergy_Info <- energyEstimation_caseWarm(data, Hz, pattern = major_SigInfo$candPattern[[ID_warmSig]], thres = min_apThres)
  }

  ### step 4] examine correlation (NOT REALIZED YET)
    
  # result frame
  # str.t <- as.character(min(data$timestamp))
  # end.t <- as.character(max(data$timestamp))
  # answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"),
  #                                           floor_date(max(data$timestamp), unit = "second"), 'secs'))
  # answer.log <- data.frame( answer.log, p= rep(0,nrow(answer.log)), q= rep(0,nrow(answer.log)))
  
  # cook
  if(ID_cookSig != 0){
    # insert cooking signals (for results)
    summary_cookTime <- cookSummary %>% mutate(cookTime = lump.start + dseconds(lump.duration/2)) %>% select(cookTime)
    valid_cookLump <- log_cookEnergy %>% mutate(validity = (SecON >= min_usageTime)) %>% bind_cols(data.frame(lumpIdx = rep(seq(nrow(summary_cookTime)), each = 2))) %>%
                      group_by(lumpIdx) %>% summarise(validity = any(validity), SecON = median(SecON), whON = median(whON), apON = median(apON)) %>% bind_cols(summary_cookTime)
    
    ### cook signals w/ proper energy consumption
    if(any(valid_cookLump$validity == TRUE)){
      proper_cookInfo <- valid_cookLump %>% filter(validity == TRUE)
      cooking.start.timestamp <- floor_date(proper_cookInfo$cookTime, unit = "second")
      cooking.start.idx <- which(answer.log$timestamp %in% cooking.start.timestamp)
      cooking.on.idx <- unique(unlist(lapply( cooking.start.idx, function(x) x + (-round(cookEnergy$startTime):(round(cookEnergy$endTime)-1)) )))
      ### avoid overflow of the idx
      cooking.on.idx <- cooking.on.idx[(cooking.on.idx >= 1) & (cooking.on.idx <= nrow(answer.log)) ]
      answer.log$p[cooking.on.idx] <- answer.log$p[cooking.on.idx] - (cookEnergy$whON*3600/(cookEnergy$startTime + cookEnergy$endTime ))
      # end
    }
    
    ### cook signals w/ minor energy consumption
    if(any(valid_cookLump$validity == FALSE)){
      proper_cookInfo <- valid_cookLump %>% filter(validity == FALSE) %>% mutate(cookTime = floor_date(cookTime, unit = "second"))
      cooking.start.idx <- sapply(proper_cookInfo$cookTime, function(x) which(x == answer.log$timestamp) )
      cooking.on.idx <- unique(unlist(mapply(function(x,sec) x + (-round(sec/2):(round(sec/2)-1)), cooking.start.idx, proper_cookInfo$SecON, SIMPLIFY = FALSE)))
      tmp_estiAP <- proper_cookInfo %>% summarise(estiAP = median(proper_cookInfo$apON)) %>% .$estiAP
      ### avoid overflow of the idx
      cooking.on.idx <- cooking.on.idx[(cooking.on.idx >= 1) & (cooking.on.idx <= nrow(answer.log)) ]
      answer.log$p[cooking.on.idx] <- answer.log$p[cooking.on.idx] - tmp_estiAP
      # end
    }

  }

  
  # warm
  if(ID_warmSig != 0){
    # insert warming signals (for results)
    new_timestamp <- floor_date(refined_warmLog$start.timestamp, unit = "second")
    
    warming.start.idx <- which(answer.log$timestamp %in% new_timestamp)
    if (round(warmEnergy_Info$sec) != 0){
      warming.on.idx <- unique(unlist(lapply( warming.start.idx, function(x) x + 0:(round(warmEnergy_Info$sec)-1) )))
    }
    
    ### avoid overflow of the idx
    warming.on.idx <- warming.on.idx[warming.on.idx <= nrow(answer.log)]
    
    if (round(warmEnergy_Info$sec) != 0){
      
      warming.usage_pwr <- warmEnergy_Info$watt
    } 
    
    answer.log$p[warming.on.idx] <- answer.log$p[warming.on.idx] + warming.usage_pwr * 2
    # end
  }
  
  return(answer.log)

}

energyEstimation_caseCook <- function(data, Hz, end.pattern, max_t.sec, thres_usageSec){
  
  # start and end point of a lump
  cookEdgeInfo <- list()
  cookEdgeInfo[[1]] <- head(end.pattern, 1)
  cookEdgeInfo[[2]] <- tail(end.pattern, 1)
  middleTime <- cookEdgeInfo[[1]]$end.timestamp + dseconds((as.numeric(cookEdgeInfo[[2]]$end.timestamp) - as.numeric(cookEdgeInfo[[1]]$end.timestamp))/2)
  
  # find the amount of energy (with its proper position)
  # for exmaple: x <- cookEdgeInfo[[2]]
  energyInfo <- ldply(cookEdgeInfo, function(x){
    cookRange <- (x$end.timestamp-dseconds(max_t.sec) ) %--% (x$end.timestamp+dseconds(max_t.sec) )
    processData <- data %>% filter( .$timestamp %within% cookRange, .$active_power >= x$ap.thres)
    effSecON <- nrow(processData) / Hz
    splitInterval <- cumsum( c(TRUE, difftime( tail(processData$timestamp,-1),
                                               head(processData$timestamp,-1), units='secs') >= thres_usageSec))
    splitLumps <- bind_cols(processData, data.frame(lumpIdx = splitInterval)) %>% group_by(lumpIdx) %>% 
                  summarise(sum=n(), lump.start = min(timestamp), lump.end = max(timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start))
    maxLumpIdx <- which.max(splitLumps$lump.duration)
    if(splitLumps$lump.duration[maxLumpIdx] < effSecON) return(data.frame(startTime = as.numeric(middleTime) - as.numeric(min(processData$timestamp)), 
                                                                          endTime = as.numeric(max(processData$timestamp)) - as.numeric(middleTime), 
                                                                          SecON = effSecON, 
                                                                          whON = effSecON*x$h2/3600, apON = x$h2)) 
    else return(data.frame(startTime = as.numeric(middleTime) - as.numeric(splitLumps$lump.start[maxLumpIdx]), 
                           endTime = as.numeric(splitLumps$lump.end[maxLumpIdx]) - as.numeric(middleTime), 
                           SecON = effSecON, 
                           whON = effSecON*x$h2/3600, apON = x$h2))
    
  })
  
  # cookEdgeInfo <- bind_rows(cookEdgeInfo) %>% bind_cols(., energyInfo)
  return(energyInfo)
}

energyEstimation_caseWarm <- function(data, Hz, pattern, thres){
  
  pattern <- pattern %>% arrange(start.idx)
  
  sigDuration <- pmax(0, tail(pattern$end.idx,-1) - head(pattern$end.idx,-1))
  
  # start(ap reference), end, and time duration (until end) for each slot
  ### case w/ rising edge
  # tbl_sigProcess <- mapply( function(s,eVec){c(s,eVec)}, as.list( head(pattern$start.idx,-1)),
  #                           mapply( function(x,l) (x+c(0:l)), head(pattern$end.idx, -1), sigDuration) )
  ### case w/ falling edge
  tbl_sigProcess <- mapply( function(x,l) (x+c(0:l)), head(pattern$end.idx, -1), sigDuration)

  # compute energy for the whole signals  
  slotInfo <- ldply( tbl_sigProcess, function(x){ 
    apDiff <- data$active_power[x[-1]] - data$active_power[x[1]]
    validPwr <- apDiff[which(apDiff > thres)]
    len <- length(validPwr)
    if(len != 0 ) W <- sum(validPwr)/len else W <- 0
    return(data.frame(length = len, energy = W))
  })
  
  slotInfo <- slotInfo %>% mutate( sec = pmax(1, floor(length/Hz)), wattHr = length * energy / Hz)
  
  med_wattHr <- median(slotInfo$wattHr)
  med_sec <- median(slotInfo$sec)
  wattResult <- med_wattHr / med_sec
  
  return(data.frame(watt = wattResult, sec = med_sec) )
}


# cookSigSearch_15Hz <- function(data, end.pattern, min_watt, max_watt, min_t.sec, max_t.sec, min_fluc.num){
# 
#   e.pattern.cooking <- end.pattern %>%
#     filter( h2 >= -max_watt ) %>%
#     filter( h2 <= -min_watt )
# 
#   if ( nrow(e.pattern.cooking) == 0 ){
#     print("There is no candidate for the cooking signals")
#     return( e.pattern.cooking )
#   }
# 
#   e.pattern.cooking <- e.pattern.cooking %>%
#     mutate( ap.thres = ( data$active_power[start.idx] +
#                            data$active_power[end.idx]) / 2 )
# 
#   # to improve speed
#   data <- data[, !duplicated(colnames(data))] %>%
#     filter( active_power >= min(e.pattern.cooking$ap.thres) ) %>%
#     filter( timestamp >= min(e.pattern.cooking$end.timestamp) - dseconds(max_t.sec) ) %>%
#     filter( timestamp <= max(e.pattern.cooking$end.timestamp) ) %>%
#     select( active_power, timestamp )
# 
#   splitInterval <- cumsum( c(TRUE, difftime( tail(e.pattern.cooking$end.timestamp,-1),
#                                              head(e.pattern.cooking$end.timestamp,-1), units='secs') > max_t.sec))
# 
#   e.pattern.cooking.split <- split.data.frame( e.pattern.cooking, splitInterval )
#   e.pattern.cooking.split.head <-
#     ldply( e.pattern.cooking.split, head, n = 1 )
# 
#   data.split <- split.data.frame( data, findInterval( data$timestamp,
#                                                       e.pattern.cooking.split.head$end.timestamp - dseconds(max_t.sec)))
# 
#   eff.energy <- mapply( function( subData, subPattern){
#     ddply( subPattern, .( end.timestamp ),
#            function(df){
#              res <- subData %>%
#                filter( active_power >= df$ap.thres ) %>%
#                filter( timestamp >= (df$end.timestamp-dseconds(max_t.sec)) ) %>%
#                filter( timestamp <= df$end.timestamp ) %>%
#                nrow
#              c('row.num' = res)
#            } )
#   }, data.split, e.pattern.cooking.split, SIMPLIFY = FALSE )
#   eff.energy <- ldply(eff.energy, .id=NULL)
# 
#   e.pattern.cooking <-
#     merge( e.pattern.cooking, eff.energy, by = 'end.timestamp' ) %>%
#     mutate( energy = row.num * abs(h2), on_sec = row.num/15, Wh = energy/15/3600 ) %>%
#     filter( on_sec >= min_t.sec )
# 
#   fluc.num <- sapply( e.pattern.cooking$start.timestamp,
#                       function( str ){
#                         strIdx <- which( e.pattern.cooking$end.timestamp >= (str-dseconds(max_t.sec) ) )
#                         endIdx <- which( e.pattern.cooking$end.timestamp <= (str) )
#                         if( length(strIdx) == 0 || length(endIdx) == 0 ) return(0)
#                         max(endIdx) - min(strIdx) + 1
#                       }, simplify = TRUE, USE.NAMES = FALSE )
# 
#   e.pattern.cooking <-
#     e.pattern.cooking %>%
#     mutate( fluctuation = fluc.num ) %>%
#     filter( fluctuation >= min_fluc.num )
# 
#   return( e.pattern.cooking )
# }
# 
# sliceTimeIdx_Quick <- function(time_stamp, str.t, end.t){
# 
#   start_idx <- which(str.t <= time_stamp)
#   end_idx <- which(time_stamp <= end.t)
# 
# 
#   if ((length(start_idx) == 0) || (length(end_idx) == 0)) return(c())
# 
#   min_start_idx <- min(start_idx)
#   max_end_idx <- max(end_idx)
# 
#   if (max_end_idx < min_start_idx) return(c())
# 
#   return(min_start_idx:max_end_idx)
# }
# 
# # cookSigSearch_15Hz_obsolete <- function(data, end.pattern, min_watt, max_watt, min_t.sec, max_t.sec, min_fluc.num){
# #
# #   e.pattern.cooking <- subset(end.pattern, h2 >= -max_watt & h2 <= -min_watt)
# #
# #   if ( nrow(e.pattern.cooking) == 0 ) {
# #     print("There is no candidate for the cooking signals")
# #     return(e.pattern.cooking)
# #   } else {
# #
# #     center.ap <- apply(e.pattern.cooking, 1, function(x) (data$active_power[as.numeric(x['start.idx'])] +
# #                                                             data$active_power[as.numeric(x['end.idx'])]   )/2 )
# #
# #     e.pattern.cooking <- data.frame(e.pattern.cooking, ap.thres = center.ap)
# #
# #     # (old one)
# #     #     eff.energy <- rbind.fill(mapply( function(t,p,d) { z <- subset(data, timestamp %within% interval(t-dseconds(max_t.sec), t ) & active_power >= p )
# #     #                                                        num <- nrow(z)
# #     #                                                        return(data.frame(row.num = num, energy = num*d ))}, e.pattern.cooking$end.timestamp, e.pattern.cooking$ap.thres,
# #     #                                                        abs(e.pattern.cooking$h2), SIMPLIFY = F))
# #
# #     # (new one: complexity reduction)
# #     eff.energy <- rbind.fill(mapply( function(t,p,d) { DataIdx <- sliceTimeIdx_Quick(data$timestamp, t-dseconds(max_t.sec), t)
# #     pwr_val <- data$active_power[DataIdx]
# #     num <- length( pwr_val[pwr_val >= p] )
# #     return(data.frame(row.num = num, energy = num*d ))}, e.pattern.cooking$end.timestamp, e.pattern.cooking$ap.thres,
# #     abs(e.pattern.cooking$h2), SIMPLIFY = F))
# #
# #     # (old one)
# #     #     fluc.num <- mapply( function(e) { z <- subset(e.pattern.cooking, end.timestamp %within% interval(e-dseconds(max_t.sec), e) )
# #     #                                       return(nrow(z))}, e.pattern.cooking$start.timestamp)
# #
# #     # (new one: complexity reduction)
# #     fluc.num <- mapply( function(e) { DataIdx <- sliceTimeIdx_Quick(e.pattern.cooking$end.timestamp, e-dseconds(max_t.sec), e)
# #     return(length(DataIdx))}, e.pattern.cooking$start.timestamp)
# #
# #     e.pattern.cooking <- data.frame(e.pattern.cooking, on_sec = eff.energy$row.num/15, Wh = eff.energy$energy/15/3600, fluctuation = fluc.num)
# #     e.pattern.cooking <- subset(e.pattern.cooking, (on_sec >= min_t.sec) & (fluctuation >= min_fluc.num))
# #
# #     return(e.pattern.cooking)
# #   }
# #
# # }
# 
# DetectGroup_15Hz_2D_RC <- function(data, resolution, main.g.name, sub.g.name, p_factor, med.time_min, med.time_max, group_size){
# 
#   findGaps <- function(x,n){
#     x <- sort(x)
#     x.diff <- data.frame( val = diff(x), idx = 1:(length(x)-1) )
#     wall.idx <- x.diff$idx[ order( x.diff$val, decreasing=T ) ][1:n] # choose first N gaps
#     result <- sapply( wall.idx, function(i) mean(x[i+c(0,1)]))
#     return( sort(result) )
#   }
# 
#   summarize_timestamp <- function( timestamp ){
#     tsDiff <- diff(as.numeric(sort(timestamp)))
#     data.frame(
#       'sum' = length(timestamp),
#       'min.t'  = min(tsDiff),
#       'min.t2' = quantile(tsDiff, .05),
#       'med.t'  = median(tsDiff),
#       'lost.sig.num' = sum( pmax( round(tsDiff / median(tsDiff)) -1 ,0)))
#   }
# 
#   logSummary <- list()
#   logDetail <- list()
#   o.list <- data
# 
#   r_idx <- 1
#   iter_idx <- 1
#   while( iter_idx <= resolution ){
# 
#     if ( nrow(o.list) +1 <= r_idx ){
#       print("stop resolution increment for searcing groups")
#       break
#     }
# 
#     xValue    <- o.list[,main.g.name] %>% unlist(use.names = FALSE)
#     yValue    <- o.list[, sub.g.name] %>% unlist(use.names = FALSE)
#     o.list$xDivision <- findInterval( xValue, findGaps( xValue, r_idx ) )
#     o.list$yDivision <- findInterval( yValue, findGaps( yValue, r_idx ) )
# 
#     removeData <- o.list %>%
#       group_by( xDivision, yDivision ) %>%
#       filter( n() < group_size )
# 
#     o.list <- o.list %>%
#       group_by( xDivision, yDivision ) %>%
#       filter( n() >= group_size )
# 
#     r_idx <- pmax( r_idx - length(group_size(removeData)), 1 )
# 
#     if (nrow(o.list) > 0){ # in case 'o.list' is empty
# 
#       group.info <- o.list %>%
#         do( summarize_timestamp(.$start.timestamp) ) %>%
#         mutate( min.med.rate  = min.t/med.t,
#                 min.med.rate2 = min.t2/med.t )
# 
#       group.info <- merge( group.info, o.list %>%
#                              summarise_each_( funs(median), c('h1','delta','sub.delta') )) %>%
#         rename( med.h1 = h1, med.d = delta, med.sub.d = sub.delta ) %>%
#         mutate( lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
#     }
# 
#     ### effective group (conservative search by default)
#     if (nrow(group.info) > 0){
# 
#       eff_group.info <- group.info %>%
#         filter( med.t >= med.time_min ) %>%
#         filter( med.t <= med.time_max ) %>%
#         filter( min.med.rate2 >= p_factor ) %>%
#         filter( xDivision != 0 )
# 
#       if( nrow(eff_group.info) > 0 ){
# 
#         logSummary[[length(logSummary)+1]] <- eff_group.info
# 
#         for(g_idx in 1:nrow(eff_group.info) ){
#           logDetail[[length(logDetail) +1]] <- o.list %>%
#             filter( xDivision == eff_group.info$xDivision[g_idx] & yDivision == eff_group.info$yDivision[g_idx] )
# 
#           o.list <- o.list %>%
#             filter( !(xDivision == eff_group.info$xDivision[g_idx] & yDivision == eff_group.info$yDivision[g_idx]) )
#         }
#       }
# 
#       # index increment for while loop
#       r_idx <- pmax( r_idx - nrow(eff_group.info) + 1, 1 )
#     } else {
#       r_idx <- r_idx + 1
#     }
#     iter_idx <- iter_idx + 1
# 
#   }
# 
#   # return information
#   if (length(logSummary) != 0) {
#     logDetail[[length(logDetail) +1]] <- rbind.fill(logSummary)
#     return(logDetail)
#   } else {
#     return(list())
#   }
# 
# }
# 
# # DetectGroup_15Hz_2D_RC_obsolete <- function(data, resolution, main.g.name, sub.g.name, p_factor, med.time_min, med.time_max, group_size){
# #
# #   logSummary <- list()
# #   logDetail <- list()
# #   o.list <- data
# #
# #   r_idx <- 1
# #   iter_idx <- 1
# #   while( iter_idx <= resolution ){
# #
# #     # print(nrow(ordered.delta))
# #     # for debug
# #     if ( nrow(o.list) +1 <= r_idx ){
# #       print("Resolution increment for searcing groups has stopped..")
# #       break
# #     }
# #
# #     # examine g1 : should be 0 to Inf!
# #     o.list <- o.list[order(o.list[[main.g.name]]),]
# #
# #     ### detect group idx (i.e., wall for each group)
# #     g1.diff <- data.frame( diff.g1 = diff(o.list[[main.g.name]]), o.idx = seq(1, nrow(o.list)-1) )
# #     o.g1.diff <- g1.diff[order(g1.diff$diff.g1, decreasing = TRUE),]
# #     wall.idx <- o.g1.diff$o.idx[1:r_idx]
# #     wall.idx <- wall.idx[order(wall.idx)]
# #
# #     ### identify group 1
# #     g1 <- cut( seq(1, nrow(o.list)), c(0, wall.idx, nrow(o.list)) ,labels= FALSE ) -1
# #     o.list <- cbind(o.list, g1)
# #
# #     # examine g2 : reactive power would be fine
# #     o.list <- o.list[order(o.list[[sub.g.name]]),]
# #
# #     ### detect group idx (i.e., wall for each group)
# #     g2.diff <- data.frame( diff.g2 = diff(o.list[[sub.g.name]]), o.idx = seq(1, nrow(o.list)-1) )
# #     o.g2.diff <- g2.diff[order(g2.diff$diff.g2, decreasing = TRUE),]
# #     wall.idx <- o.g2.diff$o.idx[1:r_idx]
# #     wall.idx <- wall.idx[order(wall.idx)]
# #
# #     ### identify group 2
# #     g2 <- cut( seq(1, nrow(o.list)), c(0, wall.idx, nrow(o.list)) ,labels= FALSE ) -1
# #     o.list <- cbind(o.list, g2)
# #     o.list <- o.list[order(o.list$start.idx), ]
# #
# #     # summary for each group
# #     group.info <- ddply(o.list, .(g1, g2), transform, tmp_g.sum = length(g1))
# #
# #     ### exclude the info from isolated elements
# #     tmp_small.group <- subset(group.info, tmp_g.sum < group_size)
# #     r_idx <- r_idx -nrow( ddply(tmp_small.group, .(g1, g2), summarize, tmp_sum = 'del') )
# #     if(r_idx < 1 ){r_idx <- 1}
# #
# #     o.list <- subset(group.info, tmp_g.sum >= group_size)
# #     group.info <- ddply(o.list, .(g1, g2), summarize, sum = length(g1), min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))),
# #                         min.med.rate = min.t/med.t, min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),  med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
# #                         lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# #
# #     ### effective group (conservative search by default)
# #     if (nrow(group.info) > 0){
# #
# #       eff_group.info <- subset(group.info, med.t >= med.time_min & med.t <= med.time_max & min.med.rate2 >= p_factor & g1 != 0 )
# #
# #       if(nrow(eff_group.info) >0){
# #
# #         logSummary[[length(logSummary)+1]] <- eff_group.info
# #
# #         for(g_idx in 1:nrow(eff_group.info) ){
# #
# #           logDetail[[length(logDetail) +1]] <- o.list[((o.list$g1 == eff_group.info$g1[g_idx]) &
# #                                                          (o.list$g2 == eff_group.info$g2[g_idx])), ]
# #           logDetail[[length(logDetail)]]$g1 <- NULL
# #           logDetail[[length(logDetail)]]$g2 <- NULL
# #           logDetail[[length(logDetail)]]$tmp_g.sum <- NULL
# #
# #
# #           o.list <- o.list[!((o.list$g1 == eff_group.info$g1[g_idx]) &
# #                                (o.list$g2 == eff_group.info$g2[g_idx])), ]
# #         }
# #       }
# #
# #       # index increment for while loop
# #       r_idx <- r_idx -nrow(eff_group.info) +1
# #       if(r_idx < 1 ){r_idx <- 1}
# #     } else {
# #
# #       r_idx <- r_idx + 1
# #     }
# #
# #     iter_idx <- iter_idx + 1
# #     o.list$g1 <- NULL
# #     o.list$g2 <- NULL
# #     o.list$tmp_g.sum <- NULL
# #   }
# #
# #   # return information
# #   if (length(logSummary) != 0) {
# #     logDetail[[length(logDetail) +1]] <- rbind.fill(logSummary)
# #     return(logDetail)
# #   } else {
# #     return(list())
# #   }
# #
# # }
# 
# parallel_warmSig_search_15Hz <- function(data, w.pattern.s, ap_h1_min, ap_h1_max, ap_delta_min, ap_delta_max, rp_delta_thres, main.g.name, sub.g.name,
#                                          periodicity, search_iter, w.med_sec.min, w.med_sec.max, eff_size, eff_group_size, w.lump.gap_sec.max,
#                                          w.lump.med_sec.margin, w.cand.DB, DB_tolerance, DB_rate_quantile, DB_timespan) {
# 
#   options(scipen=99)
# 
#   BinToDec <- function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
# 
#   chosen.w.pattern <- subset(w.pattern.s, h1 >= ap_h1_min & h1 <= ap_h1_max & delta >= ap_delta_min & delta <= ap_delta_max &
#                                abs(sub.delta) <= rp_delta_thres)
# 
#   s.info <- DetectGroup_15Hz_2D_RC(chosen.w.pattern, resolution = search_iter, main.g.name, sub.g.name,
#                                    p_factor = periodicity, med.time_min = w.med_sec.min, med.time_max = w.med_sec.max, group_size = eff_group_size)
# 
#   if (length(s.info) == 0) {
#     print("Warming mode detection for ricecooker has failed.")
#     return(list())
#   }
# 
#   print("groups for the warming mode:")
#   s.list <- s.info[[length(s.info)]]
#   rownames(s.list) <- 1:(length(s.info)-1)
# 
#   # ordering for warming candidates
#   s.list <- s.list[order(s.list$sum, decreasing = T),]
#   print(s.list)
# 
#   ID_type <- list()
#   ID_total.log <- list()
#   ID_log.summary <- list()
#   # examine candidates for the warming mode
#   for (w_idx in 1:nrow(s.list)) {
# 
#     print(w_idx)
#     chosen.s.info <- s.info[[as.numeric(row.names(s.list[w_idx,])) ]]
# 
#     lump.idx <- c(0, which(diff(as.numeric(chosen.s.info$start.timestamp)) >= w.lump.gap_sec.max), nrow(chosen.s.info))
#     chosen.lump <- which(diff(lump.idx) >= eff_size)
# 
#     if (length(chosen.lump) >0) {
# 
#       lump_log <- list()
# 
#       for (l_idx in 1:length(chosen.lump)) {
# 
#         tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
#         tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
# 
#         lump_log[[length(lump_log)+1]] <- data.frame(chosen.s.info[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx)
#       }
# 
#       lump_log.total <- rbind.fill(lump_log)
# 
#       # examine warming candidates
# 
#       ### 1) refine lumps w/ median sig. period
# 
#       lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))),
#                             min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), DB_rate_quantile)/med.t )
# 
#       ### 2) retrive DB
#       DB_result <- mapply(function(m,r) { z <- subset( lump_summary, ( abs(med.t - m) <= m*DB_tolerance ) &
#                                                          ( abs(min.med.rate2 - r) <= r*DB_tolerance ) )
#       if (nrow(z) != 0) return(T) else return(F) },
#       w.cand.DB$med.t, w.cand.DB$min.med.rate2, SIMPLIFY = T)
# 
#       if(length(which(DB_result == T)) != 0){
# 
#         # identified_type <- max(which(DB_result == T))
#         # ID_type[[length(ID_type)+1]] <- as.numeric(paste0(as.character(identified_type), paste(as.integer(DB_result), collapse = "")))
# 
#         ID_type[[length(ID_type)+1]] <- BinToDec( as.numeric( paste(as.integer(DB_result), collapse = "")))
#         cat("warming signal identified: type ", ID_type[[length(ID_type)]] , "\n")
# 
#         lump_summary <- rbind.fill(mapply(function(c,m) {if (c) return(subset(lump_summary, med.t >= (m - w.lump.med_sec.margin) &
#                                                                                 med.t <= (m + w.lump.med_sec.margin ) ))
#           else NULL }, DB_result, w.cand.DB$med.t, SIMPLIFY = F))
# 
#       } else {
# 
#         ID_type[[length(ID_type)+1]] <- 0
# 
#         tmp.med_time <- lump_summary$med.t[which.max(lump_summary$sum)]
#         lump_summary <- subset(lump_summary, med.t >= (tmp.med_time - w.lump.med_sec.margin) &
#                                  med.t <= (tmp.med_time + w.lump.med_sec.margin) )
#       }
#       ID_total.log[[length(ID_total.log)+1]] <- rbind.fill(lump_log[lump_summary$lump_idx])
#       ID_log.summary[[length(ID_log.summary)+1]] <- ddply(ID_total.log[[length(ID_total.log)]], .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
#                                                           min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
#                                                           med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                                                           lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#     }
#   }
# 
#   return(list(type = ID_type, total.log = ID_total.log, log.summary = ID_log.summary))
# }
# 
# parallel_cookSig_search_15Hz <- function(data, cooking.min_watt, cooking.max_watt, cooking.max_t.sec, cooking.min_fluc.num, periodicity,
#                                          search_iter, c.med_sec.min = w.med_sec.min, c.med_sec.max = w.med_sec.max, eff_size, c.lump.gap_sec.max,
#                                          c.lump.med_sec.margin, c.cand.DB, CDB_tolerance) {
# 
#   c.pattern.e <- DetectPattern_1Hz_new(data, position = "end", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
# 
#   e.pattern.cooking <- cookSigSearch_15Hz(data, end.pattern = c.pattern.e,
#                                           min_watt = cooking.min_watt, max_watt = cooking.max_watt, min_t.sec = -Inf, max_t.sec = cooking.max_t.sec,
#                                           min_fluc.num = cooking.min_fluc.num)
# 
#   e.pattern.cooking$delta <- -e.pattern.cooking$delta
#   e.info.cooking <- DetectGroup_15Hz_2D_RC(e.pattern.cooking, resolution = search_iter, main.g.name = "delta", sub.g.name = "sub.delta",
#                                            p_factor = periodicity, med.time_min = c.med_sec.min, med.time_max = c.med_sec.max, group_size = eff_size)
# 
#   if (length(e.info.cooking) == 0) {
#     print("group search for retrieving cook DB has failed.")
#     return(list(entire_pattern = e.pattern.cooking))
#   }
# 
#   print("groups for the cooking mode:")
#   e.list.cooking <- e.info.cooking[[length(e.info.cooking)]]
#   rownames(e.list.cooking) <- 1:(length(e.info.cooking)-1)
# 
#   # ordering for cooking candidates
#   e.list.cooking <- e.list.cooking[order(e.list.cooking$sum, decreasing = T), ]
#   print(e.list.cooking)
# 
#   CID_type <- list()
#   CID_total.log <- list()
#   CID_log.summary <- list()
#   CID_energy <- list()
# 
#   # examine candidates for the cooking mode
#   for (c_idx in 1:nrow(e.list.cooking)) {
# 
#     print(c_idx)
#     chosen.e.info.cooking <- e.info.cooking[[as.numeric(row.names(e.list.cooking[c_idx,])) ]]
# 
#     lump.idx <- c(0, which(diff(as.numeric(chosen.e.info.cooking$start.timestamp)) >= c.lump.gap_sec.max), nrow(chosen.e.info.cooking))
#     chosen.lump <- which(diff(lump.idx) >= eff_size)
# 
#     if (length(chosen.lump) >0) {
# 
#       lump_log <- list()
# 
#       for (l_idx in 1:length(chosen.lump)) {
# 
#         tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
#         tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
# 
#         lump_log[[length(lump_log)+1]] <- data.frame(chosen.e.info.cooking[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx)
#       }
# 
#       lump_log.total <- rbind.fill(lump_log)
# 
#       # examine cooking candidates
# 
#       ### 1) refine lumps w/ median sig. period
# 
#       lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t )
# 
#       tmp.med_time <- lump_summary$med.t[which.max(lump_summary$sum)]
# 
#       lump_summary <- subset(lump_summary, med.t >= (tmp.med_time * (1 - c.lump.med_sec.margin) ) &
#                                med.t <= (tmp.med_time * (1 + c.lump.med_sec.margin) ) )
# 
#       ### 2) retrive DB
#       DB_result <- mapply(function(m) { z <- subset( lump_summary, ( abs(med.t - m) <= m*CDB_tolerance ) )
#       if (nrow(z) != 0) return(T) else return(F) }, c.cand.DB$med.t, SIMPLIFY = T)
# 
#       if(length(which(DB_result == T)) != 0){
# 
#         CID_type[[length(CID_type)+1]] <- min(which(DB_result == T))
#         cat("cooking signal identified: type", CID_type[[length(CID_type)]] , "\n")
# 
#         # energy consumption estimation
#         CID_energy.esti <- rbind.fill( mapply(function(l) energyEstimation_cookDB_15Hz(data, lump_log[[l]], cooking.max_t.sec),
#                                               lump_summary$lump_idx, SIMPLIFY = F) )
# 
#         CID_energy[[length(CID_energy)+1]] <- data.frame(cook_esti_sec = median(CID_energy.esti$cook_esti_sec), cook_esti_Wh = median(CID_energy.esti$cook_esti_Wh),
#                                                          cook_esti_ap = median(CID_energy.esti$cook_esti_ap) )
# 
#       } else {
# 
#         CID_type[[length(CID_type)+1]] <- 0
#         CID_energy[[length(CID_energy)+1]] <- data.frame(cook_esti_sec = 0, cook_esti_Wh = 0)
#       }
# 
#       CID_total.log[[length(CID_total.log)+1]] <- rbind.fill(lump_log[lump_summary$lump_idx])
#       CID_log.summary[[length(CID_log.summary)+1]] <- ddply(CID_total.log[[length(CID_total.log)]], .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
#                                                             min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
#                                                             med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                                                             lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#     }
#   }
# 
#   return(list(type = CID_type, total.log = CID_total.log, log.summary = CID_log.summary, energy = CID_energy, entire_pattern = e.pattern.cooking))
# 
# }
# 
# sig.Orchestration_RC_15Hz <- function(data, warm.info_low, warm.info_high, cook.info, w.cand.DB, c.cand.DB, DB_timespan, cooking.forward_search, cooking.backward_search,
#                                       lump_merge, cooking.min_sec, margin_sub.delta, margin_delta, margin_sec) {
#   valid_info <- list()
# 
#   # warming mode selection
#   if (any(unlist(warm.info_low$type) > 0)) {
#     warm.info <- warm.info_low
#     valid_info$warm_class <- 1
#   } else if (any(unlist(warm.info_high$type) > 0)) {
#     warm.info <- warm.info_high
#     valid_info$warm_class <- 2
#   } else if (length(warm.info_low$type) != 0) {
#     warm.info <- warm.info_low
#     valid_info$warm_class <- 1
#   } else if (length(warm.info_high$type) != 0) {
#     warm.info <- warm.info_high
#     valid_info$warm_class <- 2
#   } else {
#     print("Lump search for the warming mode has failed.")
#     return(list())
#   }
# 
#   data.start_time <- min(data$timestamp)
#   data.end_time <- max(data$timestamp)
# 
#   # effect of "valid.cook_DB == T"
#   # 1) When "valid.warm_DB == T", one of those will be chosen as the base warming signal REGARDLESS OF its (short) time duration.
#   # 2) ??(NOT SURE) If the signals in "valid.cook_DB" correspond to the base warming mode, they are used for COMPENSATION.
# 
#   valid.warm_DB <- any(unlist(warm.info_low$type) > 0) | any(unlist(warm.info_high$type) > 0)
#   valid.cook_DB <- any(unlist(cook.info$type) > 0)
# 
#   cat("case: valid warming signal-->", valid.warm_DB, ", valid cooking signal-->", valid.cook_DB, "\n")
# 
#   # step 1) choose base (i.e., the most significant) warming signal
#   all.warming.duration <- sapply( 1:length(warm.info$type), function(x) sum(warm.info$log.summary[[x]]$lump.duration) )
# 
#   max_duration <- max(all.warming.duration)
#   max_duration.idx <- which.max(all.warming.duration)
# 
#   if (valid.warm_DB) {
# 
#     DB_cand.idx <- which( unlist(warm.info$type) > 0)
#     DB_cand.duration <- mapply(function(l) sum(warm.info$log.summary[[l]]$lump.duration), DB_cand.idx, SIMPLIFY = T)
# 
#     DB_chosen.idx <- DB_cand.idx[which.max(DB_cand.duration)]
#     DB_chosen.duration <- max(DB_cand.duration)
# 
#     if (valid.cook_DB || (max_duration*DB_timespan <= DB_chosen.duration) ) {
# 
#       print("DB has found a proper warming signal.")
#       chosen_warm.idx <- DB_chosen.idx
#     } else {
# 
#       print("DB has found signals for warming mode, but it was not used because the duration is too short.")
#       chosen_warm.idx <- max_duration.idx
#     }
# 
#   } else {
# 
#     print("DB does not work: choose the longest signal as the warming mode.")
#     chosen_warm.idx <- max_duration.idx
#   }
# 
#   base_warm.summary <- warm.info$log.summary[[chosen_warm.idx]]
# 
#   # step 2) (NOT REALIZED YET!) refine the identified warming signal
# 
#   # step 3) choose cooking signal
# 
#   ### step 3-1) determine search range on each (warming) lump
#   cooking.search_range <- pmin(cooking.forward_search*3600, as.numeric(base_warm.summary$lump.start) - as.numeric( c(data.start_time, base_warm.summary$lump.end[-nrow(base_warm.summary)])))
# 
#   cooking.eff_start <- lump_merge *3600 < cooking.search_range
# 
#   base_warm.summary <- data.frame(base_warm.summary, forward.search_time = cooking.search_range, eff_start = cooking.eff_start)
# 
#   base_warm.idx <- c( which(base_warm.summary$eff_start == T), length(base_warm.summary$eff_start)+1)
# 
#   ###### lump merging
#   ######### case 1
#   eff_start.idx <- base_warm.idx[ which(diff(base_warm.idx) == 1)]
# 
#   eff_warm.summary <- data.frame(forward_search = base_warm.summary$forward.search_time[eff_start.idx],
#                                  lump.start = base_warm.summary$lump.start[eff_start.idx],
#                                  lump.end = base_warm.summary$lump.end[eff_start.idx])
# 
#   ######### case 2
#   eff_start.idx <- base_warm.idx[ which(diff(base_warm.idx) >1)]
# 
#   eff_end.idx <- base_warm.idx[ which(diff(base_warm.idx) >1)+1] -1
# 
#   eff_warm.summary <-rbind(eff_warm.summary, data.frame(forward_search = base_warm.summary$forward.search_time[eff_start.idx],
#                                                         lump.start = base_warm.summary$lump.start[eff_start.idx],
#                                                         lump.end = base_warm.summary$lump.end[eff_end.idx]) )
# 
#   ######### ordering
#   eff_warm.summary <- eff_warm.summary[order(eff_warm.summary$lump.start),]
# 
#   ######### backward search range w/o overlapping
#   cooking.search_range <- pmin(cooking.backward_search*3600, as.numeric( c(eff_warm.summary$lump.start[-1]-dseconds(eff_warm.summary$forward_search[-1]), data.end_time)) -
#                                  as.numeric(eff_warm.summary$lump.start) )
# 
#   eff_warm.summary <- data.frame(eff_warm.summary, backward_search.from.start = cooking.search_range)
#   print(eff_warm.summary)
# 
#   # (1) save valid information
#   valid_info$warm_log <- warm.info$total.log[[chosen_warm.idx]]
#   valid_info$warm_summary <- base_warm.summary
#   valid_info$warm_type <- warm.info$type[[chosen_warm.idx]]
#   valid_info$warm_eff.lump.num <- nrow(eff_warm.summary)
# 
#   ### step 3-2) examine "valid.cook_DB" with the obtained "base warming mode" if possible.
#   if(valid.cook_DB) {
# 
#     CDB_total.log <- list()
#     CDB_total.lump <- list()
#     CDB_valid.lump <- list()
#     CDB_energy <- list()
# 
#     CDB_idx <- which(unlist(cook.info$type) > 0)
# 
#     for (c_idx in 1:length(CDB_idx) ) {
# 
#       tmp_cook.summary <- cook.info$log.summary[[CDB_idx[c_idx]]]
#       CDB_total.log[[length(CDB_total.log)+1]] <- cook.info$total.log[[CDB_idx[c_idx]]]
#       CDB_total.lump[[length(CDB_total.lump)+1]] <- tmp_cook.summary
#       tmp_cook.info <- mapply( function(s,f,b) z <- subset(tmp_cook.summary, lump.end %within% interval(s-dseconds(f), s+dseconds(b) )),
#                                eff_warm.summary$lump.start, eff_warm.summary$forward_search, eff_warm.summary$backward_search.from.start, SIMPLIFY = F)
# 
#       CDB_valid.lump[[length(CDB_valid.lump)+1]] <- length(which( sapply(tmp_cook.info, nrow) != 0))
#       CDB_energy[[length(CDB_energy)+1]] <- cook.info$energy[[CDB_idx[c_idx]]]
#     }
# 
#     # (2) save valid information
#     valid_info$cook_DB.log <- CDB_total.log
#     valid_info$cook_DB.lump.total <- CDB_total.lump
#     valid_info$cook_DB.lump.valid.num <- CDB_valid.lump
#     valid_info$cook_DB.esti.energy <- CDB_energy
# 
#   }
# 
#   ### step 3-3) perform the conventional (i.e., blind detection) algorithm for cooking mode search
#   if(!is.null(cook.info$entire_pattern)){
# 
#     reasonable_cook.pattern <- subset(cook.info$entire_pattern, on_sec >= cooking.min_sec)
# 
#     blind_cook.info <- mapply( function(s,f,b){ z <- subset(reasonable_cook.pattern, end.timestamp %within% interval(s-dseconds(f), s+dseconds(b) ))
#     if (nrow(z) != 0) return(z[which.max(z$on_sec),]) else return(z) },
#     eff_warm.summary$lump.start, eff_warm.summary$forward_search, eff_warm.summary$backward_search.from.start, SIMPLIFY = F)
# 
#     blind_cook.idx <- which(sapply(blind_cook.info, nrow) != 0)
# 
#     if (length(blind_cook.idx) == 0){
# 
#       cat("There is no candidate for the blind COOKING signal detection. \n")
#       cat("If you are sure on the existence of rice cookers, please adjust the cooking signal parameters. \n")
#     } else {
# 
#       blind_cook.info <- rbind.fill(blind_cook.info)
#       blind_cook.info <- blind_cook.info[duplicated(blind_cook.info$start.idx) == F,]
# 
#       ###### refine cooking information
#       new_sub.delta <- median(blind_cook.info$sub.delta)
#       new_delta <- abs(median(blind_cook.info$h2))
# 
#       new_e.pattern.cooking <- subset(reasonable_cook.pattern, sub.delta >= new_sub.delta - margin_sub.delta & sub.delta <= new_sub.delta + margin_sub.delta &
#                                         h2 >= -(new_delta + margin_delta) & h2 <= -(new_delta - margin_delta) )
# 
#       if (nrow(new_e.pattern.cooking) == 0) {
# 
#         cat("There are ambiguities in the blind COOKING signal detection: \n")
#         cat("please input more data so we can label the rice cooker in a reliable (i.e., statistical) way. \n")
#       } else {
# 
#         ######### 1) refine cooking pattern subset
#         blind_cook.info <- mapply( function(s,f,b){ z <- subset(new_e.pattern.cooking, end.timestamp %within% interval(s-dseconds(f), s+dseconds(b) ))
#         if (nrow(z) != 0) return(z[which.max(z$on_sec),]) else return(z) },
#         eff_warm.summary$lump.start, eff_warm.summary$forward_search, eff_warm.summary$backward_search.from.start, SIMPLIFY = F)
# 
#         blind_cook.idx <- which(sapply(blind_cook.info, nrow) != 0)
# 
#         if (length(blind_cook.idx) == 0) {
# 
#           print("Cooking signal refinement has failed after reducing the cooking subset.")
#         } else {
# 
#           ######### 2) refine time gap between warming and cooking
#           blind_cook.info <- rbind.fill(blind_cook.info)
#           blind_cook.info <- blind_cook.info[duplicated(blind_cook.info$start.idx) == F,]
# 
#           r.time_gap <- mapply( function(s) { value <- s - as.numeric(blind_cook.info$end.timestamp)
#           value <- value[value <= cooking.forward_search*3600 & value >= -cooking.backward_search*3600]
#           if (length(value) == 0) return(NA) else(return(value[which.min(abs(value))])) }, as.numeric(eff_warm.summary$lump.start) )
# 
#           new_time.gap <- median(r.time_gap, na.rm = T)
# 
#           ######### 3) adjust time gap between warming and cooking
#           eff_warm.summary$forward_search <- pmin(eff_warm.summary$forward_search, new_time.gap + margin_sec)
#           eff_warm.summary$backward_search.from.start <- pmin(-(new_time.gap - margin_sec), as.numeric( c(eff_warm.summary$lump.start[-1]-dseconds(eff_warm.summary$forward_search[-1]),
#                                                                                                           data.end_time)) - as.numeric(eff_warm.summary$lump.start) )
# 
#           cat("warming mode summary for blind detection: \n")
#           print(eff_warm.summary)
#           blind_cook.info <- mapply( function(s,f,b){ z <- subset(new_e.pattern.cooking, end.timestamp %within% interval(s-dseconds(f), s+dseconds(b) ))
#           if (nrow(z) != 0) return(z[which.max(z$on_sec),]) else return(z) },
#           eff_warm.summary$lump.start, eff_warm.summary$forward_search, eff_warm.summary$backward_search.from.start, SIMPLIFY = F)
# 
#           blind_cook.idx <- which(sapply(blind_cook.info, nrow) != 0)
# 
#           if (length(blind_cook.idx) == 0) {
# 
#             print("Cooking signal refinement has failed after adjusting the time gap between two modes.")
#           } else {
# 
#             ######### 4) yield information
#             blind_cook.info <- rbind.fill(blind_cook.info)
#             blind_cook.info <- blind_cook.info[duplicated(blind_cook.info$start.idx) == F,]
# 
#             cat("resultant cooking mode for blind detection: \n")
#             print(blind_cook.info)
# 
#             # (3) save valid information
#             valid_info$cook_blind.sub.delta <- median(blind_cook.info$sub.delta)
#             valid_info$cook_blind.delta <- abs(median(blind_cook.info$h2))
#             valid_info$cook_blind.fluc_num <- median(blind_cook.info$fluctuation)
#             valid_info$cook_blind.on_sec <- median(blind_cook.info$on_sec)
#             valid_info$cook_blind.time_gap <- new_time.gap
# 
#           }
#         }
# 
#       }
# 
#     }
# 
#   }
# 
# 
#   return(valid_info)
# }
# 
# JSON_result_15Hz <- function (data, information, ap_thres, min_cook.time, parameters, warm_DB_time){
# 
#   str.t <- as.character(min(data$timestamp))
#   end.t <- as.character(max(data$timestamp))
# 
#   # warming energy computation
#   w.usage_info <- energyEstimation_ricecooker(data, pattern = information$warm_log, threshold = ap_thres)
# 
#   if(length(w.usage_info) == 0){
# 
#     print("warming energy computation failure")
#     return(list())
#   }
# 
#   warm_sec <- max(1, floor(w.usage_info[1] / 15) )
#   warm_watt <- (w.usage_info[2] * w.usage_info[1]) /15 /warm_sec
# 
#   # warming information
#   general.info <- c('warm_DB' = any(information$warm_type > 0), 'cook_DB' = ("cook_DB.log" %in% names(information) ),
#                     'cook_blind' = ("cook_blind.delta" %in% names(information) ) )
# 
#   if(general.info['cook_DB'] == F & general.info['cook_blind'] == F){
# 
#     print("cooking mode detection failure")
#     return(list())
#   }
# 
#   warm_log <- information$warm_log
#   warm_med.t <- information$warm_summary$med.t[which.max(information$warm_summary$sum)]
# 
#   warm.info <- c('class' = information$warm_class, 'sample.num' = nrow(warm_log), 'eff.lump' = information$warm_eff.lump.num, 'ap.h1.min' = min(warm_log$h1),
#                  'ap.h1.max' = max(warm_log$h1), 'ap.delta.min' = min(warm_log$delta), 'ap.delta.max' = max(warm_log$delta),
#                  'rp.delta.min' = min(warm_log$sub.delta), 'rp.delta.max' = max(warm_log$sub.delta), 'type' = information$warm_type, 'med.t' = warm_med.t, 'usage.sec' = warm_sec, 'usage.watt' = warm_watt)
# 
#   meta.param <- parameters
# 
#   result <- list(`meta-version` = 1, shape_type = "ricecooker_pattern_scan", general_info = general.info,
#                  warm_info = warm.info, parameter = meta.param, generation_info = list(data_used = list(start = str.t, end = end.t, sampling = 15),
#                                                                                        computed = as.character(Sys.time())))
# 
#   # cooking information
#   ### DB-aided result
#   if (general.info['cook_DB'] == T) {
# 
#     cook_log <- information$cook_DB.log[[1]]
#     cook_summary <- information$cook_DB.lump.total[[1]]
#     cook_valid.num <- information$cook_DB.lump.valid.num[[1]]
#     cook_energy <- information$cook_DB.esti.energy[[1]]
#     information$cook_DB.esti.energy
#     cook_med.t <- cook_summary$med.t[which.max(cook_summary$sum)]
# 
#     cook.DB.info <- c('sample.num' = nrow(cook_log), 'total.lump' = nrow(cook_summary), 'valid.lump' = cook_valid.num, 'ap.h2.min' = min(cook_log$h2),
#                       'ap.h2.max' = max(cook_log$h2), 'rp.delta.min' = min(cook_log$sub.delta), 'rp.delta.max' = max(cook_log$sub.delta),
#                       'med.t' = cook_med.t, 'med.fluc' = median(cook_summary$sum), 'usage.sec' = cook_energy$cook_esti_sec, 'usage.watt' = cook_energy$cook_esti_Wh,
#                       'usage.ap' = cook_energy$cook_esti_ap )
# 
#     result$cook_DB_info <- cook.DB.info
# 
#   }
# 
#   ### blind detection result
#   if (general.info['cook_blind'] == T) {
# 
#     cook.blind.info <- c('ap.delta.med' = information$cook_blind.delta, 'rp.delta.med' = information$cook_blind.sub.delta,
#                          'time_gap' =  information$cook_blind.time_gap, 'fluc_num' = information$cook_blind.fluc_num,
#                          'on_sec' = information$cook_blind.on_sec)
# 
#     result$cook_blind_info <- cook.blind.info
#   }
# 
#   if (general.info['warm_DB'] == T) result$warm_DB_time <- warm_DB_time
# 
#   ### reliability check
# 
#   r_point <- 0.2
#   ### Case 1
#   if( general.info['cook_DB'] == TRUE) r_point <- r_point + 0.4
#   if( general.info['warm_DB'] == TRUE) r_point <- r_point + 0.2
#   result$reliability <- r_point
# 
#   return(result)
# }
# 
# meta2PatternScan.ricecooker_15Hz <- function (data, meta, usage.pwr.adjust = 2) {
# 
#   '%nin%'<- Negate('%in%')
#   options(scipen=99)
# 
#   DecToBin <- function(x) paste(rev(as.integer(intToBits(x))), collapse="")
# 
#   # result frame
#   str.t <- as.character(min(data$timestamp))
#   end.t <- as.character(max(data$timestamp))
#   answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"),
#                                             floor_date(max(data$timestamp), unit = "second"), 'secs'))
#   answer.log <- data.frame( answer.log, p= rep(0,nrow(answer.log)), q= rep(0,nrow(answer.log)))
# 
#   if(length(meta) == 0) return(answer.log)
# 
#   # parameters
#   eff_size <- meta$parameter['eff_size']
#   thres_ap.h <- meta$parameter['thres_ap.h']
#   thres_rp.delta <- meta$parameter['thres_rp.delta']
#   w.lump.gap_sec.max <- meta$parameter['w.lump.gap_sec.max']
#   w.lump.med_sec.margin <- meta$parameter['w.lump.med_sec.margin']
#   c.max_t.sec <- meta$parameter['c.max_t.sec']
#   c.lump.gap_sec.max <- meta$parameter['c.lump.gap_sec.max']
#   c.lump.med_sec.margin <- meta$parameter['c.lump.med_sec.margin']
#   c.lump_merge.max_hr <- meta$parameter['c.lump_merge.max_hr']
#   r.margin_sub.delta <- meta$parameter['r.margin_sub.delta']
#   r.margin_delta <- meta$parameter['r.margin_delta']
#   r.margin_sec <- meta$parameter['r.margin_sec']
# 
# 
#   w.class <- meta$warm_info['class']
#   w.ap.h1.min <- meta$warm_info['ap.h1.min']
#   w.ap.h1.max <- meta$warm_info['ap.h1.max']
#   w.ap.delta.min <- meta$warm_info['ap.delta.min']
#   w.ap.delta.max <- meta$warm_info['ap.delta.max']
#   w.rp.delta.min <- meta$warm_info['rp.delta.min']
#   w.rp.delta.max <- meta$warm_info['rp.delta.max']
# 
#   w.type <- meta$warm_info['type']
#   w.med.t <- meta$warm_info['med.t']
#   w.usage.sec <- meta$warm_info['usage.sec']
#   w.usage.watt <- meta$warm_info['usage.watt']
# 
#   if (meta$general_info['warm_DB']) {
#     w.cand.DB <- data.frame(med.t = meta$warm_DB_time)
#   }
# 
#   if (meta$general_info['cook_DB']) {
# 
#     c.DB.total.lump <- meta$cook_DB_info['total.lump']
#     c.DB.valid.lump <- meta$cook_DB_info['valid.lump']
#     c.DB.ap.h2.min <- meta$cook_DB_info['ap.h2.min']
#     c.DB.ap.h2.max <- meta$cook_DB_info['ap.h2.max']
#     c.DB.rp.delta.min <- meta$cook_DB_info['rp.delta.min']
#     c.DB.rp.delta.max <- meta$cook_DB_info['rp.delta.max']
#     c.DB.med.t <- meta$cook_DB_info['med.t']
#     c.DB.med.fluc <- meta$cook_DB_info['med.fluc']
#     c.DB.usage.sec <- meta$cook_DB_info['usage.sec']
#     c.DB.usage.watt <- meta$cook_DB_info['usage.watt']
#     c.DB.usage.ap <- meta$cook_DB_info['usage.ap']
#   }
# 
#   if (meta$general_info['cook_blind']) {
# 
#     c.blind.ap.delta.med <- meta$cook_blind_info['ap.delta.med']
#     c.blind.rp.delta.med <- meta$cook_blind_info['rp.delta.med']
#     c.blind.time_gap <- meta$cook_blind_info['time_gap']
#     c.blind.fluc_num <- meta$cook_blind_info['fluc_num']
#     c.blind.on_sec <- meta$cook_blind_info['on_sec']
#   }
# 
#   # examine raw data
#   # data <- data[order(data$timestamp), ]
# 
#   # (1) mark warm mode
#   warm_mode_search_flag <- F
#   valid_warm_mode_flag <- F
# 
#   if (meta$general_info['warm_DB'] || meta$general_info['cook_blind']) {
# 
#     if (meta$general_info['cook_DB']) {
# 
#       if (c.DB.valid.lump != 0) warm_mode_search_flag <- T
# 
#     } else {
# 
#       warm_mode_search_flag <- T
#     }
#   }
# 
#   if ( warm_mode_search_flag ) {
# 
#     w.pattern.s <- DetectPattern_1Hz_new(data, position = "start", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
# 
#     if (w.class == 1){
#       chosen.w.pattern <- subset(w.pattern.s, h1 >= w.ap.h1.min & h1 <= w.ap.h1.max & delta >= w.ap.delta.min & delta <= w.ap.delta.max &
#                                    abs(sub.delta) <= thres_rp.delta)
#     } else if (w.class == 2){
#       chosen.w.pattern <- subset(w.pattern.s, h1 >= w.ap.h1.min & h1 <= w.ap.h1.max & sub.delta >= w.rp.delta.min & sub.delta <= w.rp.delta.max)
#     }
# 
#     if (nrow(chosen.w.pattern) < eff_size){
# 
#       print("There's no signal for warming mode.")
#     } else {
# 
#       chosen.w.pattern <- chosen.w.pattern[order(chosen.w.pattern$start.idx), ]
#       lump.idx <- c(0, which(diff(as.numeric(chosen.w.pattern$start.timestamp)) >= w.lump.gap_sec.max ), nrow(chosen.w.pattern))
#       chosen.lump <- which(diff(lump.idx) >= eff_size)
# 
#       if(length(chosen.lump)== 0){
# 
#         print("There's no lump for warming mode.")
#       } else {
# 
#         lump_log <- list()
# 
#         # examine each lump
#         for (l_idx in 1:length(chosen.lump)) {
# 
#           tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
#           tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
# 
#           lump_log[[length(lump_log)+1]] <- data.frame(chosen.w.pattern[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx)
#         }
# 
#         lump_log.total <- rbind.fill(lump_log)
# 
#         lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))) )
# 
#         if( w.type > 0){
# 
#           DB_result <- as.logical(as.numeric(tail(strsplit(DecToBin(w.type), "")[[1]], length(w.cand.DB$med.t))))
# 
#           # < another option >
#           # max_type_num <- as.numeric(head(strsplit(as.character(w.type), "")[[1]], 1))
#           # DB_result <- c(rep(TRUE, max_type_num), rep(FALSE, length(w.cand.DB$med.t) - max_type_num))
# 
#           lump_summary <- rbind.fill(mapply(function(c,m) {if (c) return(subset(lump_summary, med.t >= (m - w.lump.med_sec.margin) &
#                                                                                   med.t <= (m + w.lump.med_sec.margin ) ))
#             else NULL }, DB_result, w.cand.DB$med.t, SIMPLIFY = F))
#         } else {
# 
#           lump_summary <- subset(lump_summary, med.t >= (w.med.t - w.lump.med_sec.margin) &
#                                    med.t <= (w.med.t + w.lump.med_sec.margin) )
#         }
# 
#         if(nrow(lump_summary) == 0) {
# 
#           print("medium time is too chaotic in warming mode.")
#         } else {
# 
#           lump_log.total <- rbind.fill(lump_log[lump_summary$lump_idx])
#           valid_warm_mode_flag <- TRUE
# 
#           ### insert warming signals (for results)
#           new_timestamp <- floor_date(lump_log.total$start.timestamp, unit = "second")
# 
#           warming.start.idx <- which(answer.log$timestamp %in% new_timestamp)
# 
#           if(round(w.usage.sec) <= 0 || w.usage.watt <= 0){
# 
#             print("energy estimation failure on warming mode")
#           } else {
# 
#             warming.on.idx <- unique(unlist(lapply( warming.start.idx, function(x) x + 0:(round(w.usage.sec)-1) )))
# 
#             ###### avoid overflow of the idx
#             warming.on.idx <- warming.on.idx[warming.on.idx <= nrow(answer.log)]
# 
#             answer.log$p[warming.on.idx] <- w.usage.watt * usage.pwr.adjust
#             ### end
#             print("Warming mode is done.")
#           }
# 
#         }
#       }
#     }
# 
#   }
# 
#   # (2) mark cook mode
#   if (meta$general_info['cook_DB']) {
# 
#     c.pattern.e <- DetectPattern_1Hz_new(data, position = "end", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
# 
#     chosen.e.info.cooking <- subset( c.pattern.e, (c.DB.ap.h2.min-r.margin_delta) <= h2 & (c.DB.ap.h2.max+r.margin_delta) >= h2 & (c.DB.rp.delta.min - r.margin_sub.delta) <= sub.delta &
#                                        (c.DB.rp.delta.max + r.margin_sub.delta) >= sub.delta)
# 
#     lump.idx <- c(0, which(diff(as.numeric(chosen.e.info.cooking$start.timestamp)) >= c.lump.gap_sec.max), nrow(chosen.e.info.cooking))
#     chosen.lump <- which(diff(lump.idx) >= eff_size)
# 
#     if (length(chosen.lump) >0) {
# 
#       lump_log <- list()
# 
#       for (l_idx in 1:length(chosen.lump)) {
# 
#         tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
#         tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
# 
#         lump_log[[length(lump_log)+1]] <- data.frame(chosen.e.info.cooking[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx)
#       }
# 
#       lump_log.total <- rbind.fill(lump_log)
# 
#       # examine cooking candidates
# 
#       ### 1) refine lumps w/ median sig. period
# 
#       lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t )
# 
#       lump_summary <- subset(lump_summary, med.t >= (c.DB.med.t * (1 - c.lump.med_sec.margin) ) &
#                                med.t <= (c.DB.med.t * (1 + c.lump.med_sec.margin) ) )
# 
#       if (nrow(lump_summary) == 0) {
#         print('cook mode identification failure with DB')
#       } else {
# 
#         total_log <- rbind.fill(lump_log[lump_summary$lump_idx])
#         log_summary <- ddply(total_log, .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.med = median(start.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
#                              min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
#                              med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                              lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#         print(log_summary)
# 
#         # (old one)
#         #         CDB_duration <- mapply( function(s,e) c.time <- (as.numeric(s) + as.numeric(e))/2 , log_summary$lump.start, log_summary$lump.end, SIMPLIFY = T )
#         #         CDB_duration <- as.POSIXct(CDB_duration ,origin = "1970-01-01", tz = "Asia/Seoul")
# 
#         # (new one: avg. -> med.)
#         CDB_duration <- log_summary$lump.med
# 
#         # insert cooking signals (for results)
#         cooking.start.timestamp <- floor_date(CDB_duration, unit = "second")
#         cooking.start.idx <- which(answer.log$timestamp %in% cooking.start.timestamp)
#         cooking.on.idx <- unique(unlist(lapply( cooking.start.idx, function(x) x + -round(c.DB.usage.sec/2):(round(c.DB.usage.sec/2)-1) )))
#         ### avoid overflow of the idx
#         cooking.on.idx <- cooking.on.idx[(cooking.on.idx >= 1) & (cooking.on.idx <= nrow(answer.log)) ]
#         answer.log$p[cooking.on.idx] <- answer.log$p[cooking.on.idx] - c.DB.usage.ap
#         # end
# 
#       }
#     }
#   } else if (meta$general_info['cook_blind'] && valid_warm_mode_flag) {
# 
#     base_warm.summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
#                                min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
#                                med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                                lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#     ###### 1) choose cooking signal
# 
#     ######### 1-1) determine search range on each (warming) lump
# 
#     cooking.search_range <- pmin(as.numeric(base_warm.summary$lump.start) - as.numeric( c(min(data$timestamp), base_warm.summary$lump.end[-nrow(base_warm.summary)])))
#     cooking.eff_start <- c.lump_merge.max_hr *3600 < cooking.search_range
#     base_warm.summary <- data.frame(base_warm.summary, forward.search_time = cooking.search_range, eff_start = cooking.eff_start)
#     base_warm.idx <- c( which(base_warm.summary$eff_start == T), length(base_warm.summary$eff_start)+1)
# 
#     ###### 2) lump merging
#     ######### case 1
#     eff_start.idx <- base_warm.idx[ which(diff(base_warm.idx) == 1)]
# 
#     eff_warm.summary <- data.frame(forward_search = base_warm.summary$forward.search_time[eff_start.idx],
#                                    lump.start = base_warm.summary$lump.start[eff_start.idx],
#                                    lump.end = base_warm.summary$lump.end[eff_start.idx])
# 
#     ######### case 2
#     eff_start.idx <- base_warm.idx[ which(diff(base_warm.idx) >1)]
# 
#     eff_end.idx <- base_warm.idx[ which(diff(base_warm.idx) >1)+1] -1
# 
#     eff_warm.summary <-rbind(eff_warm.summary, data.frame(forward_search = base_warm.summary$forward.search_time[eff_start.idx],
#                                                           lump.start = base_warm.summary$lump.start[eff_start.idx],
#                                                           lump.end = base_warm.summary$lump.end[eff_end.idx]) )
# 
#     ######### ordering
#     eff_warm.summary <- eff_warm.summary[order(eff_warm.summary$lump.start),]
# 
#     ######### refine forward searching time
#     eff_warm.summary$forward_search <- pmin(eff_warm.summary$forward_search, c.blind.time_gap + r.margin_sec)
#     eff_warm.summary$backward_search.from.start <- pmin(-(c.blind.time_gap - r.margin_sec), as.numeric( c(eff_warm.summary$lump.start[-1]-dseconds(eff_warm.summary$forward_search[-1]),
#                                                                                                           max(data$timestamp))) - as.numeric(eff_warm.summary$lump.start) )
# 
#     cat("warming mode summary for blind detection: \n")
#     print(eff_warm.summary)
# 
#     ######### search cooking signal
# 
#     c.pattern.e <- DetectPattern_1Hz_new(data, position = "end", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
# 
#     e.pattern.cooking <- cookSigSearch_15Hz(data, end.pattern = c.pattern.e,
#                                             min_watt = c.blind.ap.delta.med-r.margin_delta, max_watt = c.blind.ap.delta.med+r.margin_delta,
#                                             min_t.sec = -Inf, max_t.sec = c.max_t.sec, min_fluc.num = 0)
# 
#     e.pattern.cooking <- subset(e.pattern.cooking, (sub.delta >= (c.blind.rp.delta.med - r.margin_sub.delta)) &
#                                   (sub.delta <= (c.blind.rp.delta.med + r.margin_sub.delta)))
# 
#     blind_cook.info <- mapply( function(s,f,b){ z <- subset(e.pattern.cooking, end.timestamp %within% interval(s-dseconds(f), s+dseconds(b) ))
#     if (nrow(z) != 0) return(z[which.max(z$on_sec),]) else return(z) },
#     eff_warm.summary$lump.start, eff_warm.summary$forward_search, eff_warm.summary$backward_search.from.start, SIMPLIFY = F)
# 
#     blind_cook.idx <- which(sapply(blind_cook.info, nrow) != 0)
# 
#     if (length(blind_cook.idx) == 0) {
# 
#       print("Cooking signal refinement has failed after adjusting the time gap between two modes.")
#     } else {
# 
#       ######### 4) yield information
#       blind_cook.info <- rbind.fill(blind_cook.info)
#       blind_cook.info <- blind_cook.info[duplicated(blind_cook.info$start.idx) == F,]
# 
#       cat("resultant cooking mode for blind detection: \n")
#       print(blind_cook.info)
# 
#       # insert cooking signals (for results)
#       cooking.end.timestamp <- floor_date(blind_cook.info$end.timestamp, unit = "second")
#       cooking.end.idx <- which(answer.log$timestamp %in% cooking.end.timestamp)
#       cooking.on.idx <- unique(unlist(lapply( cooking.end.idx, function(x) x - floor(c.blind.on_sec):0)))
#       ### avoid overflow of the idx
#       cooking.on.idx <- cooking.on.idx[(cooking.on.idx >= 1) & (cooking.on.idx <= nrow(answer.log)) ]
#       answer.log$p[cooking.on.idx] <- answer.log$p[cooking.on.idx] + c.blind.ap.delta.med
#       # end
#     }
# 
#   }
# 
#   return(answer.log)
# 
# }
# 
# energyEstimation_cookDB_15Hz <- function(data, end.pattern, max_t.sec){
# 
#   e.pattern.cooking <- end.pattern[c(1, nrow(end.pattern)), ]
# 
#   # (old one)
#   #   eff.energy1 <- rbind.fill(mapply( function(t,p,d) { z <- subset(data, timestamp %within% interval(t-dseconds(max_t.sec), t+dseconds(max_t.sec) ) & active_power >= p )
#   #                                                      num <- nrow(z)
#   #                                                      return(data.frame(row.num = num, energy = num*d))}, e.pattern.cooking$end.timestamp, e.pattern.cooking$ap.thres,
#   #                                                      abs(e.pattern.cooking$h2), SIMPLIFY = F))
# 
#   # (new one: complexity reduction)
#   eff.energy <- rbind.fill(mapply( function(t,p,d) { DataIdx <- sliceTimeIdx_Quick(data$timestamp, t-dseconds(max_t.sec), t+dseconds(max_t.sec))
#   pwr_val <- data$active_power[DataIdx]
#   num <- length( pwr_val[pwr_val >= p] )
#   return(data.frame(row.num = num, energy = num*d ))}, e.pattern.cooking$end.timestamp, e.pattern.cooking$ap.thres,
#   abs(e.pattern.cooking$h2), SIMPLIFY = F))
# 
#   e.pattern.cooking <- data.frame(e.pattern.cooking, cook.esti_on_sec = eff.energy$row.num/15, cook.esti_Wh = eff.energy$energy/15/3600)
# 
#   e.pattern.cooking <- e.pattern.cooking[which.max(e.pattern.cooking$cook.esti_on_sec),]
# 
#   return( data.frame(cook_esti_sec = e.pattern.cooking$cook.esti_on_sec, cook_esti_Wh = e.pattern.cooking$cook.esti_Wh, cook_esti_ap = e.pattern.cooking$h2) )
# }
# 
# generate.PatternScan.meta.ricecooker_15Hz <- function (data, eff_size = 5, eff_group_size = 40, periodicity = 0.1, thres_ap.h = 15, thres_rp.delta = 20){
# 
#   # options(scipen = 15)
# 
#   # parameters
#   search_iter <- 600
# 
#   # warming mode search
#   w.thres_ap.max <- 170
# 
#   w.med_sec.min <- 10
#   w.med_sec.max <- 240
# 
#   # w.lump (basic parameters)
#   w.lump.gap_sec.max <- 270
#   w.lump.med_sec.margin <- 1
# 
#   # w.lump (candidate DB)
#   w.cand.DB <- data.frame(med.t = c(16, 30, 32, 48, 64, 80, 96, 112, 128, 144, 160), #30 is from the brand 'Coochan'
#                           min.med.rate2 = c(1, 1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10 ) )
#   # magnitude DB (off the record): 1(106, 70, (91.5, 103.5), 42)
#   # magnitude DB (off the record): 2(40.5), 3(100), 4(60)
# 
#   DB_tolerance <- .02 # 2%
#   DB_rate_quantile <- .08 # 8%
# 
#   DB_timespan <- .1 # compared to maximum duration
# 
# 
#   # cooking mode search
#   cooking.min_watt <- 750
#   cooking.max_watt <- 1600
#   cooking.min_t.sec <- 7*60
#   cooking.max_t.sec <- 40*60
#   cooking.min_fluc.num <- 0
# 
#   # c.lump (basic parameters)
#   c.lump.gap_sec.max <- 3600
#   c.lump.med_sec.margin <- 0.0625 # 1/16 : In case of 16 sec, the margin is 1 sec.
# 
#   # cooking repetition DB
#   c.cand.DB <- data.frame(med.t = c(16, 32))
# 
#   CDB_tolerance <- .02 # 2%
# 
#   # connection between warming & cooking
#   wnc.forward.time_diff.max_hr <- 4
#   wnc.backward.time_diff.max_hr <- 0.5
#   lump_merge.max_hr <- 1/3
# 
# 
#   # cooking mode refine
#   margin_sub.delta <- 20
#   margin_delta <- 50 # watt
#   margin_sec <- 120*60 #second
# 
#   # cooking.med_watt <- 1200
# 
#   # warming & cooking signal search (in parallel)
#   w.pattern.s <- DetectPattern_1Hz_new(data, position = "start", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
# 
#   warm.info_low <- parallel_warmSig_search_15Hz(data, w.pattern.s, ap_h1_min = thres_ap.h, ap_h1_max = w.thres_ap.max, ap_delta_min = thres_ap.h, ap_delta_max = w.thres_ap.max,
#                                                 rp_delta_thres = thres_rp.delta, main.g.name = "h1", sub.g.name = "delta", periodicity, search_iter, w.med_sec.min, w.med_sec.max,
#                                                 eff_size, eff_group_size, w.lump.gap_sec.max, w.lump.med_sec.margin, w.cand.DB, DB_tolerance, DB_rate_quantile, DB_timespan)
# 
#   warm.info_high <- list()
#   if (!any(unlist(warm.info_low$type) > 0)){
# 
#     warm.info_high <- parallel_warmSig_search_15Hz(data, w.pattern.s, ap_h1_min = w.thres_ap.max, ap_h1_max = cooking.min_watt, ap_delta_min = -Inf, ap_delta_max = Inf,
#                                                    rp_delta_thres = Inf, main.g.name = "h1", sub.g.name = "sub.delta", periodicity, search_iter, w.med_sec.min, w.med_sec.max,
#                                                    eff_size, eff_group_size, w.lump.gap_sec.max, w.lump.med_sec.margin, w.cand.DB, DB_tolerance, DB_rate_quantile, DB_timespan)
#   }
# 
#   cook.info <- parallel_cookSig_search_15Hz(data, cooking.min_watt, cooking.max_watt, cooking.max_t.sec, cooking.min_fluc.num, periodicity,
#                                             search_iter, c.med_sec.min = w.med_sec.min, c.med_sec.max = w.med_sec.max, eff_size, c.lump.gap_sec.max,
#                                             c.lump.med_sec.margin, c.cand.DB, CDB_tolerance)
# 
#   valid.info <- sig.Orchestration_RC_15Hz(data, warm.info_low, warm.info_high, cook.info, w.cand.DB, c.cand.DB, DB_timespan,
#                                           cooking.forward_search = wnc.forward.time_diff.max_hr, cooking.backward_search = wnc.backward.time_diff.max_hr,
#                                           lump_merge = lump_merge.max_hr, cooking.min_sec = cooking.min_t.sec,
#                                           margin_sub.delta, margin_delta, margin_sec)
# 
#   if (length(valid.info) == 0) {
#     print("Valid signal search for a rice cooker has failed.")
#     return(list())
#   }
# 
#   # make parameter group to return
#   parameters <- c('eff_size' = eff_size, 'thres_ap.h' = thres_ap.h, 'thres_rp.delta' = thres_rp.delta, 'w.lump.gap_sec.max' = w.lump.gap_sec.max,
#                   'w.lump.med_sec.margin' = w.lump.med_sec.margin, 'c.max_t.sec' = cooking.max_t.sec, 'c.lump.gap_sec.max' = c.lump.gap_sec.max,
#                   'c.lump.med_sec.margin' = c.lump.med_sec.margin, 'c.lump_merge.max_hr' = lump_merge.max_hr, 'r.margin_sub.delta' = margin_sub.delta,
#                   'r.margin_delta' = margin_delta, 'r.margin_sec' = margin_sec)
# 
#   # JSON result
#   return(JSON_result_15Hz(data, information = valid.info, ap_thres = thres_ap.h, min_cook.time = cooking.min_t.sec, parameters,
#                           warm_DB_time = w.cand.DB$med.t))
# 
# }
# 
# 
# ###############################################################################################
# # Legacy Function (for 1Hz)
# # These are useless when 15Hz functions are completed
# ###############################################################################################
# 
# generate.PatternScan.meta.ricecooker_1Hz <- function (data, eff_size, periodicity, thres_ap.h, thres_rp.delta){
# 
#   # parameters
#   search_iter <- 1000
#   warming.min_minute <- 4.5 # examine each lump to check 'med.t'
#   warming.lump.min_gap.minute <- 30 # lumps within x min. are merged when searching cooking (Also, this value should be larger than or equal to 'cooking.max_t.minute' to avoid an overlapped power estimation.)
#   warming.sec_margin <- 3 # x second margin on median time
# 
#   wnc_max.hr.gap <- 4 # time interval between warming & cooking
# 
# 
#   cooking.min_watt <- 800
#   cooking.med_watt <- 1200
#   cooking.max_watt <- 1600
# 
# 
#   # cooking.min_Wh <- 30 # Wh
#   cooking.min_t.minute <- 4 # to choose the right candidate
#   cooking.min_fluctuation.num <- 0 # to choose the right candidate
# 
#   cooking.max_t.minute <- 30 # this value should be smaller than or equal to 'cooking.max_t.minute' to avoid an overlapped power estimation.
# 
#   cooking.cand_split.watt <- 400
# 
#   # cooking information refinement (based on median values)
#   r.sub.delta_margin <- 45
#   r.delta_margin <- 150 # watt
#   r.time.gap.minute_margin <- 60
#   r.cooking.minute_margin <- 8
# 
#   # examine raw data
#   data <- data[order(data$timestamp),]
#   #   time.diff <- diff(as.numeric(data$timestamp))
#   #   cat("strange time difference percentage in data:", length(time.diff[time.diff<0.9])/length(time.diff) * 100  , "\n")
# 
# 
#   s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
# 
#   s.pattern <- s.pattern[order(s.pattern$start.idx),]
# 
#   # medium-power signals for warming
#   s.pattern.warming1 <- subset(s.pattern, (h1 >= thres_ap.h) & (abs(sub.delta) <= thres_rp.delta) & !is.na(min.slope) ) # original set
#   s.pattern.warming2 <- subset(s.pattern, (h1 >= thres_ap.h) & is.na(min.slope) ) # added set
# 
#   s.pattern.warming <- rbind(s.pattern.warming1, s.pattern.warming2)
#   s.pattern.warming <- s.pattern.warming[order(s.pattern.warming$start.idx),]
# 
#   s.info <- list()
#   s.info <- DetectGroup_15Hz_1D(s.pattern.warming, resolution = search_iter, p_factor = periodicity, group_size = eff_size )
# 
#   json.candidate <- list()
#   json.result <- list()
#   str.t <- as.character(min(data$timestamp))
#   end.t <- as.character(max(data$timestamp))
# 
#   stored_lump.summary <- list()
#   stored_eff.cooking_info <- list()
#   stored_lump.log <- list()
#   stored_med_time <- list()
# 
#   if (length(s.info) == 0) {
#     print("Detection for ricecooker has failed.")
#     return(list())
#   } else {
# 
#     s.list <- s.info[[length(s.info)]]
#     s.list <- subset(s.list, med.t <= warming.min_minute*60)
# 
#     if (nrow(s.list) == 0) {
#       print("Detection for ricecooker has failed.")
#       return(list())
# 
#     } else {
#       print("groups for warming mode:")
#       print(s.list)
# 
#       # high-power signals for cooking (naive)
#       s.pattern.cooking <- cooking.candidate.search(data, start.pattern = s.pattern, min_watt = cooking.min_watt, max_watt = cooking.max_watt,
#                                                     min_t.minute = cooking.min_t.minute, max_t.minute = cooking.max_t.minute, min_fluc.num = cooking.min_fluctuation.num)
# 
#       # examine candidates for warming mode
#       for (w_idx in 1:nrow(s.list)) {
# 
#         print(w_idx)
#         chosen.s.info <- s.info[[as.numeric(row.names(s.list[w_idx,])) ]]
#         lump.idx <- c(0, which(diff(as.numeric(chosen.s.info$start.timestamp)) >= warming.min_minute*60 ), nrow(chosen.s.info))
#         chosen.lump <- which(diff(lump.idx) >= eff_size)
# 
#         if (length(chosen.lump) >0) {
# 
#           lump_log <- list()
# 
#           for (l_idx in 1:length(chosen.lump)) {
# 
#             tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
#             tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
# 
#             #             tmp.start.timestamp <- chosen.s.info$start.timestamp[tmp.start.idx]
#             #             tmp.end.timestamp <- chosen.s.info$end.timestamp[tmp.end.idx]
#             #             print(tmp.start.timestamp %--% tmp.end.timestamp)
# 
#             lump_log[[length(lump_log)+1]] <- data.frame(chosen.s.info[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx )
# 
#           }
# 
#           lump_log.total <- rbind.fill(lump_log)
# 
#           lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))) )
# 
#           tmp.med_time <- lump_summary$med.t[which.max(lump_summary$sum)]
# 
#           lump_summary <- subset(lump_summary, med.t >= (tmp.med_time - warming.sec_margin) &
#                                    med.t <= (tmp.med_time + warming.sec_margin) )
# 
#           lump_log.total <- rbind.fill(lump_log[lump_summary$lump_idx])
# 
#           lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
#                                 min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
#                                 med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                                 lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#           if (nrow(lump_summary) != 1) {
# 
#             lump_eff.searching.time <- pmin(wnc_max.hr.gap*3600, (as.numeric(lump_summary$lump.start) - as.numeric(c(min(data$timestamp) ,lump_summary$lump.end[-nrow(lump_summary)]) ) ) )
#           } else {
# 
#             lump_eff.searching.time <- pmin(wnc_max.hr.gap*3600, (as.numeric(lump_summary$lump.start) - as.numeric(min(data$timestamp) ) ) )
#           }
#           lump_eff.start <- warming.lump.min_gap.minute* 60 < lump_eff.searching.time
#           lump_summary <- data.frame(lump_summary, searching.time = lump_eff.searching.time, eff.start = lump_eff.start)
# 
#           print(lump_summary)
#           eff.lump_summary <- subset(lump_summary, eff.start == T)
# 
#           eff.cooking_info <- mapply(function(s,t) { z <- subset(s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#           if (nrow(z) != 0) return(z[which.max(z$est.Wh),]) else return(z) },
#           eff.lump_summary$lump.start, eff.lump_summary$searching.time, SIMPLIFY = F)
# 
#           cooking_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) != 0)
# 
#           if (length(cooking_eff.lump.idx) != 0) {
# 
#             eff.cooking_info <- rbind.fill(eff.cooking_info)
#             #             wnc_gap.minute <- (as.numeric(eff.lump_summary$lump.start[cooking_eff.lump.idx]) - as.numeric(eff.cooking_info[,'start.timestamp']) )/60
#             #             eff.cooking_info <- cbind(eff.cooking_info, time_gap = wnc_gap.minute)
#             stored_lump.log[[w_idx]] <- lump_log.total
#             stored_lump.summary[[w_idx]] <- lump_summary
#             stored_eff.cooking_info[[w_idx]] <- eff.cooking_info
#             stored_med_time[[w_idx]] <- tmp.med_time
# 
#             print(eff.lump_summary)
#             print(eff.cooking_info)
#           }
#         }
#       }
# 
#       total_eff.cooking_info <- rbind.fill(stored_eff.cooking_info)
# 
#       if (length(total_eff.cooking_info) == 0) {
# 
#         cat("There is no candidate for COOKING signal \n")
#         cat("If you are sure on the existence of rice cookers, please adjust the cooking signal power range set from ", cooking.min_watt ,"to ", cooking.max_watt ," watt. \n")
#         return(json.result)
#       } else {
# 
#         # remove duplicated candidate
#         watt.ordered <- total_eff.cooking_info[duplicated(total_eff.cooking_info$start.idx) == F,]
#         watt.ordered <- watt.ordered[order(watt.ordered$delta),]
#         watt_idx <- which(diff(watt.ordered$delta) >= cooking.cand_split.watt)
# 
#         if (length(watt_idx) != 0) {
#           cat("There are", length(watt_idx)+1, "groups for cooking candidates." , "\n")
#           cooking.watt_s.idx <- c(0, watt_idx) + 1
#           cooking.watt_e.idx <- c(watt_idx, nrow(watt.ordered))
# 
#           cooking.group_med.watt <- mapply(function(s,e) median(watt.ordered$delta[s:e]),
#                                            cooking.watt_s.idx, cooking.watt_e.idx)
# 
#           chosen_cooking.cand.idx <- which.min(abs(cooking.group_med.watt-cooking.med_watt))
#           chosen_cooking.cand <- watt.ordered[cooking.watt_s.idx[chosen_cooking.cand.idx]:cooking.watt_e.idx[chosen_cooking.cand.idx],]
# 
#         } else {
# 
#           print("There is no ambiguity in choosing cooking candidates.")
#           chosen_cooking.cand <- watt.ordered
#         }
# 
#         # refine cooking information
#         tmp_med.sub.delta <- median(chosen_cooking.cand$sub.delta)
#         tmp_med.delta <- median(chosen_cooking.cand$delta)
#         tmp_med.est.fluc.num <- median(chosen_cooking.cand$est.fluc.num)
# 
#         tmp_med.est.on.sec <- median(chosen_cooking.cand$est.on.sec)
#         tmp_min.est.on.sec <- min(chosen_cooking.cand$est.on.sec)
# 
#         tmp_s.pattern <- subset(s.pattern, sub.delta >= tmp_med.sub.delta - r.sub.delta_margin &
#                                   sub.delta <= tmp_med.sub.delta + r.sub.delta_margin )
# 
#         refined_s.pattern.cooking <- cooking.candidate.search(data, start.pattern = tmp_s.pattern, min_watt = (tmp_med.delta - r.delta_margin), max_watt = (tmp_med.delta + r.delta_margin),
#                                                               min_t.minute = 0,
#                                                               max_t.minute = min(cooking.max_t.minute, (tmp_med.est.on.sec/60) + r.cooking.minute_margin),
#                                                               min_fluc.num = 0)
# 
#         if (nrow(refined_s.pattern.cooking) == 0) {
# 
#           print("There are ambiguities in finding reasonable COOKING signals: please input more data so we can label the rice cooker in a reliable (i.e., statistical) way.")
#         } else {
# 
#           # recalculation of 'lump summary': step 1) w/ refined cooking s.pattern
#           for (w_idx in 1:length(stored_lump.summary)) {
#             if ( !is.null(stored_lump.summary[[w_idx]]) ) {
# 
#               lump_summary <- stored_lump.summary[[w_idx]]
#               chosen.s.info <- stored_lump.log[[w_idx]]
#               eff.lump_summary <- subset(lump_summary, eff.start == T)
# 
#               eff.cooking_info <- mapply(function(s,t) { z <- subset(refined_s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#               if (nrow(z) != 0) return(z[which.max(z$est.Wh),]) else return(z) },
#               eff.lump_summary$lump.start, eff.lump_summary$searching.time, SIMPLIFY = F)
# 
#               cooking_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) != 0)
# 
#               # recalculation of 'lump summary': step 2) w/ refined time gap between warming & cooking (search again)
#               if (length(cooking_eff.lump.idx) != 0) {
#                 eff.cooking_info <- rbind.fill(eff.cooking_info)
#                 wnc_gap.minute <- (as.numeric(eff.lump_summary$lump.start[cooking_eff.lump.idx]) - as.numeric(eff.cooking_info[,'start.timestamp']) )/60
#                 eff.cooking_info <- cbind(eff.cooking_info, time_gap = wnc_gap.minute)
#                 tmp_med.time_gap <- median(eff.cooking_info$time_gap)
# 
#                 if (nrow(lump_summary) != 1) {
# 
#                   lump_summary$searching.time <- pmin((tmp_med.time_gap + r.time.gap.minute_margin)*60, (as.numeric(lump_summary$lump.start) - as.numeric(c(min(data$timestamp) ,lump_summary$lump.end[-nrow(lump_summary)]) ) ) )
#                 } else {
# 
#                   lump_summary$searching.time <- pmin((tmp_med.time_gap + r.time.gap.minute_margin)*60, (as.numeric(lump_summary$lump.start) - as.numeric(min(data$timestamp) ) ) )
#                 }
# 
#                 lump_summary$eff.start <- min(cooking.max_t.minute*60, tmp_med.est.on.sec + r.cooking.minute_margin*60) < lump_summary$searching.time
#                 eff.lump_summary <- subset(lump_summary, eff.start == T)
# 
#                 eff.cooking_info <- mapply(function(s,t) { z <- subset(refined_s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#                 if (nrow(z) != 0) return(z[which.max(z$est.Wh),]) else return(z) },
#                 eff.lump_summary$lump.start, eff.lump_summary$searching.time, SIMPLIFY = F)
# 
#                 cooking_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) != 0)
#                 warming_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) == 0)
# 
#                 if (length(cooking_eff.lump.idx) != 0) {
# 
#                   eff.cooking_info <- rbind.fill(eff.cooking_info)
#                   wnc_gap.minute <- (as.numeric(eff.lump_summary$lump.start[cooking_eff.lump.idx]) - as.numeric(eff.cooking_info[,'start.timestamp']) )/60
#                   eff.cooking_info <- cbind(eff.cooking_info, time_gap = wnc_gap.minute)
# 
#                   print(eff.lump_summary)
#                   print(eff.cooking_info)
# 
#                   # refine cooking information
#                   tmp_med.sub.delta <- median(eff.cooking_info$sub.delta)
#                   tmp_med.delta <- median(eff.cooking_info$delta)
#                   tmp_med.est.fluc.num <- median(eff.cooking_info$est.fluc.num)
#                   tmp_med.time_gap <- median(eff.cooking_info$time_gap)
# 
#                   tmp_med.est.on.sec <- median(eff.cooking_info$est.on.sec)
#                   tmp_min.est.on.sec <- min(eff.cooking_info$est.on.sec)
# 
#                   # energy computation
#                   usage_info <- energyEstimation_ricecooker(data, pattern = chosen.s.info, threshold = thres_ap.h)
# 
#                   if (length(usage_info) != 0){
#                     warming.esti.sec <- usage_info[1]
#                     warming.esti.watt <- usage_info[2]
#                   } else {
#                     warming.esti.sec <- 0
#                     warming.esti.watt <- 0
#                   }
# 
#                   # computation for json information
# 
#                   if (cooking.min_t.minute*60 <= tmp_med.est.on.sec) {
# 
#                     json.candidate[[length(json.candidate)+1]] <- data.frame(cand_idx = w_idx, cooking_num = length(cooking_eff.lump.idx), warming_num = length(warming_eff.lump.idx),
#                                                                              total_sig.num = sum(lump_summary$sum),
#                                                                              ap.h1.min = min(chosen.s.info$h1), ap.h1.max = max(chosen.s.info$h1), ap.h1.med = median(chosen.s.info$h1),
#                                                                              rp.delta.med = median(chosen.s.info$sub.delta), min.t = median(lump_summary$min.t), med.t = stored_med_time[[w_idx]],
#                                                                              min.med.rate2 = median(lump_summary$min.med.rate2), usage_pwr_med = median(lump_summary$med.d),
#                                                                              usage_pwr_max = max(lump_summary$med.d), esti_usage.sec = warming.esti.sec, esti_usage.watt = warming.esti.watt,
#                                                                              c.rp.med_delta = tmp_med.sub.delta, c.ap.med_delta = tmp_med.delta, c.ap.med_fluc.num = tmp_med.est.fluc.num,
#                                                                              c.med_time.gap = tmp_med.time_gap, c.med_cooking.time = tmp_med.est.on.sec, c.min_cooking.time = tmp_min.est.on.sec)
# 
#                   } else {
# 
#                     print("Detected cooking information is not reasonable.")
#                   }
#                 }
#               }
#             }
#           }
#         }
#       }
# 
# 
#       #       for (w_idx in 1:nrow(s.list)) {
#       #
#       #
#       #             # refine cooking information
#       #             tmp_med.sub.delta <- median(eff.cooking_info$sub.delta)
#       #             tmp_med.delta <- median(eff.cooking_info$delta)
#       #             tmp_med.est.fluc.num <- median(eff.cooking_info$est.fluc.num)
#       #             tmp_med.time_gap <- median(eff.cooking_info$time_gap)
#       #
#       #             tmp_med.est.on.sec <- median(eff.cooking_info$est.on.sec)
#       #             tmp_min.est.on.sec <- min(eff.cooking_info$est.on.sec)
#       #
#       #             print(data.frame(length(cooking_eff.lump.idx), length(warming_eff.lump.idx), tmp_med.sub.delta,
#       #                              tmp_med.delta, tmp_med.est.on.sec, tmp_min.est.on.sec,
#       #                              tmp_med.est.fluc.num, tmp_med.time_gap))
#       #
#       #             tmp_s.pattern <- subset(s.pattern, sub.delta >= tmp_med.sub.delta - r.sub.delta_margin &
#       #                                       sub.delta <= tmp_med.sub.delta + r.sub.delta_margin )
#       #             refined_s.pattern.cooking <- cooking.candidate.search(data, start.pattern = tmp_s.pattern, min_watt = (tmp_med.delta - r.delta_margin), max_watt = (tmp_med.delta + r.delta_margin),
#       #                                                                   min_t.minute = 0,
#       #                                                                   max_t.minute = min(cooking.max_t.minute, (tmp_med.est.on.sec/60) + r.cooking.minute_margin),
#       #                                                                   min_fluc.num = 0)
#       #
#       #             if (nrow(refined_s.pattern.cooking) == 0) {
#       #               print("There are ambiguities in finding reasonable cooking signals.")
#       #
#       #             } else {
#       #
#       #               # recalculation of 'lump summary'
#       #               if (nrow(lump_summary) != 1) {
#       #
#       #                 lump_summary$searching.time <- pmin((tmp_med.time_gap + r.time.gap.minute_margin)*60, (as.numeric(lump_summary$lump.start) - as.numeric(c(min(data$timestamp) ,lump_summary$lump.end[-nrow(lump_summary)]) ) ) )
#       #               } else {
#       #
#       #                 lump_summary$searching.time <- pmin((tmp_med.time_gap + r.time.gap.minute_margin)*60, (as.numeric(lump_summary$lump.start) - as.numeric(min(data$timestamp) ) ) )
#       #               }
#       #
#       #               lump_summary$eff.start <- min(cooking.max_t.minute*60, tmp_med.est.on.sec + r.cooking.minute_margin*60) < lump_summary$searching.time
#       #               eff.lump_summary <- subset(lump_summary, eff.start == T)
#       #
#       #               eff.cooking_info <- mapply(function(s,t) { z <- subset(refined_s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#       #                                                          if (nrow(z) != 0) return(z[which.max(z$est.Wh),]) else return(z) },
#       #                                          eff.lump_summary$lump.start, eff.lump_summary$searching.time, SIMPLIFY = F)
#       #
#       #               cooking_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) != 0)
#       #               warming_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) == 0)
#       #
#       #               if (length(cooking_eff.lump.idx) != 0) {
#       #
#       #                 eff.cooking_info <- rbind.fill(eff.cooking_info)
#       #                 wnc_gap.minute <- (as.numeric(eff.lump_summary$lump.start[cooking_eff.lump.idx]) - as.numeric(eff.cooking_info[,'start.timestamp']) )/60
#       #                 eff.cooking_info <- cbind(eff.cooking_info, time_gap = wnc_gap.minute)
#       #
#       #                 print(eff.lump_summary)
#       #                 print(eff.cooking_info)
#       #
#       #                 # refine cooking information
#       #                 tmp_med.sub.delta <- median(eff.cooking_info$sub.delta)
#       #                 tmp_med.delta <- median(eff.cooking_info$delta)
#       #                 tmp_med.est.fluc.num <- median(eff.cooking_info$est.fluc.num)
#       #                 tmp_med.time_gap <- median(eff.cooking_info$time_gap)
#       #
#       #                 tmp_med.est.on.sec <- median(eff.cooking_info$est.on.sec)
#       #                 tmp_min.est.on.sec <- min(eff.cooking_info$est.on.sec)
#       #
#       #                 # computation for json information
#       #
#       #                 if (cooking.min_t.minute*60 <= tmp_med.est.on.sec) {
#       #
#       #                   json.candidate[[length(json.candidate)+1]] <- data.frame(cand_idx = w_idx, cooking_num = length(cooking_eff.lump.idx), warming_num = length(warming_eff.lump.idx),
#       #                                                                            total_sig.num = sum(lump_summary$sum),
#       #                                                                            ap.h1.min = min(chosen.s.info$h1), ap.h1.max = max(chosen.s.info$h1), ap.h1.med = median(chosen.s.info$h1),
#       #                                                                            rp.delta.med = median(chosen.s.info$sub.delta), min.t = median(lump_summary$min.t) ,med.t = tmp.med_time,
#       #                                                                            min.med.rate2 = median(lump_summary$min.med.rate2), usage_pwr_med = median(lump_summary$med.d),
#       #                                                                            usage_pwr_max = max(lump_summary$med.d),
#       #                                                                            c.rp.med_delta = tmp_med.sub.delta, c.ap.med_delta = tmp_med.delta, c.ap.med_fluc.num = tmp_med.est.fluc.num,
#       #                                                                            c.med_time.gap = tmp_med.time_gap, c.med_cooking.time = tmp_med.est.on.sec, c.min_cooking.time = tmp_min.est.on.sec)
#       #
#       #                 } else {
#       #
#       #                   print("Detected cooking information is not reasonable.")
#       #                 }
#       #
#       #
#       #
#       #               }
#       #
#       #             }
#       #
#       #           }
#       #
#       #         }
#       #       }
# 
# 
#       # make JSON meta file
#       if (length(json.candidate) != 0){
# 
#         json.candidate.list <- rbind.fill(json.candidate)
# 
#         print(json.candidate.list)
# 
#         max.cooking.num <- max(json.candidate.list$cooking_num)
# 
#         js <- subset(json.candidate.list, cooking_num == max.cooking.num)
# 
#         if ( nrow(js) > 1 ) {
#           js <- js[which.max(js$total_sig.num),]
#         }
# 
#         general.json.info <- c('cooking_num' = js$cooking_num, 'warming_num' = js$warming_num)
# 
#         warming.json.info <- c('eff_sample.num' = js$total_sig.num, 'ap.h1.min' = js$ap.h1.min, 'ap.h1.max' = js$ap.h1.max, 'ap.h1.med' = js$ap.h1.med,
#                                'rp.delta.med' = js$rp.delta.med, 'min.t' = js$min.t, 'med.t' = js$med.t, 'min.med.rate2' = js$min.med.rate2,
#                                'usage_pwr_med' = js$usage_pwr_med, 'usage_pwr_max' = js$usage_pwr_max, 'esti_usage.sec' = js$esti_usage.sec, 'esti_usage.watt' = js$esti_usage.watt)
# 
#         cooking.json.info <- c('c.rp.med_delta' = js$c.rp.med_delta, 'c.ap.med_delta' = js$c.ap.med_delta, 'c.ap.med_fluc.num' = js$c.ap.med_fluc.num,
#                                'c.med_time.gap' = js$c.med_time.gap, 'c.med_cooking.time' = js$c.med_cooking.time, 'c.min_cooking.time' = js$c.min_cooking.time)
# 
#         meta.parameters <- c('eff_size' = eff_size, 'periodicity' = periodicity, 'thres_ap.h' = thres_ap.h, 'thres_rp.delta' = thres_rp.delta,
#                              'search_iter' = search_iter, 'warming.min_minute' = warming.min_minute, 'warming.lump.min_gap.minute' = warming.lump.min_gap.minute,
#                              'warming.sec_margin' = warming.sec_margin, 'wnc_max.hr.gap' = wnc_max.hr.gap, 'cooking.min_watt' = cooking.min_watt,
#                              'cooking.max_watt' = cooking.max_watt, 'cooking.min_t.minute' = cooking.min_t.minute, 'cooking.min_fluctuation.num' = cooking.min_fluctuation.num,
#                              'cooking.max_t.minute' = cooking.max_t.minute, 'r.sub.delta_margin' = r.sub.delta_margin, 'r.delta_margin' = r.delta_margin,
#                              'r.time.gap.minute_margin' = r.time.gap.minute_margin, 'r.cooking.minute_margin' = r.cooking.minute_margin)
# 
#         json.result <- list(`meta-version` = 1, shape_type = "ricecooker_pattern_scan", general_info = general.json.info,
#                             warming_info = warming.json.info,
#                             cooking_info = cooking.json.info,
#                             parameters = meta.parameters,
#                             generation_info = list(data_used = list(start = str.t, end = end.t, sampling = 1),
#                                                    computed = as.character(Sys.time())) )
#       } else {
# 
#         print("There is no candidate for rice cooker.")
#       }
#     }
#   }
# 
#   return(json.result)
# }
# 
# meta2PatternScan.ricecooker_1Hz <- function (data, meta, usage.pwr.adjust = 1.5, compensation_min.minute = 4) {
# 
#   # parameters
#   str.t <- as.character(min(data$timestamp))
#   end.t <- as.character(max(data$timestamp))
#   answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"),
#                                             floor_date(max(data$timestamp), unit = "second"), 'secs'))
#   answer.log <- data.frame( answer.log, p= rep(0,nrow(answer.log)), q= rep(0,nrow(answer.log)))
# 
#   '%nin%'<- Negate('%in%')
# 
#   # meta$warming_info['']
#   ap.h1.min <- meta$warming_info['ap.h1.min']
#   ap.h1.max <- meta$warming_info['ap.h1.max']
#   meta$warming_info['ap.h1.med']
#   meta$warming_info['rp.delta.med']
#   min.t <- meta$warming_info['min.t']
#   w.med.t <- meta$warming_info['med.t']
#   meta$warming_info['min.med.rate2']
#   usage_pwr_med <- meta$warming_info['usage_pwr_med']
#   usage_pwr_max <- meta$warming_info['usage_pwr_max']
# 
#   esti_usage.sec <- meta$warming_info['esti_usage.sec']
#   esti_usage.watt <- meta$warming_info['esti_usage.watt']
# 
#   # meta$cooking_info['']
#   c.rp.med_delta <- meta$cooking_info['c.rp.med_delta']
#   c.ap.med_delta <- meta$cooking_info['c.ap.med_delta']
#   c.med_time.gap <- meta$cooking_info['c.med_time.gap']
#   c.med_cooking.time <- meta$cooking_info['c.med_cooking.time']
# 
# 
#   # meta$parameters['']
#   eff_size <- meta$parameters['eff_size']
#   meta$parameters['periodicity']
#   thres_ap.h <- meta$parameters['thres_ap.h']
#   thres_rp.delta <- meta$parameters['thres_rp.delta']
#   warming.min_minute <- meta$parameters['warming.min_minute']
#   #  warming.lump.min_gap.minute <- meta$parameters['warming.lump.min_gap.minute']
#   warming.sec_margin <- meta$parameters['warming.sec_margin']
#   #   meta$parameters['']
#   #   meta$parameters['']
#   r.sub.delta_margin <- meta$parameters['r.sub.delta_margin']
#   r.delta_margin <- meta$parameters['r.delta_margin']
#   cooking.max_t.minute <- meta$parameters['cooking.max_t.minute']
#   r.time.gap.minute_margin <- meta$parameters['r.time.gap.minute_margin']
#   #   meta$parameters['']
#   #   meta$parameters['']
# 
#   # examine raw data
#   data <- data[order(data$timestamp),]
# 
#   s.pattern <- DetectPattern_1Hz_new_modifiedforRicecooker(data, position = "start", main_type = "active", sub_type = "reactive")
#   s.pattern <- s.pattern[order(s.pattern$start.idx),]
# 
#   # medium-power signals for warming
#   s.pattern.warming1 <- subset(s.pattern, (h1 >= ap.h1.min) & (h1 <= ap.h1.max) & (abs(sub.delta) <= thres_rp.delta) & !is.na(min.slope) ) # original set
#   s.pattern.warming2 <- subset(s.pattern, (h1 >= ap.h1.min) & (h1 <= ap.h1.max) & is.na(min.slope) ) # added set
# 
#   s.pattern.warming <- rbind(s.pattern.warming1, s.pattern.warming2)
# 
#   if (nrow(s.pattern.warming) < eff_size){
# 
#     stop("There's no signal for warming mode.")
#   } else {
# 
#     chosen.s.info <- s.pattern.warming[order(s.pattern.warming$start.idx),]
#     lump.idx <- c(0, which(diff(as.numeric(chosen.s.info$start.timestamp)) >= warming.min_minute*60 ), nrow(chosen.s.info))
#     chosen.lump <- which(diff(lump.idx) >= eff_size)
# 
#     if (length(chosen.lump)  == 0) {
# 
#       stop("There's no lump for warming mode.")
#     } else {
# 
#       lump_log <- list()
# 
#       # examine each lump
#       for (l_idx in 1:length(chosen.lump)) {
# 
#         tmp.start.idx <- lump.idx[chosen.lump[l_idx]] +1
#         tmp.end.idx <- lump.idx[chosen.lump[l_idx]+1]
# 
#         lump_log[[length(lump_log)+1]] <- data.frame(chosen.s.info[tmp.start.idx:tmp.end.idx,], lump_idx = l_idx )
#       }
# 
#       lump_log.total <- rbind.fill(lump_log)
# 
#       lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), med.t = median(diff(as.numeric(start.timestamp))) )
# 
#       lump_summary <- subset(lump_summary, med.t >= (w.med.t - warming.sec_margin) &
#                                med.t <= (w.med.t + warming.sec_margin) )
# 
#       if (nrow(lump_summary) == 0) {
# 
#         stop("Not enough value for warming.sec_margin.")
#       } else {
# 
#         lump_log.total <- rbind.fill(lump_log[lump_summary$lump_idx])
# 
#         # insert warming signals (for results)
#         new_timestamp <- floor_date(lump_log.total$start.timestamp, unit = "second")
# 
#         warming.start.idx <- which(answer.log$timestamp %in% new_timestamp)
#         if (round(esti_usage.sec) != 0){
#           warming.on.idx <- unique(unlist(lapply( warming.start.idx, function(x) x + 0:(round(esti_usage.sec)-1) )))
#         } else {
#           print("Warming energy estimation is not valid.")
#           warming.on.idx <- unique(unlist(lapply( warming.start.idx, function(x) x + 0:floor(min.t))))
#         }
# 
#         ### avoid overflow of the idx
#         warming.on.idx <- warming.on.idx[warming.on.idx <= nrow(answer.log)]
# 
#         if (round(esti_usage.sec) != 0){
# 
#           warming.usage_pwr <- esti_usage.watt
#         } else {
# 
#           if (usage_pwr_med < 0) {
#             warming.usage_pwr <- usage_pwr_max
#           } else {
#             warming.usage_pwr <- (usage_pwr_med + usage_pwr_max)/2
#           }
#         }
# 
#         answer.log$p[warming.on.idx] <- warming.usage_pwr * usage.pwr.adjust
#         # end
# 
#         lump_summary <- ddply(lump_log.total, .(lump_idx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
#                               min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
#                               med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                               lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#         if (nrow(lump_summary) != 1) {
# 
#           lump_searching.time <- pmin((c.med_time.gap + r.time.gap.minute_margin)*60, (as.numeric(lump_summary$lump.start) - as.numeric(c(min(data$timestamp) ,lump_summary$lump.end[-nrow(lump_summary)]) ) ) )
#         } else {
# 
#           lump_searching.time <- pmin((c.med_time.gap + r.time.gap.minute_margin)*60, (as.numeric(lump_summary$lump.start) - as.numeric(min(data$timestamp) ) ) )
#         }
# 
#         lump_eff.start <- cooking.max_t.minute*60 < lump_searching.time
#         lump_summary <- data.frame(lump_summary, searching.time = lump_searching.time, eff.start = lump_eff.start)
# 
#         print(lump_summary)
# 
#         eff.lump_summary <- subset(lump_summary, eff.start == T)
# 
#         if (nrow(eff.lump_summary) == 0) {
# 
#           print("Not enough lumps for searching cooking signals")
#         } else {
# 
#           # cooking signal computation
#           cooking_s.pattern <- subset(s.pattern, sub.delta >= c.rp.med_delta - r.sub.delta_margin &
#                                         sub.delta <= c.rp.med_delta + r.sub.delta_margin )
#           s.pattern.cooking <- cooking.candidate.search(data, start.pattern = cooking_s.pattern, min_watt = (c.ap.med_delta - r.delta_margin), max_watt = (c.ap.med_delta + r.delta_margin),
#                                                         min_t.minute = 0,
#                                                         max_t.minute = cooking.max_t.minute,
#                                                         min_fluc.num = 0)
# 
#           eff.cooking_info <- mapply(function(s,t) { z <- subset(s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#           if (nrow(z) != 0) return(z[which.max(z$est.Wh),]) else return(z) },
#           eff.lump_summary$lump.start, eff.lump_summary$searching.time, SIMPLIFY = F)
# 
#           cooking_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) != 0)
#           warming_eff.lump.idx <- which(sapply(eff.cooking_info, nrow) == 0)
# 
#           print(eff.lump_summary)
# 
#           if (length(cooking_eff.lump.idx) != 0) {
# 
#             eff.cooking_info <- rbind.fill(eff.cooking_info)
#             wnc_gap.minute <- (as.numeric(eff.lump_summary$lump.start[cooking_eff.lump.idx]) - as.numeric(eff.cooking_info[,'start.timestamp']) )/60
#             eff.cooking_info <- cbind(eff.cooking_info, time_gap = wnc_gap.minute)
#             print(eff.cooking_info)
# 
#             # insert cooking signals (for results)
#             cooking.start.timestamp <- floor_date(eff.cooking_info$start.timestamp, unit = "second")
#             cooking.start.idx <- which(answer.log$timestamp %in% cooking.start.timestamp)
#             cooking.on.idx <- unique(unlist(lapply( cooking.start.idx, function(x) x + 0:floor(c.med_cooking.time))))
#             ### avoid overflow of the idx
#             cooking.on.idx <- cooking.on.idx[cooking.on.idx <= nrow(answer.log)]
#             answer.log$p[cooking.on.idx] <- c.ap.med_delta
#             # end
# 
#             # cooking signal compensation
#             cooking.lump_start.idx <- c(0, which(diff(as.numeric(s.pattern.cooking$start.timestamp)) > cooking.max_t.minute*60 )) + 1
#             cooking.lump.start.log <- s.pattern.cooking[cooking.lump_start.idx,]
# 
#             detected.start.idx <- mapply(function(s) { z <- subset(cooking.lump.start.log, start.timestamp %within% interval(s-dminutes(cooking.max_t.minute),s+dminutes(cooking.max_t.minute) ) )
#             return(z) }, eff.cooking_info$start.timestamp, SIMPLIFY = F)
#             detected.start.idx <- rbind.fill(detected.start.idx)
# 
#             print("compensated cooking mode:")
#             cooking_compensated.start.log <- cooking.lump.start.log[which(cooking.lump.start.log$start.idx %nin% detected.start.idx$start.idx),]
# 
#           } else {
#             # cooking signal compensation
#             cooking.lump_start.idx <- c(0, which(diff(as.numeric(s.pattern.cooking$start.timestamp)) > cooking.max_t.minute*60 )) + 1
#             cooking.lump.start.log <- s.pattern.cooking[cooking.lump_start.idx,]
# 
#             print("compensated cooking mode:")
#             cooking_compensated.start.log <- cooking.lump.start.log
#           }
# 
# 
# 
#           cooking_compensated.start.log <- subset(cooking_compensated.start.log, est.on.sec >= compensation_min.minute*60 )
# 
#           print(cooking_compensated.start.log)
# 
#           if (nrow(cooking_compensated.start.log) != 0) {
# 
#             # insert cooking signals (for results)
#             cooking.start.timestamp <- floor_date(cooking_compensated.start.log$start.timestamp, unit = "second")
#             cooking.start.idx <- which(answer.log$timestamp %in% cooking.start.timestamp)
#             cooking.on.idx <- unique(unlist(mapply( function(x,d) x + 0:floor(d), cooking.start.idx, pmin(cooking_compensated.start.log$est.on.sec,c.med_cooking.time))))
# 
# 
#             ### avoid overflow of the idx
#             cooking.on.idx <- cooking.on.idx[cooking.on.idx <= nrow(answer.log)]
#             answer.log$p[cooking.on.idx] <- c.ap.med_delta
#             # end
# 
#           }
# 
#         }
# 
#       }
# 
#     }
# 
#   }
# 
#   return(answer.log)
# 
# }
# 
# DetectGroup_15Hz_1D <- function(data, resolution, p_factor, group_size){
# 
#   logSummary <- list()
#   logDetail <- list()
#   ordered.h <- data
# 
#   r_idx <- 1
#   iter_idx <- 1
#   while( iter_idx <= resolution ){
# 
#     # print(nrow(ordered.delta))
#     # for debug
#     if ( nrow(ordered.h) +1 <= r_idx ){
#       print("Resolution increment for searcing groups has stopped..")
#       break
#     }
# 
#     # examine delta
#     ordered.h <- ordered.h[order(ordered.h$h1),]
# 
#     # examine h
#     ordered.h <- ordered.h[order(ordered.h$h1),]
# 
#     ### detect group idx (i.e., wall for each group)
#     h.diff <- data.frame( diff.h = diff(ordered.h$h1), ordered.idx = seq(1, nrow(ordered.h)-1) )
#     ordered.h.diff <- h.diff[order(h.diff$diff.h, decreasing = TRUE),]
#     wall.idx <- ordered.h.diff$ordered.idx[1:r_idx]
#     wall.idx <- wall.idx[order(wall.idx)]
# 
#     ### identify group 1
#     g1 <- cut( seq(1, nrow(ordered.h)), c(0, wall.idx, nrow(ordered.h)) ,labels= FALSE ) -1
#     ordered.h <- cbind(ordered.h, g1)
#     ordered.h <- ordered.h[order(ordered.h$start.idx), ]
# 
#     # summary for each group
# 
#     group.info <- ddply(ordered.h, .(g1), transform, tmp_g.sum = length(g1))
# 
#     ### exclude the info from isolated elements
#     tmp_small.group <- subset(group.info, tmp_g.sum < group_size)
#     r_idx <- r_idx -nrow( ddply(tmp_small.group, .(g1), summarize, tmp_sum = 'del') )
#     if(r_idx < 1 ){r_idx <- 1}
# 
#     ordered.h <- subset(group.info, tmp_g.sum >= group_size)
#     group.info <- ddply(ordered.h, .(g1), summarize, sum = length(g1), min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))),
#                         min.med.rate = min.t/med.t, min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),  med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                         lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#     ### effective group (conservative search by default)
#     if (nrow(group.info) > 0){
# 
#       eff_group.info <- subset(group.info, (min.med.rate2 >= p_factor) & !(g1 == 0) )
# 
#       if(nrow(eff_group.info) >0){
# 
#         logSummary[[length(logSummary)+1]] <- eff_group.info
# 
#         for(g_idx in 1:nrow(eff_group.info) ){
# 
#           logDetail[[length(logDetail) +1]] <- ordered.h[(ordered.h$g1 == eff_group.info$g1[g_idx]), ]
#           logDetail[[length(logDetail)]]$g1 <- NULL
#           logDetail[[length(logDetail)]]$tmp_g.sum <- NULL
# 
#           ordered.h <- ordered.h[!(ordered.h$g1 == eff_group.info$g1[g_idx]) , ]
#         }
#       }
# 
#       # index increment for while loop
#       r_idx <- r_idx -nrow(eff_group.info) +1
#       if(r_idx < 1 ){r_idx <- 1}
#     } else {
# 
#       r_idx <- r_idx + 1
#     }
# 
#     iter_idx <- iter_idx + 1
#     ordered.h$g1 <- NULL
#     ordered.h$tmp_g.sum <- NULL
#   }
# 
#   # return information
#   if (length(logSummary) != 0) {
#     logDetail[[length(logDetail) +1]] <- rbind.fill(logSummary)
#     return(logDetail)
#   } else {
#     return(list())
#   }
# 
# }
# 
# DetectGroup_1Hz_onedim <- function(data, resolution, p_factor, group_size){
# 
#   logSummary <- list()
#   logDetail <- list()
#   ordered.delta <- data
# 
#   r_idx <- 1
#   iter_idx <- 1
#   while( iter_idx <= resolution ){
# 
#     # print(nrow(ordered.delta))
#     # for debug
#     if ( nrow(ordered.delta) +1 <= r_idx ){
#       print("Resolution has been reduced..")
#       break
#     }
# 
#     # examine delta
#     ordered.delta <- ordered.delta[order(ordered.delta$delta),]
# 
#     ### detect group idx (i.e., wall for each group)
#     delta.diff <- data.frame( diff.delta = diff(ordered.delta$delta), ordered.idx = seq(1, nrow(ordered.delta)-1) )
#     ordered.delta.diff <- delta.diff[order(delta.diff$diff.delta, decreasing = TRUE),]
#     wall.idx <- ordered.delta.diff$ordered.idx[1:r_idx]
#     wall.idx <- wall.idx[order(wall.idx)]
# 
#     ### identify group 1
#     g1 <- cut( seq(1, nrow(ordered.delta)), c(0, wall.idx, nrow(ordered.delta)) ,labels= FALSE ) -1
#     ordered.delta <- cbind(ordered.delta, g1)
#     ordered.delta <- ordered.delta[order(ordered.delta$start.idx), ]
#     # summary for each group
# 
#     group.info <- ddply(ordered.delta, .(g1), transform, tmp_g.sum = length(g1))
# 
#     ### exclude the info from isolated elements
#     tmp_small.group <- subset(group.info, tmp_g.sum < group_size)
#     r_idx <- r_idx -nrow( ddply(tmp_small.group, .(g1), summarize, tmp_sum = 'del') )
#     if(r_idx < 1 ){r_idx <- 1}
# 
#     ordered.delta <- subset(group.info, tmp_g.sum >= group_size)
#     group.info <- ddply(ordered.delta, .(g1), summarize, sum = length(g1), min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))),
#                         min.med.rate = min.t/med.t, min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.d = median(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
#                         lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
# 
#     ### effective group (conservative search by default)
#     if (nrow(group.info) > 0){
# 
#       eff_group.info <- subset(group.info, (min.med.rate >= p_factor) & !(g1 == 0) )
# 
#       if(nrow(eff_group.info) >0){
# 
#         logSummary[[length(logSummary)+1]] <- eff_group.info
# 
#         for(g_idx in 1:nrow(eff_group.info) ){
# 
#           logDetail[[length(logDetail) +1]] <- ordered.delta[(ordered.delta$g1 == eff_group.info$g1[g_idx]), ]
# 
#           ordered.delta <- ordered.delta[!(ordered.delta$g1 == eff_group.info$g1[g_idx]) , ]
#         }
#       }
# 
#       # index increment for while loop
#       r_idx <- r_idx -nrow(eff_group.info) +1
#       if(r_idx < 1 ){r_idx <- 1}
#     } else {
# 
#       r_idx <- r_idx + 1
#     }
# 
#     iter_idx <- iter_idx + 1
#     ordered.delta$g1 <- NULL
#     ordered.delta$tmp_g.sum <- NULL
#   }
# 
#   # return information
#   if (length(logSummary) != 0) {
#     logDetail[[length(logDetail) +1]] <- rbind.fill(logSummary)
#     return(logDetail)
#   } else {
#     return(list())
#   }
# 
# }
# 
# #============================================================================================
# # no restriction on delta and sub.delta !!!
# DetectPattern_1Hz_new_modifiedforRicecooker <- function ( data, position = c("start", "end")
#                                                           , main_type = c("reactive", "active")
#                                                           , sub_type = c("active", "reactive")
#                                                           , c_factor = 0
#                                                           , debug.mode = F ){
# 
#   mainLog <- switch( main_type,
#                      'reactive' = data.frame( timestamp = data$timestamp, value = data$reactive_power),
#                      'active' = data.frame( timestamp = data$timestamp, value = data$active_power) )
#   subLog <- switch( sub_type,
#                     'reactive' = data.frame( timestamp = data$timestamp, value = data$reactive_power),
#                     'active' = data.frame( timestamp = data$timestamp, value = data$active_power) )
# 
#   # for debug
#   if ( nrow(mainLog) != nrow(subLog) ){
#     cat("Log mismatch between active and reactive power!", nrow(mainLog), " and ", nrow(subLog), "\n")
#     stop("Error")
#   }
# 
#   slope <- get_speed_new(mainLog$timestamp, mainLog$value)
#   state <- scan_state_new(slope)
#   united.state <- stateCompression( state )
# 
#   edgeType <- match.arg(position)
#   if ( edgeType == "start"){
# 
#     rising.edge.pattern <- list( c(1,2,3,4), c(1,3,4), c(1,2,3,1), c(1,3,1) )
#     resultant1 <- ldply( rising.edge.pattern,
#                          function( pattern ){ patternIdx <- compressedPatternMatching_RC( united.state, pattern )
#                          get_patternFeature_RC( mainLog, subLog, slope, patternIdx )} )
# 
#     rising.edge.pattern <- list(  c(1,2), c(1,3) )
#     resultant2 <- ldply( rising.edge.pattern,
#                          function( pattern ){ patternIdx <- compressedPatternMatching_RC( united.state, pattern )
#                          get_patternFeature_RC( mainLog, subLog, slope, patternIdx )} )
#     resultant2 <- subset( resultant2, !(start.idx %in% resultant1$start.idx))
#     if( nrow(resultant2) > 0 ) resultant1 <- rbind.fill( resultant1, resultant2 )
# 
#     resultant <- resultant1
#     if( is.null(resultant) ) stop('There is no proper rising edge')
# 
#     if(debug.mode){
#       par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
#       plot(x=data$timestamp, y=data$active_power, ann=FALSE, xaxt="n", type='l')
#       abline( v=resultant$start.timestamp, col='red' )
#       plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
#       title(paste('Step1 : rising edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)),
#             outer=TRUE)
#       mtext("Timestamp", 1, 1, outer=TRUE)
#       mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
#       par(mfrow=c(1,1))
#     }
#     return( resultant )
# 
#   }else if (edgeType == "end"){
# 
#     falling.edge.pattern <- list( c(3,4), c(3,1), c(3,2) )
#     resultant <- ldply( falling.edge.pattern,
#                         function( pattern ){ patternIdx <- compressedPatternMatching_RC( united.state, pattern )
#                         get_patternFeature_RC( mainLog, subLog, slope, patternIdx, c_factor )} )
#     if( is.null(resultant) ) stop('There is no proper falling edge')
# 
#     if(debug.mode){
#       par(mfrow=c(2,1), oma=c(4, 4, 4, 2.5), mar=rep(.1, 4), cex=1, las=1)
#       plot(x=data$timestamp, y=data$active_power, ann=FALSE, xaxt="n", type='l')
#       abline( v=resultant$start.timestamp, col='blue' )
#       plot(x=data$timestamp, y=data$reactive_power, ann=FALSE, xaxt="n", type='l')
#       title(paste('Step1 : falling edge detection 결과\n', min(data$timestamp), '--', max(data$timestamp)),
#             outer=TRUE)
#       mtext("Timestamp", 1, 1, outer=TRUE)
#       mtext("Reactive / Active power", 2, 3, outer=TRUE, las=0)
#       par(mfrow=c(1,1))
#     }
# 
#     resultant[,sapply(resultant, is.numeric)] <- abs(resultant[,sapply(resultant, is.numeric)])
#     return( resultant )
#   }
# 
# }
# 
# compressedPatternMatching_RC <- function( compressedState, targetPattern ){
# 
#   pattern.length <- length(targetPattern)
# 
#   if( targetPattern[1] %in% c(1,2) ){ # rising edge
# 
#     if( any( targetPattern %in% c(3,4) & shift(targetPattern,-1) %in% c(1,2) ) ){
#       # power가 증가 -> 감소 -> 증가로 정의되면 마지막 증가는 패턴에서 제외
#       pattern.length <- max(which(targetPattern %in% c(3,4) & shift(targetPattern,-1) %in% c(1,2) ))
#     }
# 
#     matchedPattern.str <- grepl.pattern( compressedState$val, targetPattern )
#     matchedPattern.end <- matchedPattern.str + pattern.length - 1
#     maxSlope.pos <- matchedPattern.str + which( targetPattern == 1 & shift(targetPattern,-1) %in% c(2,3)) - 1
#     minSlope.pos <- matchedPattern.str + which( targetPattern == 3 & shift(targetPattern,-1) %in% c(4,1)) - 1
#     peak.pos <- matchedPattern.str + which( targetPattern %in% c(1,2) & shift(targetPattern,-1) %in% c(3,4)) - 1
# 
#     result <- data.frame(str = compressedState$str[ matchedPattern.str ],
#                          end = compressedState$end[ matchedPattern.end ])
# 
#     if( length(maxSlope.pos) != 0 ) result <- cbind( result, data.frame(maxSlope = compressedState$end[ maxSlope.pos ]))
#     if( length(minSlope.pos) != 0 ) result <- cbind( result, data.frame(minSlope = compressedState$end[ minSlope.pos ]))
#     if( length(peak.pos) == 0 ) peak.pos <- matchedPattern.end
# 
#     result <- cbind( result, data.frame(peak = compressedState$end[ peak.pos ] + 1))
#     return( result )
# 
#   }else if( targetPattern[1] %in% c(3,4) ){ # falling edge
# 
#     matchedPattern.str <- grepl.pattern( compressedState$val, targetPattern )
#     matchedPattern.end <- matchedPattern.str + pattern.length - 1
#     minSlope.pos <- matchedPattern.str + which( targetPattern == 3 & shift(targetPattern,-1) %in% c(4,1,2)) - 1
#     peak.pos <- matchedPattern.str
# 
#     return( data.frame(str = compressedState$str[ matchedPattern.str ],
#                        end = compressedState$end[ matchedPattern.end ],
#                        minSlope = compressedState$end[ minSlope.pos ],
#                        peak = compressedState$str[ peak.pos ]))
#   }
# }
# 
# #============================================================================================
# # h1 : peak - left (main)
# # h2 : rite - peak (main)
# # delta : rite - left = h1 + h2 (main)
# # sub.delta : rite - left (sub)
# get_patternFeature_RC <- function(data, sub.data, speed, patternIndex, c_factor, c_flag = F){
# 
#   timestamp <- data$timestamp
#   value <- data$value
#   sub.value <- sub.data$value
# 
#   patternFeature <- patternIndex[,c('str','end')]
#   patternFeature$len <- patternIndex$end - patternIndex$str + 1
#   patternFeature$str.t <- timestamp[patternIndex$str ]
#   patternFeature$end.t <- timestamp[patternIndex$end+1]
# 
#   if('maxSlope' %in% names(patternIndex)) patternFeature$max.slope <- speed[ patternIndex$maxSlope ]
#   if('minSlope' %in% names(patternIndex)) patternFeature$min.slope <- speed[ patternIndex$minSlope ]
# 
#   if('peak' %in% names(patternIndex)){
#     peak.value.set <- value[ patternIndex$peak ]
#     patternFeature$h1 <- peak.value.set - value[ patternIndex$str ]
#     patternFeature$h2 <- value[ patternIndex$end+1 ] - peak.value.set
#   }
# 
# 
#   patternFeature$sub.delta <- sub.value[patternIndex$end+1] - sub.value[patternIndex$str]
#   patternFeature$delta <- value[patternIndex$end+1] - value[patternIndex$str]
# 
#   if( !missing(c_factor) && c_factor >= 1 ){
# 
#     patternFeature$org_h2 <- patternFeature$h2
#     patternFeature$org_sub.delta <- patternFeature$sub.delta
# 
#     patternFeature <- subset( patternFeature, (str > c_factor*len) & (end < length(value)-c_factor*len) )
#     if( nrow(patternFeature) == 0 ) return(NULL)
# 
#     if (c_flag == T) {
# 
#       consistency_flag <- rep(1, nrow(patternFeature))
#       consistency_flag[ which(patternFeature$org_h2 > 0 | patternFeature$org_sub.delta > 0) ] <- 0
# 
#     }
# 
#     for(c_idx in 1:c_factor) {
# 
#       left.str <- patternFeature$str - patternFeature$len * c_idx
#       left.end <- patternFeature$str - patternFeature$len * (c_idx-1) - 1
#       rite.str <- patternFeature$end + patternFeature$len * (c_idx-1) + 2
#       rite.end <- patternFeature$end + patternFeature$len * c_idx + 1
# 
#       tmp_sign_delta <- mapply( function(ls,le,rs,re) mean(value[rs:re]) - mean(value[ls:le]),
#                                 left.str, left.end, rite.str, rite.end )
# 
#       tmp_sign_sub.delta <- mapply( function(ls,le,rs,re) mean(sub.value[rs:re]) - mean(sub.value[ls:le]),
#                                     left.str, left.end, rite.str, rite.end)
# 
#       if (c_flag == T) {
#         consistency_flag[ which(tmp_sign_delta > 0 | tmp_sign_sub.delta > 0) ] <- 0
#       }
# 
#       patternFeature$h2 <- patternFeature$h2 + tmp_sign_delta
#       patternFeature$delta <- patternFeature$delta + tmp_sign_delta
#       patternFeature$sub.delta <- patternFeature$sub.delta + tmp_sign_sub.delta
# 
#     }
# 
#     if (c_flag == T) {
#       patternFeature <- patternFeature[which(consistency_flag == 1),]
#     }
# 
#     patternFeature$h2 <- patternFeature$h2 / (c_factor+1)
#     patternFeature$delta <- patternFeature$delta / (c_factor+1)
#     patternFeature$sub.delta <- patternFeature$sub.delta / (c_factor+1)
#   }
#   names(patternFeature) <- mapvalues( names(patternFeature), warn_missing=F,
#                                       c('str','end','len','str.t','end.t'),
#                                       c('start.idx','end.idx','pattern.size',
#                                         'start.timestamp','end.timestamp') )
# 
#   return(patternFeature)
# }
# 
# cooking.candidate.search <- function(data, start.pattern, min_watt, max_watt, min_t.minute, max_t.minute, min_fluc.num) {
# 
#   s.pattern.cooking <- subset(start.pattern, delta >= min_watt & delta <= max_watt)
# 
#   if ( nrow(s.pattern.cooking) == 0 ) {
#     print("There is no candidate for cooking signals")
#     return(s.pattern.cooking)
#   } else {
# 
#     s.pattern.cooking_actual.center.ap <- apply(s.pattern.cooking, 1, function(x) ((data$active_power[as.numeric(x['start.idx'])] +
#                                                                                       data$active_power[as.numeric(x['end.idx'])]   )/2) )
# 
#     s.pattern.cooking <- data.frame(s.pattern.cooking, center.ap = s.pattern.cooking_actual.center.ap)
# 
#     #         s.pattern.cooking_eff.time.period <- apply(s.pattern.cooking, 1, function(x) { z <- subset(data[as.numeric(x['end.idx']):(as.numeric(x['end.idx'])+cooking.max_t.minute*60),],
#     #                                                                                                    active_power < as.numeric(x['center.ap']))
#     #                                                                                        if (nrow(z) != 0) return(min(z$timestamp)) else return(data$timestamp[(as.numeric(x['end.idx'])+cooking.max_t.minute*60)]) })
# 
#     s.pattern.cooking_eff.energy <- data.matrix( mapply(function(t,p,d) {z <- subset(data, timestamp %within% interval(t, t+dminutes(max_t.minute) ) &
#                                                                                        active_power > p )
#     num <- nrow(z)
#     return(c(num, num*d ))}, s.pattern.cooking$end.timestamp, s.pattern.cooking$center.ap, s.pattern.cooking$delta ) )
# 
#     s.pattern.cooking_fluc.num <- mapply( function(e) {z <- subset(s.pattern.cooking, start.timestamp %within% interval(e, e+dminutes(max_t.minute) ) )
#     return(nrow(z))}, s.pattern.cooking$end.timestamp)
# 
#     s.pattern.cooking <- data.frame(s.pattern.cooking, # time.on.sec = as.numeric(s.pattern.cooking_eff.time.period) - as.numeric(s.pattern.cooking$end.timestamp),
#                                     est.on.sec = s.pattern.cooking_eff.energy[1,], est.Wh = s.pattern.cooking_eff.energy[2,]/3600, est.fluc.num = s.pattern.cooking_fluc.num)
# 
#     s.pattern.cooking <- subset(s.pattern.cooking, (est.on.sec >= min_t.minute*60) & (est.fluc.num >= min_fluc.num))
# 
#     #     s.pattern.cooking_eff.fluc.num <- apply(s.pattern.cooking_naive, 1, function(x) { z <- subset(s.pattern.cooking_naive, start.timestamp.numeric >= as.numeric(x['start.timestamp.numeric']) &
#     #                                                                                                     start.timestamp.numeric <= (as.numeric(x['start.timestamp.numeric'])+1800))
#     #
#     #                                                                                       return(nrow(z)) })
# 
#     #           apply(s.pattern.cooking, 1, function(x) { z <- subset(data[as.numeric(x['end.idx']):(as.numeric(x['end.idx'])+cooking.max_t.minute*60),],
#     #                                                                                               active_power > as.numeric(x['center.ap']))
#     #                                                                                   return(nrow(z))} )
#     #
#     #         cooking.max.watt <- t(data.matrix( mapply(function(s,t) { z <- subset(s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#     #                                                                   if (nrow(z) != 0) return(z[which.max(z$delta),]) else return(z) },
#     #                                                                   lump_summary$lump.start, lump_summary$searching.time) ))
# 
#     return(s.pattern.cooking)
# 
#   }
# 
# }
# 
# cooking.candidate.search_15Hz <- function(data, start.pattern, min_watt, max_watt, min_t.minute, max_t.minute, min_fluc.num, Hertz = 15) {
# 
#   s.pattern.cooking <- subset(start.pattern, delta >= min_watt & delta <= max_watt)
# 
#   if ( nrow(s.pattern.cooking) == 0 ) {
#     print("There is no candidate for cooking signals")
#     return(s.pattern.cooking)
#   } else {
# 
#     s.pattern.cooking_actual.center.ap <- apply(s.pattern.cooking, 1, function(x) ((data$active_power[as.numeric(x['start.idx'])] +
#                                                                                       data$active_power[as.numeric(x['end.idx'])]   )/2) )
# 
#     s.pattern.cooking <- data.frame(s.pattern.cooking, center.ap = s.pattern.cooking_actual.center.ap)
# 
#     #         s.pattern.cooking_eff.time.period <- apply(s.pattern.cooking, 1, function(x) { z <- subset(data[as.numeric(x['end.idx']):(as.numeric(x['end.idx'])+cooking.max_t.minute*60),],
#     #                                                                                                    active_power < as.numeric(x['center.ap']))
#     #                                                                                        if (nrow(z) != 0) return(min(z$timestamp)) else return(data$timestamp[(as.numeric(x['end.idx'])+cooking.max_t.minute*60)]) })
# 
#     s.pattern.cooking_eff.energy <- data.matrix( mapply(function(t,p,d) {z <- subset(data, timestamp %within% interval(t, t+dminutes(max_t.minute) ) &
#                                                                                        active_power > p )
#     num <- nrow(z) / Hertz
#     return(c(num, num*d ))}, s.pattern.cooking$end.timestamp, s.pattern.cooking$center.ap, s.pattern.cooking$delta ) )
# 
#     s.pattern.cooking_fluc.num <- mapply( function(e) {z <- subset(s.pattern.cooking, start.timestamp %within% interval(e, e+dminutes(max_t.minute) ) )
#     return(nrow(z))}, s.pattern.cooking$end.timestamp)
# 
#     s.pattern.cooking <- data.frame(s.pattern.cooking, # time.on.sec = as.numeric(s.pattern.cooking_eff.time.period) - as.numeric(s.pattern.cooking$end.timestamp),
#                                     est.on.sec = s.pattern.cooking_eff.energy[1,], est.Wh = s.pattern.cooking_eff.energy[2,]/3600, est.fluc.num = s.pattern.cooking_fluc.num)
# 
#     s.pattern.cooking <- subset(s.pattern.cooking, (est.on.sec >= min_t.minute*60) & (est.fluc.num >= min_fluc.num))
# 
#     #     s.pattern.cooking_eff.fluc.num <- apply(s.pattern.cooking_naive, 1, function(x) { z <- subset(s.pattern.cooking_naive, start.timestamp.numeric >= as.numeric(x['start.timestamp.numeric']) &
#     #                                                                                                     start.timestamp.numeric <= (as.numeric(x['start.timestamp.numeric'])+1800))
#     #
#     #                                                                                       return(nrow(z)) })
# 
#     #           apply(s.pattern.cooking, 1, function(x) { z <- subset(data[as.numeric(x['end.idx']):(as.numeric(x['end.idx'])+cooking.max_t.minute*60),],
#     #                                                                                               active_power > as.numeric(x['center.ap']))
#     #                                                                                   return(nrow(z))} )
#     #
#     #         cooking.max.watt <- t(data.matrix( mapply(function(s,t) { z <- subset(s.pattern.cooking, start.timestamp %within% interval(s-dseconds(t),s))
#     #                                                                   if (nrow(z) != 0) return(z[which.max(z$delta),]) else return(z) },
#     #                                                                   lump_summary$lump.start, lump_summary$searching.time) ))
# 
#     return(s.pattern.cooking)
# 
#   }
# 
# }
# 
# 
# energyEstimation_ricecooker <- function(data, pattern, threshold){
# 
#   o_pattern <- pattern[ order(pattern$start.idx), ]
# 
#   slot_length <- pmax(0, o_pattern$start.idx[-1] - o_pattern$end.idx[-nrow(o_pattern)])
# 
#   slot.idx_list <- mapply( function(s,eVec){c(s,eVec)}, as.list( o_pattern$start.idx[-nrow(o_pattern)] ),
#                            mapply( function(x,l) (x+c(0:l)), o_pattern$end.idx[-nrow(o_pattern)], slot_length) )
# 
#   chosen.slot.info <- lapply(slot.idx_list, function(x){ signal.diff <- data$active_power[x[-1]] - data$active_power[x[1]]
#   if( any(abs(signal.diff) < threshold, na.rm=T)){
#     loc <- min( which( abs(signal.diff) < threshold ))
#     len <- x[-1][loc] - x[1] # +1 -1
#     W <- sum(pmax(0,data$active_power[x[1]:x[-1][loc]] - data$active_power[x[1]] ))/len
#     return(data.frame(length = len, watt = W))
#   }else return(NULL) } )
# 
#   if (length(chosen.slot.info[!sapply(chosen.slot.info, is.null )]) != 0) {
#     chosen.slot.info <- rbind.fill(chosen.slot.info)
#     return(c(median(chosen.slot.info$length), median(chosen.slot.info$watt)))
# 
#   } else {
#     print("energy estimation for warming mode has failed.")
#     return(c())
# 
#   }
# 
# }
# 


