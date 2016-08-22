###############################################################################################
# Functions for JPN
###############################################################################################
sigSearch_majorCase <- function(patterns, min_watt, max_watt, group_periodicity, group_quantProb, group_searchIter, group_medSec_min,
                                group_medSec_max, effective_size, lump_timeGapThres, timeTolerance, DB_majorSig, USE.DB = TRUE){
  
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
  
  if(!USE.DB) return(major_group_info)
  
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
    } else group.info <- data.frame()
    
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


sigOrchestration_riceCooker <- function(data, Hz, answer.log, major_SigInfo, major_DBInfo, thresTime, range_energyEsti, min_usageTime, min_apThres, annihilation_cookNumThres, annihilation_pattern,
                                        annihilation_spanSec, annihilation_apMinThres, annihilation_apMaxThres, annihilation_searchIter, annihilation_effSize, annihilation_medSec_min, 
                                        annihilation_medSec_max, annihilation_thresProb, annihilation_maxApGap, annihilation_coverageProb, lump_timeGapThres, annihilation_particleLumpSize){
  
  if( length(major_SigInfo) == 0){
    print("detection failure: no major signal")
    return(list())
  }
  
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  
  # possible properties of signals for the orchestration
  # 1) flag in DB 1) time duration # 1) active power # 1) sum of signatures # 1) number of lumps
  
  # step 1] classify the signals: short(0) or long(1)
  sigClass <- lapply(major_SigInfo$candSummary, function(summaryInfo){
    if( max(summaryInfo$lump.duration) >= thresTime) return(1) else return(0)
  })
  
  # step 2] process short(0) signals to identify the cooking mode
  ID_cookSig <- 0
  cookSummary <- data.frame()
  # cookEnergy <- data.frame(startTime = 0, endTime = 0,  SecON = 0, whON = 0, apON = 0, minSigNum = 0) 
  if(!any(sigClass == 0)) {
    print('no cooking candidates in major signals')
    
  } else {
    print('using proper short signals for cooking candidates')
    shortCand <- which(sigClass == 0)
    
    # details of the energy consumption
    energyEsti <- lapply(shortCand, function(sigIdx){
      lumpSummary <- major_SigInfo$candSummary[[sigIdx]]
      lumpLog <- major_SigInfo$candPattern[[sigIdx]] 
      
      estiOutput <- ldply(seq(nrow(lumpSummary)), function(rowIdx){
        tmp_lumpIdx <- lumpSummary$lumpIdx[rowIdx]
        tmp_LumpPattern <- lumpLog %>% filter(lumpIdx == tmp_lumpIdx)
        return(energyEstimation_caseCook(data, Hz, tmp_LumpPattern, range_energyEsti, min_usageTime))
      }) # %>% ldply(., .id= NULL)
      
      return(cbind(estiOutput, sigNum = rep(lumpSummary$sum, each = 2)))
    })
    
    # examine proper metric
    validMetric <- ldply(energyEsti, function(energyData){
      
      properUsage <- energyData %>% filter(SecON >= min_usageTime)
      if(nrow(properUsage) == 0) return(data.frame(sigNum_org = 0,sigNum_proper = 0,nRow_org = 0,nRow_proper = 0) ) else{
        return(data.frame(sigNum_org = min(energyData$sigNum), sigNum_proper = min(properUsage$sigNum), LumpNum_org = nrow(energyData),
                          LumpNum_proper = nrow(properUsage), sigNum_max = max(properUsage$sigNum), 
                          sigNum_min = min(properUsage$sigNum), secON_median = median(properUsage$SecON), secON_sd = sd(properUsage$SecON), whON_median = median(properUsage$whON), 
                          nRow_org = nrow(energyData)/2, nRow_proper = nrow(properUsage)/2) )
      }
    })
    
    # candMetric <- validMetric %>% mutate(chosenMetric = sigNum_proper * nRow_proper) %>% .$chosenMetric 
    candMetric <- validMetric %>% mutate(chosenMetric = sigNum_proper) %>% .$chosenMetric 
    
    # print("ricecooker sigNum:")    
    # print(max(candMetric)) # min. signal numer among the proper usages: 10 to 51
    
    if(all(candMetric == 0)){
      print('no cooking candidates with proper energy usage')
      
    } else {
      tmp_ID <- which.max(candMetric)
      proper_lumpNum <- validMetric$nRow_proper[tmp_ID]
      ID_cookSig <- shortCand[tmp_ID]
      cookEnergy <- energyEsti[[tmp_ID]] %>% filter(SecON >= min_usageTime) %>% 
        summarise(startTime = median(startTime), endTime = median(endTime), SecON = median(SecON), whON = median(whON), apON = median(apON), minSigNum = min(sigNum))
      
      # to split the cook signals into two categories (later in consumption output)
      cookSummary <- major_SigInfo$candSummary[[ID_cookSig]]
      cookActLog <- major_SigInfo$candPattern[[ID_cookSig]]
    }
  }
  
  # step 3] process long(1) signals + residual short(0) signals to detect warming mode
  ### 3-1] choose candidates
  ID_warmSig <- 0
  flag_warmType <- 0
  if(!any(sigClass == 1)) {
    # try 1] prototype
    # considerCand <- which(sigClass == 0)
    # considerCand <- considerCand[considerCand != ID_cookSig]
    # 
    # if(length(considerCand) == 0){ # in case of no possible candidate
    #   print('no warming candidates in major signals')
    #   longCand <- NULL
    # } else { 
    #   print('using residual signals for warming candidates')
    #   longCand <- considerCand
    # }
    
    # try 2] conservative approach, i.e., w/o reusing shortCand
    # print('no warming candidates in major signals')
    # longCand <- NULL
    
    # try 3] general signal investigation  
    # function name would be 'broadSearch_warmSig()'
    if( (ID_cookSig != 0) & (nrow(cookSummary) >= annihilation_cookNumThres)){
      
      longCand <- NULL
      print("annihilation start")
      
      ### reduce pattern number
      annihilation_pattern <- annihilation_pattern %>% filter(h2 < -annihilation_apMinThres, h2 > -annihilation_apMaxThres) 
      
      ### obtatin particle
      tmp_period <- cookSummary %>% mutate(startTime = lump.end, endTime = lump.end + dseconds(annihilation_spanSec) ) %>% select(startTime, endTime)
      particle_info <- mapply( function(sTime,eTime, Idx) return(annihilation_pattern %>% filter( end.timestamp %within% (sTime %--% eTime) ) %>% 
                                                                   rowwise() %>% mutate(particlePosition = as.numeric(end.timestamp) - as.numeric(sTime), lumpIdx = Idx ) ), 
                               tmp_period$startTime, tmp_period$endTime, seq(nrow(tmp_period)), SIMPLIFY = FALSE)
      
      ### obtain anti-particle
      tmp_period <- cookSummary %>% mutate(endTime = lump.start, startTime = lump.start - dseconds(annihilation_spanSec)) %>% select(startTime, endTime)
      antiParticle_info <- mapply( function(sTime,eTime, Idx) return(annihilation_pattern %>% filter( end.timestamp %within% (sTime %--% eTime) ) %>% 
                                                                       rowwise() %>% mutate(particlePosition = as.numeric(end.timestamp) - as.numeric(eTime), lumpIdx = Idx) ), 
                                   tmp_period$startTime, tmp_period$endTime, seq(nrow(tmp_period)), SIMPLIFY = FALSE)
      
      ### merge
      entireParticle <- bind_rows(particle_info, antiParticle_info) %>% arrange(start.timestamp)
      annihilationInfo <- detectGroup_pairAnnihilation(data = entireParticle, resolution = annihilation_searchIter, main.g.name = "delta", sub.g.name = "sub.delta", property.name = "particlePosition", 
                                                       a_factor = annihilation_thresProb, AP_maxGap = annihilation_maxApGap, med.time_min = annihilation_medSec_min, med.time_max = annihilation_medSec_max, group_size = annihilation_effSize)
      
      if(length(annihilationInfo) != 0){
        
        print("groups with annihilation:")
        group_list <- annihilationInfo[[length(annihilationInfo)]] %>% arrange(desc(sum))
        print(group_list)
        refinedAnnihilationLog <- chooseRefinedCand(annihilationInfo, entireCookLump = nrow(cookSummary), entirePattern = annihilation_pattern, 
                                                    annihilation_particleLumpSize, annihilation_coverageProb)
        
        if(nrow(refinedAnnihilationLog) != 0) {
          
          # remove cook signal
          compare_cookLog <- cookActLog %>% ungroup() %>% select(start.idx, end.idx, start.timestamp, end.timestamp, h2, sub.delta)
          compare_refinedLog <- refinedAnnihilationLog %>% select(start.idx, end.idx, start.timestamp, end.timestamp, h2, sub.delta)
          
          refined_warmLog <- dplyr::setdiff(compare_refinedLog, compare_cookLog)
          
          # consider lump
          if(nrow(refined_warmLog) != 0) {
            
            # split the chosen rows into several lumps depending on time gap
            splitInterval <- cumsum( c(TRUE, difftime( tail(refined_warmLog$start.timestamp,-1),
                                                       head(refined_warmLog$start.timestamp,-1), units='secs') >= lump_timeGapThres))
            splitLumps <- bind_cols(refined_warmLog, data.frame(lumpIdx = splitInterval))
            # splitLumps <- split.data.frame( chosen_rows, splitInterval) %>% bind_rows(.id = 'lumpIdx')
            
            # summarize information
            lumpSummary <- ddply(splitLumps, .(lumpIdx), summarize, sum = length(start.timestamp))
            # lumpSummary <- ddply(splitLumps, .(lumpIdx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
            #                      min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), max.t = max(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t)
            
            lumpSummary <- lumpSummary %>% filter(sum >= annihilation_particleLumpSize)
            refined_warmLog <- splitLumps %>% filter( lumpIdx %in% lumpSummary$lumpIdx)
            
            if(nrow(refined_warmLog) != 0){
              ID_warmSig <- length(sigClass) + 1
              flag_warmType <- 1
              warmLumpNum <- nrow(lumpSummary)
              # energy estimation (w/ original pattern log)
              warmEnergy_Info <- energyEstimation_caseWarm(data, Hz, pattern = refined_warmLog, thres = min_apThres)
            }
          }
        }
      }
      
    } else {
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
    warmLumpNum <- nrow(major_SigInfo$candSummary[[ID_warmSig]])
    # energy estimation (w/ original pattern log)
    warmEnergy_Info <- energyEstimation_caseWarm(data, Hz, pattern = major_SigInfo$candPattern[[ID_warmSig]], thres = min_apThres)
  }
  
  # step 4] create valid meta information
  generalInfo <- c("Cook" =  (ID_cookSig > 0), "Warm" = (ID_warmSig > 0))
  
  if( !(generalInfo[["Cook"]] | generalInfo[["Warm"]]) ){
    print("detection failure: no valid signal")
    return(list())
  }
  
  result <- list(meta_version = "0.1.0",
                 shape_type = "ricecooker_pattern_scan",
                 generation_info = list(
                   data_used = list(
                     start = str.t, 
                     end = end.t, 
                     sampling = 10
                   ),
                   computed = as.character(Sys.time())
                 ),
                 general_info = generalInfo)
  
  if (generalInfo[["Cook"]]){
    
    cookInfo <- c("lumps" = nrow(cookSummary),
                  "properLumps" = proper_lumpNum,
                  "min_ap.h2" = min(cookActLog$h2),
                  "max_ap.h2" = max(cookActLog$h2),
                  "min_rp.delta" = min(cookActLog$sub.delta),
                  "max_rp.delta" = max(cookActLog$sub.delta),
                  "DB_TimeSec" = major_SigInfo$candType[[ID_cookSig]],
                  "minFlucs_properLump" = cookEnergy$minSigNum,
                  "usage_startTime" = cookEnergy$startTime,
                  "usage_endTime" = cookEnergy$endTime,
                  "usage_secON" = cookEnergy$SecON,
                  "usage_whON" = cookEnergy$whON,
                  "usage_apON" = cookEnergy$apON
    )
    
    result$cook_DB_info <- cookInfo
  }
  
  if (generalInfo[["Warm"]]){
    
    if(flag_warmType == 0) DB_TimeSec <- major_SigInfo$candType[[ID_warmSig]] else DB_TimeSec <- 0
    
    warmInfo <- c("type" = flag_warmType, # 0: longCand // 1: annihilation
                  "samples" = nrow(refined_warmLog),
                  "lumps" = warmLumpNum,
                  "min_ap.h2" = min(refined_warmLog$h2),
                  "max_ap.h2" = max(refined_warmLog$h2),
                  "min_rp.delta" = min(refined_warmLog$sub.delta),
                  "max_rp.delta" = max(refined_warmLog$sub.delta),
                  "DB_TimeSec" = DB_TimeSec,
                  "usage_secON" = warmEnergy_Info$sec,
                  "usage_watt" = warmEnergy_Info$watt
    )
    
    result$warm_info <- warmInfo
  }
  
  return(result)  
  
}

energyEstimation_caseCook <- function(data, Hz, end.pattern, max_t.sec, thres_usageSec){
  
  # start and end point of a lump
  end.pattern <- end.pattern %>% mutate( ap.thres = (data$active_power[start.idx] + data$active_power[end.idx]) / 2 )
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

detectGroup_pairAnnihilation <- function(data, resolution, main.g.name, sub.g.name, property.name, a_factor, AP_maxGap, med.time_min, med.time_max, group_size){
  
  findGaps <- function(x,n){
    x <- sort(x)
    x.diff <- data.frame( val = diff(x), idx = 1:(length(x)-1) )
    wall.idx <- x.diff$idx[ order( x.diff$val, decreasing=T ) ][1:n] # choose first N gaps
    result <- sapply( wall.idx, function(i) mean(x[i+c(0,1)]))
    return( sort(result) )
  }
  
  summarize_particleInfo <- function( groupInfo, colName){
    
    tsDiff <- diff(as.numeric(sort(groupInfo$start.timestamp)))
    colParticle <- groupInfo[,colName]
    data.frame( 
      'sum' = length(groupInfo$start.timestamp), 
      'min.t'  = min(tsDiff),
      'med.t'  = median(tsDiff),
      'min.h2' = min(groupInfo$h2),
      'max.h2' = max(groupInfo$h2),
      'lost.sig.num' = sum( pmax( round(tsDiff / median(tsDiff)) -1 ,0)),
      'anti.particle' = nrow(colParticle %>% filter(. < 0)),
      'particle' = nrow(colParticle %>% filter(. >= 0)),
      'lump.coverage' = length(unique(groupInfo$lumpIdx))
    )
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
        do( summarize_particleInfo(., colName = property.name) ) %>%
        mutate( min.med.rate  = min.t/med.t, 
                particle.prob = particle/sum)
      
      group.info <- merge( group.info, o.list %>% 
                             summarise_each_( funs(median), c('h1','delta','sub.delta') )) %>%
        dplyr::rename( med.h1 = h1, med.d = delta, med.sub.d = sub.delta ) %>%
        mutate( lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
    } else group.info <- data.frame()
    
    ### effective group (conservative search by default)
    if (nrow(group.info) > 0){
      
      eff_group.info <- group.info %>%
        filter( med.t >= med.time_min ) %>%
        filter( med.t <= med.time_max ) %>%
        filter( abs(min.h2 - max.h2) <= AP_maxGap ) %>%
        filter( particle.prob >= a_factor ) %>%
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

chooseRefinedCand <- function(annihilationInfo, entireCookLump, entirePattern, annihilation_particleLumpSize, annihilation_coverageProb){
  
  antiParticle_reduction <- function(df, summary, group){
    # step 1] adjust h2
    if(df$delta < summary$med.d) refined_group <- group %>% filter(delta > df$delta) else refined_group <- group %>% filter(delta < df$delta)
    # step 2] adjust sub delta
    if(df$sub.delta < summary$med.sub.d) refined_group <- refined_group %>% filter(sub.delta > df$sub.delta) else refined_group <- refined_group %>% filter(sub.delta < df$sub.delta)
    
    return(refined_group)
  }
  
  compute_removingCost <- function(df, summary, group){
    
    refined_group <- antiParticle_reduction(df, summary, group)
    
    # step 3] yield result
    refinedParticleNum <- refined_group %>% filter(particlePosition >= 0) %>% nrow(.)
    removedParticleNum <- summary$particle - refinedParticleNum
    refinedAntiParticleNum <- refined_group %>% filter(particlePosition < 0) %>% nrow(.)
    removedAntiParticleNum <- summary$anti.particle - refinedAntiParticleNum
    
    return(data.frame(removedParticle = removedParticleNum, removedAntiParticle = removedAntiParticleNum, removingCost = removedParticleNum/removedAntiParticleNum))
  }
  
  eff_lumpCoverage <- sapply(seq(length(annihilationInfo)-1), function(x) { return(annihilationInfo[[x]] %>% filter(particlePosition >= 0) %>% group_by(lumpIdx) %>% 
                                                                                     summarise(particleSum = n()) %>% filter(particleSum >= annihilation_particleLumpSize) %>% nrow(.)) })
  
  list_validGroup <- bind_cols(annihilationInfo[[length(annihilationInfo)]], data.frame(orderNum = seq( length(annihilationInfo)-1 ), eff.lump.coverage = eff_lumpCoverage )) %>%  
    filter(eff.lump.coverage >= (entireCookLump* annihilation_coverageProb))
  
  if(nrow(list_validGroup) == 0){
    print("no valid group from annihilation")
    return(list_validGroup)
  }
  
  # step 1] remove anti particles by assessing cost
  improvedResult <- lapply(seq(nrow(list_validGroup)), function(candIdx){
    
    chosen_summary <- list_validGroup[candIdx, ]
    chosen_group <- annihilationInfo[[ chosen_summary$orderNum ]]
    
    chosen_antiParticle <- chosen_group %>% filter(particlePosition < 0)
    if(nrow(chosen_antiParticle) >= annihilation_particleLumpSize){ # not high quality information
      tmp_antiParticle_chkMetric <- chosen_antiParticle %>% rowwise() %>% do( compute_removingCost(df = ., summary = chosen_summary, group = chosen_group)) %>%
        bind_cols(chosen_antiParticle, .) %>% arrange(removingCost)
      antiParticle_chkMetric <- tmp_antiParticle_chkMetric %>% filter(removedAntiParticle >= (nrow(chosen_antiParticle)-annihilation_particleLumpSize) ) 
      
      if(nrow(antiParticle_chkMetric) != 0){
        antiParticle_chkMetric <- antiParticle_chkMetric[1,]
      } else{
        maxRemoveNum <- max(tmp_antiParticle_chkMetric$removedAntiParticle)
        antiParticle_chkMetric <- tmp_antiParticle_chkMetric[tmp_antiParticle_chkMetric$removedAntiParticle == maxRemoveNum, ] %>% arrange(removingCost) %>% .[1,]
      }
      cat("anti particle elimination at ", candIdx, "\n")
      print(antiParticle_chkMetric)
      return(antiParticle_reduction(df = antiParticle_chkMetric, summary = chosen_summary, group = chosen_group))
    }
    return(chosen_group)
  })
  
  # step 2] choose one list
  chosenResult <- sapply(improvedResult, function(x) return(nrow(x))) %>% which.max(.)
  
  # step 3] particle extension
  chosenInput <- improvedResult[[chosenResult]] %>% .[!duplicated(.$start.idx), ] %>% arrange(start.idx)
  extendedPattern <- entirePattern %>% filter(delta >= min(chosenInput$delta), delta <= max(chosenInput$delta), 
                                              sub.delta >= min(chosenInput$sub.delta), sub.delta <= max(chosenInput$sub.delta) )
  
  return(extendedPattern)
  
}


forceRiceCooker_jp <- function (data, prior_meta = NULL){
  
  # apply algorithm to each channel
  ### channel 1
  x.df <- data %>% filter(channel == 1) %>% arrange(timestamp)
  meta1 <- riceCookerCore_jp(data = x.df, prior_meta = prior_meta)
  if(!is.null(meta1)) meta1$channel <- 1
  
  ### channel 2
  x.df <- data %>% filter(channel == 2) %>% arrange(timestamp)
  meta2 <- riceCookerCore_jp(data = x.df, prior_meta = prior_meta)
  if(!is.null(meta2)) meta2$channel <- 2
  
  # compare metas
  if(is.null(meta1) && is.null(meta2)){
    print("no rice cooker at home: reuse prior meta")
    return(prior_meta)
  } else if(is.null(meta1) && !is.null(meta2)){
    return(meta2)
  } else if(!is.null(meta1) && is.null(meta2)){
    return(meta1)
  } else {
    print("meta orchestration")
    # step 1] examine cook signal
    if(meta1$general[["Cook"]] && !meta2$general[["Cook"]]){
      return(meta1)
    } else if(!meta1$general[["Cook"]] && meta2$general[["Cook"]]){
      return(meta2)
    } else {
      # step 2] examine warm signal
      if(meta1$general[["Warm"]] && !meta2$general[["Warm"]]){
        return(meta1)
      } else if(!meta1$general[["Warm"]] && meta2$general[["Warm"]]){
        return(meta2)
      } else {
        # step 3] examine cook lumps
        if(meta1$general[["Cook"]]){ # cook signals exist
          
          if(meta1$cook[['minFlucs_properLump']] >= meta2$cook[['minFlucs_properLump']]){ # meta2$cook[["properLumps"]]
            return(meta1)
          } else return(meta2)
        } else { # warm signals exist
          if(meta1$warm[["lumps"]] >= meta2$warm[["lumps"]]){
            return(meta1)
          } else return(meta2)
        }
      }
    }
  }
}


predict.forceRiceCooker_jp <- function(object, meta, apMargin = 50, rpMargin = 20, usage_powerAdjust = 2) {
  data <- object
  
  # options(scipen = 15)
  answer.log <- data.frame( timestamp = seq(floor_date(min(data$timestamp), unit = "second"),
                                            floor_date(max(data$timestamp), unit = "second"), 'secs'))
  answer.log <- data.frame( answer.log, p= rep(0,nrow(answer.log)), q= rep(0,nrow(answer.log)))
  cookCnt <- 0
  cookLog_forAnnihilation <- NULL
  
  if(length(meta) == 0) return(list(usage = answer.log, confidence = 0, count = cookCnt, samplingRate = 1.))
  
  data <- data %>% filter(channel == meta$channel) %>% arrange(timestamp)
  
  # parameters
  existCook <- meta$general_info[["Cook"]]
  existWarm <- meta$general_info[["Warm"]]
  metaHz <- meta$parameter[["Hz"]]
  timeTolerance <- meta$parameter[["DB_toleranceSec"]]
  major_lumpGap <- meta$parameter[["major_lumpGap"]]
  eff_size <- meta$parameter[["eff_size"]]
  
  #########################################################################
  ### DetectPattern_1Hz_new:: cowork with SJLEE
  e_pattern <- DetectPattern_1Hz_new(data, position = "end", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
  #########################################################################
  
  
  # support detection
  if(existCook || existWarm){
    
    periodicity <- meta$parameter[["periodicity"]]
    search_iter <- meta$parameter[["search_iter"]]
    quantileProb <- meta$parameter[["quantileProb"]]
    major_min_watt <- meta$parameter[["major_min_watt"]]  
    major_max_watt <- meta$parameter[["major_max_watt"]]
    major_medSec_min <- meta$parameter[["major_medSec_min"]]
    major_medSec_max <- meta$parameter[["major_medSec_max"]]
    DB_majorSig <- data.frame(time = meta$DBTime)
    
    sigGroupInfo <- sigSearch_majorCase( patterns = e_pattern, min_watt = major_min_watt, max_watt = major_max_watt, group_periodicity = periodicity, group_quantProb = quantileProb, 
                                         group_searchIter = search_iter, group_medSec_min = major_medSec_min, group_medSec_max = major_medSec_max, effective_size = eff_size, 
                                         lump_timeGapThres = major_lumpGap, timeTolerance, DB_majorSig, USE.DB = FALSE)
  }
  
  # cook signal search
  if(existCook){
    cook_DBTime <- meta$cook_DB_info[["DB_TimeSec"]]
    cook_minTimeON <- meta$parameter[["min_cookOn"]]
    cook_usageStartTime <- meta$cook_DB_info[["usage_startTime"]]
    cook_usageEndTime <- meta$cook_DB_info[["usage_endTime"]]
    cook_usageWhON <- meta$cook_DB_info[["usage_whON"]]
    # cook_usageSecON <- meta$cook_DB_info[["usage_secON"]]
    # cook_usageApON <- meta$cook_DB_info[["usage_apON"]]
    
    # set reduce
    cookCandLog <- e_pattern %>% filter( h2 >= meta$cook_DB_info[["min_ap.h2"]]-apMargin, h2 <= meta$cook_DB_info[["max_ap.h2"]]+apMargin, 
                                         sub.delta >= meta$cook_DB_info[["min_rp.delta"]]-rpMargin, sub.delta <= meta$cook_DB_info[["max_rp.delta"]]+rpMargin)
    
    if (nrow(cookCandLog) >= 2){
      # check time gap for each signal  
      timeGap <- difftime( tail(cookCandLog$start.timestamp,-1), head(cookCandLog$start.timestamp,-1), units='secs')
      
      ### cook Log processing
      # option 1] skip obstacle, i.e., foreseeing (to refine)
      # GapLen <- length(timeGap)
      # Gap_cumSum <- cumsum(as.numeric(timeGap))
      # chosenCandIdx <- sapply(seq(GapLen), function(x) {
      #   if(x == 1) tmp_Gap <- Gap_cumSum[seq(x, GapLen)] else tmp_Gap <- (Gap_cumSum[seq(x, GapLen)] - Gap_cumSum[x-1])
      #   checkCumTime <- ((cook_DBTime+timeTolerance) >= tmp_Gap) & ((cook_DBTime-timeTolerance) <= tmp_Gap)
      #   if(any(checkCumTime)) return(c(x, which(checkCumTime)+x )) else return(NULL)
      # }, simplify = FALSE)
      # chosenCandIdx <- unique(unlist(chosenCandIdx))
      
      # option 2] look only adjacent signal
      chosenCandIdx <- which( ( (cook_DBTime-timeTolerance) <= timeGap) & ((cook_DBTime+timeTolerance) >= timeGap) )
      chosenCandIdx <- unique( c(chosenCandIdx, chosenCandIdx+1) )
      
      if(length(chosenCandIdx) >= 2){
        chosenCand <- cookCandLog[sort(chosenCandIdx),]
        
        # split interval
        splitInterval <- cumsum( c(TRUE, difftime( tail(chosenCand$start.timestamp,-1),
                                                   head(chosenCand$start.timestamp,-1), units='secs') >= major_lumpGap))
        splitLumps <- bind_cols(chosenCand, data.frame(lumpIdx = splitInterval)) 
        lumpSummary <- ddply(splitLumps, .(lumpIdx), summarize, sum = length(start.timestamp), lump.start = min(start.timestamp), lump.end = max(end.timestamp), lump.duration = as.numeric(lump.end) - as.numeric(lump.start),
                             min.t = min(diff(as.numeric(start.timestamp))), med.t = median(diff(as.numeric(start.timestamp))), max.t = max(diff(as.numeric(start.timestamp))), min.med.rate2 = quantile(diff(as.numeric(start.timestamp)), .05)/med.t, med.h1 = median(h1),
                             min.d = min(delta), med.d = median(delta), max.d = max(delta), sd.d = sd(delta), med.sub.d = median(sub.delta), lost.sig.num = sum( pmax( round( diff(as.numeric(start.timestamp) ) / med.t ) -1 ,0)),
                             lost.sig.rate =  lost.sig.num*100/(sum+lost.sig.num) )
        lumpSummary <- lumpSummary %>% filter(sum >= eff_size) 
        if(nrow(lumpSummary) != 0){
          lumpLog <- splitLumps %>% filter(lumpIdx %in% lumpSummary$lumpIdx)
          # for pair annihilation
          if (existWarm) if(meta$warm_info[["type"]] == 1) cookLog_forAnnihilation <- lumpLog
          
          
          # energy estimation
          estiOutput <- ldply(seq(nrow(lumpSummary)), function(rowIdx){
            tmp_lumpIdx <- lumpSummary$lumpIdx[rowIdx]
            tmp_LumpPattern <- lumpLog %>% filter(lumpIdx == tmp_lumpIdx)
            return(energyEstimation_caseCook(data, metaHz, tmp_LumpPattern, meta$parameter[["proper_cookTime"]], cook_minTimeON))
          })
          # put sig fluc. num.: cbind(estiOutput, sigNum = rep(lumpSummary$sum, each = 2))
          
          # mark cook signal
          summary_cookTime <- lumpSummary %>% mutate(cookTime = lump.start + dseconds(lump.duration/2)) %>% select(cookTime)
          valid_cookLump <- estiOutput %>% mutate(validity = (SecON >= cook_minTimeON)) %>% bind_cols(data.frame(lumpIdx = rep(seq(nrow(summary_cookTime)), each = 2))) %>%
            group_by(lumpIdx) %>% summarise(validity = any(validity), SecON = median(SecON), whON = median(whON), apON = median(apON)) %>% bind_cols(summary_cookTime)
          
          ### cook signals w/ proper energy consumption
          if(any(valid_cookLump$validity == TRUE)){
            proper_cookInfo <- valid_cookLump %>% filter(validity == TRUE)
            cooking.start.timestamp <- floor_date(proper_cookInfo$cookTime, unit = "second")
            cooking.start.idx <- which(answer.log$timestamp %in% cooking.start.timestamp)
            
            # cook number count
            cookCnt <- length(cooking.start.timestamp)
            
            cooking.on.idx <- unique(unlist(lapply( cooking.start.idx, function(x) x + (-round(cook_usageStartTime):(round(cook_usageEndTime)-1)) )))
            ### avoid overflow of the idx
            cooking.on.idx <- cooking.on.idx[(cooking.on.idx >= 1) & (cooking.on.idx <= nrow(answer.log)) ]
            answer.log$p[cooking.on.idx] <- answer.log$p[cooking.on.idx] - (cook_usageWhON*3600/(cook_usageStartTime + cook_usageEndTime ))
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
      }
    }
  }
  
  # warm signal search
  if(existWarm){
    
    warm_DBTime <- meta$warm_info[["DB_TimeSec"]]
    warm_type <- meta$warm_info[["type"]]
    particleSize <- meta$parameter[["annihilation_particleLumpSize"]]
    # set reduce
    warmCandLog <- e_pattern %>% filter( h2 >= meta$warm_info[["min_ap.h2"]], h2 <= meta$warm_info[["max_ap.h2"]], 
                                         sub.delta >= meta$warm_info[["min_rp.delta"]], sub.delta <= meta$warm_info[["max_rp.delta"]])
    
    if(nrow(warmCandLog) >= 2){
      
      refinedCandLog <- data.frame()
      # step 1] consider time constraint 
      if(warm_type == 0){
        ### check time gap for each signal
        timeGap <- difftime( tail(warmCandLog$start.timestamp,-1), head(warmCandLog$start.timestamp,-1), units='secs')
        chosenCandIdx <- which( ( (warm_DBTime-timeTolerance) <= timeGap) & ((warm_DBTime+timeTolerance) >= timeGap) )
        chosenCandIdx <- unique( c(chosenCandIdx, chosenCandIdx+1) )
        if( length(chosenCandIdx) >= 2){
          refinedCandLog <- warmCandLog[sort(chosenCandIdx),]
        }
      } else if(warm_type == 1){
        refinedCandLog <- warmCandLog
      }
      
      # step 2] identify lumps
      if(nrow(refinedCandLog) != 0){
        
        splitInterval <- cumsum( c(TRUE, difftime( tail(refinedCandLog$start.timestamp,-1),
                                                   head(refinedCandLog$start.timestamp,-1), units='secs') >= major_lumpGap))
        splitLumps <- bind_cols(refinedCandLog, data.frame(lumpIdx = splitInterval))
        lumpSummary <- ddply(splitLumps, .(lumpIdx), summarize, sum = length(start.timestamp))
        if(warm_type == 0){
          eff_lumpSize <- eff_size
        } else if(warm_type == 1){
          eff_lumpSize <- particleSize
        }
        lumpSummary <- lumpSummary %>% filter(sum >= eff_lumpSize)
        refinedCandLog <- splitLumps %>% filter( lumpIdx %in% lumpSummary$lumpIdx)
      }
      
      # step 3] signal extension
      if( (nrow(refinedCandLog) != 0) && (length(sigGroupInfo) != 0) ){
        
        h2Min <- min(refinedCandLog$h2)
        h2Max <- max(refinedCandLog$h2)
        sub.deltaMin <- min(refinedCandLog$sub.delta)
        sub.deltaMax <- max(refinedCandLog$sub.delta)
        
        extCand <- lapply(seq(length(sigGroupInfo) -1), function(x) return(sigGroupInfo[[x]] %>% filter(h2 >= h2Min, h2 <= h2Max, sub.delta >= sub.deltaMin, sub.delta <= sub.deltaMax)))
        extCandNum <- sapply(extCand, function(x) return(nrow(x)))
        
        if(max(extCandNum) != 0){
          chosenExtCand <- extCand[[which.max(extCandNum)]]
          lump_timeIdx <- refinedCandLog %>% group_by(lumpIdx) %>% summarise(min_startIdx = min(start.idx), max_endIdx = max(end.idx))
          newOutput <- apply(lump_timeIdx, 1, function(x) return(chosenExtCand %>% filter(start.idx >= as.numeric(x['min_startIdx']), end.idx <= as.numeric(x['max_endIdx'])) )) %>% bind_rows() %>% 
            ungroup() %>% select(-xDivision, -yDivision)
          
          finWarmLog <- refinedCandLog %>% select(-lumpIdx) %>% bind_rows(newOutput) 
          refinedCandLog <- finWarmLog[!duplicated(finWarmLog$start.idx),] %>% arrange(start.idx)
        }
      }
      
      # step 4] remove cook signal for case of pair annihilation (if possible)
      if( (nrow(refinedCandLog) != 0) && !is.null(cookLog_forAnnihilation)){
        compareCookLog <- cookLog_forAnnihilation %>% ungroup() %>% select(start.idx, end.idx, start.timestamp, end.timestamp, h2, sub.delta)
        compareWarmLog <- refinedCandLog %>% select(start.idx, end.idx, start.timestamp, end.timestamp, h2, sub.delta)
        refinedCandLog <- dplyr::setdiff(compareWarmLog, compareCookLog)
      }
      
      # step 5] insert signals (for results)
      if(nrow(refinedCandLog) != 0){
        
        new_timestamp <- floor_date(refinedCandLog$start.timestamp, unit = "second")
        warming.start.idx <- which(answer.log$timestamp %in% new_timestamp)
        warming.on.idx <- unique(unlist(lapply( warming.start.idx, function(x) x + 0:(round(meta$warm_info[["usage_secON"]])-1) )))
        # avoid overflow of the idx
        warming.on.idx <- warming.on.idx[warming.on.idx <= nrow(answer.log)]
        answer.log$p[warming.on.idx] <- answer.log$p[warming.on.idx] + meta$warm_info[["usage_watt"]] * usage_powerAdjust
        
      }
    }
  }
  
  return(list(usage = answer.log, confidence = meta$reliability, count = cookCnt, samplingRate = 1.))
}


riceCookerCore_jp <- function (data, prior_meta = NULL, Hz = 10){
  
  # options(scipen = 15)
  
  # major parameters
  eff_size = 7
  periodicity = 0.2
  # min_cookFluc = 9
  
  # signal search: major case
  search_iter <- 650 # increased by 50 due to the JPN site '10012085'JUNE DATA
  quantileProb <- 0.25
  major_min_watt <- 250
  major_max_watt <- 1350
  major_medSec_min <- 10
  major_medSec_max <- 700 # 600 -> 700 for SITE NO. 10012078, JULY DATA 
  major_lumpGap <- 3600
  
  # DB 
  # flag - 0: cooking || 1: warming || 2: both (NOT USING YET!)
  DB_major <- data.frame(time = c(15, 16, 20, 32, 45, 50, 60, 75, 128), 
                         flag = c( 0,  2,  2,  1,  1,  1,  1,  1,   1))
  DB_toleranceSec <- 0.4
  
  # signal orchestration
  thres_timeDuration <- 3600 * 2
  proper_cookTime <- 40 * 60
  min_cookOn <- 6 * 60
  ap_ignoreThres <- 5 
  
  # signal annihilation
  annihilation_cookNumThres <- 2
  annihilation_spanSec <- 3*3600
  annihilation_thresProb <- 0.9
  annihilation_maxApGap <- 300 # set loosely 
  annihilation_apMinThres <- 50
  annihilation_apMaxThres <- major_max_watt # set to 1600
  annihilation_particleLumpSize <- 8 # might be larger than eff_size
  annihilation_coverageProb <- 0.75
  
  #########################################################################
  ### DetectPattern_1Hz_new:: cowork with SJLEE
  e_pattern <- DetectPattern_1Hz_new(data, position = "end", main_type = "active", sub_type = "reactive", useConsistencyFlag = F)
  #########################################################################
  
  majorSig_info <- sigSearch_majorCase( patterns = e_pattern, min_watt = major_min_watt, max_watt = major_max_watt, group_periodicity = periodicity, group_quantProb = quantileProb, 
                                        group_searchIter = search_iter, group_medSec_min = major_medSec_min, group_medSec_max = major_medSec_max, effective_size = eff_size, 
                                        lump_timeGapThres = major_lumpGap, timeTolerance = DB_toleranceSec, DB_majorSig = DB_major)
  
  validMeta <- sigOrchestration_riceCooker(data, Hz, answer.log, major_SigInfo = majorSig_info, major_DBInfo = DB_major, thresTime = thres_timeDuration, 
                                           range_energyEsti = proper_cookTime, min_usageTime = min_cookOn, min_apThres = ap_ignoreThres, annihilation_cookNumThres, annihilation_pattern = e_pattern, 
                                           annihilation_spanSec, annihilation_apMinThres, annihilation_apMaxThres, annihilation_searchIter = search_iter, annihilation_effSize = eff_size, annihilation_medSec_min = major_medSec_min,
                                           annihilation_medSec_max = major_medSec_max, annihilation_thresProb, annihilation_maxApGap, annihilation_coverageProb, lump_timeGapThres = major_lumpGap, annihilation_particleLumpSize)
  
  if (length(validMeta) == 0) {
    print("no valid signal: reusing prior_meta")
    return(prior_meta)
  }
  
  # if (validMeta$general[["Cook"]]){
  #   if(validMeta$cook[["minFlucs_properLump"]] < min_cookFluc){
  #     print("too small cook fluctuation: reusing prior_meta")
  #     return(prior_meta)
  #   }
  # }
  
  
  # parameter setup
  validMeta$parameter <- c("Hz" = Hz, "eff_size" = eff_size, "periodicity" = periodicity, "search_iter" = search_iter, "quantileProb" = quantileProb, "major_min_watt" = major_min_watt,
                           "major_max_watt" = major_max_watt, "major_medSec_min" = major_medSec_min, "major_medSec_max" = major_medSec_max, "DB_toleranceSec" = DB_toleranceSec, 
                           "major_lumpGap" = major_lumpGap, "proper_cookTime" = proper_cookTime, "min_cookOn" = min_cookOn, "annihilation_particleLumpSize" = annihilation_particleLumpSize)
  
  validMeta$DBTime <- DB_major$time
  
  # compute reliability
  metaReliability <- 0.2
  if(validMeta$general[["Cook"]] == TRUE) metaReliability <- metaReliability + 0.4
  if(validMeta$general[["Warm"]] == TRUE) metaReliability <- metaReliability + 0.2
  validMeta$reliability <- metaReliability
  
  return(validMeta)
}

