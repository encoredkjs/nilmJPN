convert_15Hz_to_1Hz <- function(dt){
  ddply(dt, .(appliance, timestamp = lubridate::floor_date(timestamp,'second')), 
        summarize, active_power = median(active_power), reactive_power = median(reactive_power))
}


kimchi.usage <- function(testData.1Hz, solutionData.1Hz, debug.mode = F){
  
  #testData <- convert_15Hz_to_1Hz( testData )
  
  meta.json <- "[\n {\n \"meta-version\":      1,\n\"shape_type\": \"cyclic_box\",\n\"rising_edge\": {\n \"Delta.reactive_power_height\": 241.84,\n\"Delta.reactive_power_sigma\": 2.5658,\n\"Delta.reactive_power_min\": 237.48,\n\"Delta.reactive_power_max\": 245.65,\n\"Delta.reactive_power_med\": 242.07,\n\"Delta.active_power_height\": 250.49,\n\"Delta.active_power_sigma\": 4.9715,\n\"Delta.active_power_min\": 243.59,\n\"Delta.active_power_max\": 259.18,\n\"Delta.active_power_med\": 249.72 \n},\n\"falling_edge\": {\n \"Delta.reactive_power_height\": -236.87,\n\"Delta.reactive_power_sigma\": 2.6681,\n\"Delta.reactive_power_min\": -241.39,\n\"Delta.reactive_power_max\": -232.36,\n\"Delta.reactive_power_med\": -236.78,\n\"Delta.active_power_height\": -196.39,\n\"Delta.active_power_sigma\": 1.9725,\n\"Delta.active_power_min\": -199.28,\n\"Delta.active_power_max\": -193.5,\n\"Delta.active_power_med\": -196.75 \n},\n\"cycle\": {\n \"working_time\":    122,\n\"working_time_sigma\": 5.1062,\n\"working_time_min\":    106,\n\"working_time_max\":    127,\n\"duty_cycle\": 0.10164,\n\"tot_time\":   1145 \n},\n\"box.no\": {\n \"n.pt\":     71 \n},\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-14\",\n\"end\": \"2015-11-14 23:59:55\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-20 14:45:39\" \n},\n\"level\":      3 \n} \n]"
  result <- multipleMeta2cyclicBoxExtend( fromJSON(meta.json), testData.1Hz, keepSingleEdge = T, 
                                minMag.reactive = 10, postProcessing = T )[[1]]
  result$p <- result$p + 4.559
#  base_power  <- 0.142 * nrow(data)
  
  if( !missing(solutionData.1Hz) ){
    
    #solutionData <- convert_15Hz_to_1Hz( solutionData )

    print( c("Estimated val = ", sum(result$p)))
    print( c("    Exact val = ", sum(solutionData.1Hz$active_power)))
    print( c("  Usage ratio = ", sum(result$p) / sum(solutionData.1Hz$active_power)))
    
    if( debug.mode ){
      plot( solutionData.1Hz$timestamp, solutionData.1Hz$active_power, type='l' )
      lines( result$timestamp, result$p, col='red' )
    }
  }
  
  return(result)
}


refrigerator.usage.15Hz <- function(testData.15Hz, solutionData.15Hz, debug.mode = F){
  
  #testData <- convert_15Hz_to_1Hz( testData )
  #solutionData <- convert_15Hz_to_1Hz( solutionData )

   meta.json <- "[\n {\n \"meta-version\":      1,\n\"shape_type\": \"cyclic_box\",\n\"rising_edge\": {\n \"Delta.reactive_power_height\": 71.095,\n\"Delta.reactive_power_sigma\": 5.4489,\n\"Delta.reactive_power_min\": 67.713,\n\"Delta.reactive_power_max\": 79.841,\n\"Delta.reactive_power_med\": 69.357,\n\"Delta.active_power_height\": 164.44,\n\"Delta.active_power_sigma\": 29.823,\n\"Delta.active_power_min\": 148.77,\n\"Delta.active_power_max\": 214.83,\n\"Delta.active_power_med\": 154.99 \n},\n\"falling_edge\": {\n \"Delta.reactive_power_height\": -63.559,\n\"Delta.reactive_power_sigma\": 0.7018,\n\"Delta.reactive_power_min\": -64.743,\n\"Delta.reactive_power_max\": -62.644,\n\"Delta.reactive_power_med\": -63.506,\n\"Delta.active_power_height\": -117.86,\n\"Delta.active_power_sigma\": 0.75555,\n\"Delta.active_power_min\": -119.17,\n\"Delta.active_power_max\": -116.99,\n\"Delta.active_power_med\": -117.72 \n},\n\"cycle\": {\n \"working_time\":  39834,\n\"working_time_sigma\": 2699.3,\n\"working_time_min\":  37002,\n\"working_time_max\":  45116,\n\"duty_cycle\": 0.44949,\n\"tot_time\":  88611 \n},\n\"box.no\": {\n \"n.pt\":     14 \n},\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-14 00:00:00\",\n\"end\": \"2015-11-14 23:59:59\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 10:39:02\" \n},\n\"level\":      3 \n} \n]"
  
   result <- multipleMeta2cyclicBoxExtend( fromJSON(meta.json), testData.15Hz, keepSingleEdge = T, 
                                           minMag.reactive = 10, postProcessing = T )[[1]]
   result$p <- result$p + 0.114 # 0.796
   
   if( !missing(solutionData.15Hz) ){
     print( c("Estimated val = ", sum(result$p)))
     print( c("    Exact val = ", sum(solutionData.15Hz$active_power)))
     print( c("  Usage ratio = ", sum(result$p) / sum(solutionData.15Hz$active_power)))
     
     if( debug.mode ){
       plot( solutionData.15Hz$timestamp, solutionData.15Hz$active_power, type='l' )
       lines( result$timestamp, result$p, col='red' )
     }
   }
   return(result)
}

purifier.usage.15Hz <- function(testData.1Hz, solutionData.1Hz, testData.15Hz, solutionData.15Hz, debug.mode = F){
  
  # heater
  #meta.json1 <- "[\n {\n \"meta-version\":      1,\n\"shape_type\": \"cyclic_box\",\n\"rising_edge\": {\n \"Delta.active_power_height\": 525.39,\n\"Delta.active_power_sigma\": 6.1226,\n\"Delta.active_power_min\": 515.37,\n\"Delta.active_power_max\":  534.2,\n\"Delta.active_power_med\": 524.91,\n\"Delta.reactive_power_height\": 0.52871,\n\"Delta.reactive_power_sigma\": 0.34328,\n\"Delta.reactive_power_min\":  0.491,\n\"Delta.reactive_power_max\": 0.6242,\n\"Delta.reactive_power_med\":  0.564 \n},\n\"falling_edge\": {\n \"Delta.active_power_height\": -524.86,\n\"Delta.active_power_sigma\": 5.9155,\n\"Delta.active_power_min\": -533.78,\n\"Delta.active_power_max\": -515.44,\n\"Delta.active_power_med\": -524.58,\n\"Delta.reactive_power_height\": -0.55272,\n\"Delta.reactive_power_sigma\": 0.37746,\n\"Delta.reactive_power_min\": -0.6445,\n\"Delta.reactive_power_max\": -0.5284,\n\"Delta.reactive_power_med\": -0.599 \n},\n\"cycle\": {\n \"working_time\":    415,\n\"working_time_sigma\":  16.15,\n\"working_time_min\":    375,\n\"working_time_max\":    446,\n\"duty_cycle\": 0.084543,\n\"tot_time\": 4828.5 \n},\n\"box.no\": {\n \"n.pt\":     85 \n},\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11\",\n\"end\": \"2015-11-17\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 13:20:01\" \n},\n\"level\":      4 \n} \n]"
  meta.json1 <- "[\n {\n \"meta-version\":      1,\n\"shape_type\": \"cyclic_box\",\n\"rising_edge\": {\n \"Delta.active_power_height\": 525.39,\n\"Delta.active_power_sigma\": 6.1226,\n\"Delta.active_power_min\": 515.37,\n\"Delta.active_power_max\":  534.2,\n\"Delta.active_power_med\": 524.91,\n\"Delta.reactive_power_height\": 0.52871,\n\"Delta.reactive_power_sigma\": 0.34328,\n\"Delta.reactive_power_min\":  0.491,\n\"Delta.reactive_power_max\": 0.6242,\n\"Delta.reactive_power_med\":  0.564 \n},\n\"falling_edge\": {\n \"Delta.active_power_height\": -524.86,\n\"Delta.active_power_sigma\": 5.9155,\n\"Delta.active_power_min\": -533.78,\n\"Delta.active_power_max\": -515.44,\n\"Delta.active_power_med\": -524.58,\n\"Delta.reactive_power_height\": -0.55272,\n\"Delta.reactive_power_sigma\": 0.37746,\n\"Delta.reactive_power_min\": -0.6445,\n\"Delta.reactive_power_max\": -0.5284,\n\"Delta.reactive_power_med\": -0.599 \n},\n\"cycle\": {\n \"working_time\":   6225,\n\"working_time_sigma\":  16.15,\n\"working_time_min\":   5625,\n\"working_time_max\":   6690,\n\"duty_cycle\": 0.084543,\n\"tot_time\": 4828.5 \n},\n\"box.no\": {\n \"n.pt\":     85 \n},\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11\",\n\"end\": \"2015-11-17\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 13:20:01\" \n},\n\"level\":      4 \n} \n]"
  
  result1 <- multipleMeta2cyclicBoxExtend( fromJSON(meta.json1), testData.15Hz, useSoftMatching = F, extension.p = .1, keepSingleEdge = T,   
                                           ignoreRP = T, postProcessing = T )[[1]]
  
  #testData <- convert_15Hz_to_1Hz( testData )
  
  # motor
  #meta.json2 <- "[\n {\n \"meta-version\":      1,\n\"shape_type\": \"cyclic_box\",\n\"rising_edge\": {\n \"Delta.reactive_power_height\": 167.62,\n\"Delta.reactive_power_sigma\": 2.4371,\n\"Delta.reactive_power_min\": 163.78,\n\"Delta.reactive_power_max\": 171.49,\n\"Delta.reactive_power_med\": 167.44,\n\"Delta.active_power_height\": 138.34,\n\"Delta.active_power_sigma\": 8.6997,\n\"Delta.active_power_min\": 125.56,\n\"Delta.active_power_max\": 151.14,\n\"Delta.active_power_med\": 139.03 \n},\n\"falling_edge\": {\n \"Delta.reactive_power_height\": -162.44,\n\"Delta.reactive_power_sigma\": 2.5571,\n\"Delta.reactive_power_min\": -166.33,\n\"Delta.reactive_power_max\": -158.82,\n\"Delta.reactive_power_med\": -162.27,\n\"Delta.active_power_height\": -139.05,\n\"Delta.active_power_sigma\": 1.7021,\n\"Delta.active_power_min\": -140.97,\n\"Delta.active_power_max\": -136.05,\n\"Delta.active_power_med\": -139.31 \n},\n\"cycle\": {\n \"working_time\":   22410,\n\"working_time_sigma\": 62.354,\n\"working_time_min\":   20730,\n\"working_time_max\":   24375,\n\"duty_cycle\": 0.075759,\n\"tot_time\":  11354 \n},\n\"box.no\": {\n \"n.pt\":     44 \n},\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11\",\n\"end\": \"2015-11-17\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 13:19:59\" \n},\n\"level\":      3 \n} \n]"
  meta.json2 <- "[\n {\n \"meta-version\":      1,\n\"shape_type\": \"cyclic_box\",\n\"rising_edge\": {\n \"Delta.reactive_power_height\": 167.62,\n\"Delta.reactive_power_sigma\": 2.4371,\n\"Delta.reactive_power_min\": 163.78,\n\"Delta.reactive_power_max\": 171.49,\n\"Delta.reactive_power_med\": 167.44,\n\"Delta.active_power_height\": 138.34,\n\"Delta.active_power_sigma\": 8.6997,\n\"Delta.active_power_min\": 125.56,\n\"Delta.active_power_max\": 151.14,\n\"Delta.active_power_med\": 139.03 \n},\n\"falling_edge\": {\n \"Delta.reactive_power_height\": -162.44,\n\"Delta.reactive_power_sigma\": 2.5571,\n\"Delta.reactive_power_min\": -166.33,\n\"Delta.reactive_power_max\": -158.82,\n\"Delta.reactive_power_med\": -162.27,\n\"Delta.active_power_height\": -139.05,\n\"Delta.active_power_sigma\": 1.7021,\n\"Delta.active_power_min\": -140.97,\n\"Delta.active_power_max\": -136.05,\n\"Delta.active_power_med\": -139.31 \n},\n\"cycle\": {\n \"working_time\":   1494,\n\"working_time_sigma\": 62.354,\n\"working_time_min\":   1382,\n\"working_time_max\":   1625,\n\"duty_cycle\": 0.075759,\n\"tot_time\":  11354 \n},\n\"box.no\": {\n \"n.pt\":     44 \n},\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11\",\n\"end\": \"2015-11-17\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 13:19:59\" \n},\n\"level\":      3 \n} \n]"
  
  result2 <- multipleMeta2cyclicBoxExtend( fromJSON(meta.json2), testData.1Hz, keepSingleEdge = T, 
                                           postProcessing = T )[[1]]

#  result1 <- ddply(result1, .(timestamp = lubridate::floor_date(timestamp,'second')), summarize, p = median(p))
#  result2 <- ddply(result2, .(timestamp = lubridate::floor_date(timestamp,'second')), summarize, p = median(p))

#  result <- merge( result1, result2, by='timestamp' )
#  result$p <- result$p.x + result$p.y
  
    if( !missing(solutionData.1Hz) & !missing(solutionData.15Hz)  ){
  
    #solutionData <- convert_15Hz_to_1Hz( solutionData )
  
    print( c("Estimated val = ", sum(result1$p) + sum(result2$p)))
    print( c("    Exact val = ", sum(solutionData.15Hz$active_power)))
    print( c("  Usage ratio = ", (sum(result1$p) + sum(result2$p)) / sum(solutionData.15Hz$active_power)))
    
    if( debug.mode ){
      plot( solutionData.15Hz$timestamp, solutionData.15Hz$active_power, type='l' )
      lines( result1$timestamp, result1$p, col='red' )
      lines( result2$timestamp, result2$p, col='red' )
    }
    }
#  return(result)
  return(list(result1,result2))
}


iron.usage.15Hz <- function(testData.15Hz, solutionData.15Hz, debug.mode = F){
  
  #testData <- convert_15Hz_to_1Hz( testData )
  
  #meta.json <- "{\n \"meta-version\":      1,\n\"shape_type\": \"pattern_scan_heavy\",\n\"summit_flag\":     -1,\n\"rising_edge\": {\n \"ap_min\":   1769,\n\"ap_median\": 1814.4,\n\"ap_max\": 1863.4,\n\"rp_min\": -5.1635,\n\"rp_median\": -4.849,\n\"rp_max\": -4.642 \n},\n\"falling_edge\": {\n \"ap_min\":  -1855,\n\"ap_median\": -1796.2,\n\"ap_max\": -1756.7,\n\"rp_min\":  4.717,\n\"rp_median\": 4.8625,\n\"rp_max\":  5.162 \n},\n\"duration\": {\n \"duration_min\":     12,\n\"duration_median\":     14,\n\"duration_max\":     65 \n},\n\"n.pt\":     36,\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11 13:39:46\",\n\"end\": \"2015-11-15 20:35:57\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 17:29:15\" \n} \n}"
  meta.json <- "{\n \"meta-version\":      1,\n\"shape_type\": \"pattern_scan_heavy\",\n\"summit_flag\":     2,\n\"rising_edge\": {\n \"ap_min\":   1769,\n\"ap_median\": 1814.4,\n\"ap_max\": 1863.4,\n\"rp_min\": -5.1635,\n\"rp_median\": -4.849,\n\"rp_max\": -4.642 \n},\n\"falling_edge\": {\n \"ap_min\":  -1855,\n\"ap_median\": -1796.2,\n\"ap_max\": -1756.7,\n\"rp_min\":  4.717,\n\"rp_median\": 4.8625,\n\"rp_max\":  5.162 \n},\n\"duration\": {\n \"duration_min\":     12,\n\"duration_median\":     14,\n\"duration_max\":     65 \n},\n\"n.pt\":     36,\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11 13:39:46\",\n\"end\": \"2015-11-15 20:35:57\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 17:29:15\" \n} \n}"
  meta.json <- "{\n \"meta-version\":      1,\n\"shape_type\": \"pattern_scan_heavy\",\n\"summit_flag\":      2,\n\"rising_edge\": {\n \"ap_min\": 1789.4,\n\"ap_median\": 1819.8,\n\"ap_max\": 1871.2,\n\"rp_min\": -1.2474,\n\"rp_median\": 0.71477,\n\"rp_max\": 1.1564 \n},\n\"falling_edge\": {\n \"ap_min\": -1827.2,\n\"ap_median\": -1791.6,\n\"ap_max\": -1776.8,\n\"rp_min\": -1.0332,\n\"rp_median\": -0.70513,\n\"rp_max\": 1.2638 \n},\n\"duration\": {\n \"duration_min\":     14,\n\"duration_median\":   31.5,\n\"duration_max\":     70 \n},\n\"n.pt\":      4,\n\"generation_info\": {\n \"data_used\": {\n \"start\": \"2015-11-11 13:39:46\",\n\"end\": \"2015-11-15 20:35:57\",\n\"sampling\":      1 \n},\n\"computed\": \"2015-11-23 19:04:40\" \n} \n}"
  result <- meta2PatternScanHeavy_1Hz( testData.15Hz, fromJSON(meta.json), debug.mode = T, ignoreRP = T, extension.p = .2 )
  boxSummary <- series.to.box.lists( result$timestamp, result$p, 0 )[[1]]
  
  if( any(boxSummary$duration < 10) ){
    boxSummary <- subset( boxSummary, duration < 10 )
    for( i in 1:nrow(boxSummary) ){
      result$p[boxSummary$str[i]:boxSummary$end[i]] <- 0
    }
  }
  
  if( !missing(solutionData.15Hz) ){
    #solutionData <- convert_15Hz_to_1Hz( solutionData )
    
    print( c("Estimated val = ", sum(result$p)))
    print( c("    Exact val = ", sum(solutionData.15Hz$active_power)))
    print( c("  Usage ratio = ", sum(result$p) / sum(solutionData.15Hz$active_power)))
    
    if( debug.mode ){
      plot( solutionData.15Hz$timestamp, solutionData.15Hz$active_power, type='l' )
      lines( result$timestamp, result$p, col='red' )
    }
  }
  return( result )
}