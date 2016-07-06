from.para.to.JSON <- function( data, r.edge, f.edge, on.off, box.no, locality = -1 ){
  
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  
  json.result <- list( 'meta-version' = 1, 
                       'shape_type'   = 'cyclic_box', 
                       'rising_edge'  = r.edge, 
                       'falling_edge' = f.edge, 
                       'cycle'        = on.off, 
                       'box.no'       = box.no, 
                       'locality'     = locality, 
                       'generation_info' = list('data_used' = list( 'start' = str.t, 'end' = end.t, 'sampling' = 1),
                                                'computed' = as.character(Sys.time())))
  return(json.result)
}

from.para.to.JSON2 <- function( data, r.edge, f.edge, box.no ){
  
  if( missing(f.edge) ){
    f.edge <- r.edge
    f.edge['ap_height'] <- -r.edge['ap_height']
    f.edge['rp_height'] <- -r.edge['rp_height']
  }
  
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))

  json.result <- list( 'meta-version' = 1, 
                       'shape_type'   = 'high_power', 
                       'rising_edge'  = r.edge, 
                       'falling_edge' = f.edge, 
                       'box.no'       = box.no, 
                       'generation_info' = list('data_used' = list( 'start' = str.t, 'end' = end.t, 'sampling' = 1),
                                                'computed' = as.character(Sys.time())))
  return(json.result)
}

from.para.to.JSON3 <- function( data, StandbyPowerValue ){
  
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  
  json.result <- list( 'meta-version' = 1, 
                       'shape_type'   = 'StandbyPower', 
                       'value'        = StandbyPowerValue, 
                       'generation_info' = list('data_used' = list( 'start' = str.t, 'end' = end.t, 'sampling' = 1),
                                                'computed' = as.character(Sys.time())))
  return(json.result)
}

from.para.to.JSON4 <- function( data, StandbyPowerValue ){
  
  str.t <- as.character(min(data$timestamp))
  end.t <- as.character(max(data$timestamp))
  
  json.result <- list( 'meta-version' = 1, 
                       'shape_type'   = 'StandbyPower', 
                       'value'        = StandbyPowerValue, 
                       'generation_info' = list('data_used' = list( 'start' = str.t, 'end' = end.t, 'sampling' = 1),
                                                'computed' = as.character(Sys.time())))
  
  tmp_start.json <- c( summit_flag = sign_flag.b
                       , rp_min = min(s.data.log$delta)
                       , rp_max = max(s.data.log$delta)
                       , rp_med = median(s.data.log$delta)
                       , ap_min = min(s.data.log$sub.delta)
                       , ap_max = max(s.data.log$sub.delta)
                       , ap_med = median(s.data.log$sub.delta)
                       , sample.num = s.list$sum[chosen.e.idx[m_idx]]
                       , min.t = s.list$min.t[chosen.e.idx[m_idx]]
                       , med.t = s.list$med.t[chosen.e.idx[m_idx]]
                       , min.med.rate = s.list$min.med.rate[chosen.e.idx[m_idx]]
                       , lost.sample.num = s.list$lost.sig.num[chosen.e.idx[m_idx]]
                       , lost.sample.rate = s.list$lost.sig.rate[chosen.e.idx[m_idx]])
  
  # end data info
  tmp_end.json <- c( rp_min = min(e.data.log$h2), 
                     rp_max = max(e.data.log$h2),
                     ap_min = min(e.data.log$sub.delta), 
                     ap_max = max(e.data.log$sub.delta),
                     ap_med = median(e.data.log$sub.delta), 
                     slotNum.zero = e.data.summary$slot_zeros, 
                     slotNum.one  = e.data.summary$slot_ones, 
                     ZerotoOneratio = e.data.summary$metric.info, 
                     EffTimeOn.med = e.data.summary$eff.TimeOn.med, 
                     EffTimeOn.min = e.data.summary$eff.TimeOn.min, 
                     EffTimeOn.max = e.data.summary$eff.TimeOn.max,
                     EffTimeOn.sd = e.data.summary$eff.TimeOn.sd, 
                     EffPwrDrop.med = e.data.summary$eff.PowerDrop.med, 
                     EffPwrDrop.min = e.data.summary$eff.PowerDrop.min, 
                     EffPwrDrop.max = e.data.summary$eff.PowerDrop.max, 
                     EffPwrDrop.sd = e.data.summary$eff.PowerDrop.sd, 
                     EffRP_Drop.med = e.data.summary$eff.rpDrop.med, 
                     EffRP_Drop.min = e.data.summary$eff.rpDrop.min, 
                     EffRP_Drop.max = e.data.summary$eff.rpDrop.max, 
                     EffRP_Drop.sd = e.data.summary$eff.rpDrop.sd)
  
  return(json.result)
}
