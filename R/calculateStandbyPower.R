calculateStandbyPower <- function(data)
  min(envelopeDetector(data)$UpperEnvelope) * length(data$active_power)