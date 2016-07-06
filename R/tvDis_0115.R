#' @title tvDis
#' @author CH
#' @description disaggregate TV from total data
#' @update 2015.12.03
#' @input c(timestamp, active_power, reactive_power)
#' @output c(timestamp, p, q) : tv timestamp, active_power, reactive_power
#' procedure : feature extraction -> classification -> reconstruction
#'

#source(file = 'tvFeatureExtractor.R')
#source(file = 'connectMethod.R')
#source(file = 'tvReconstructor.R')
#source(file="edgeDetector.R")
#source(file="cutBlocks.R")

library(stringr)
library(xgboost)
library(pracma)

# load trained classifier
# load("/home/nilm/CloudNilm/libs/tvClassifier.rda")

# Training Mode (Fit learning model)
generate.meta.tv <- function(X = NULL, ...){

  ###
  # dummy function
  # copied and modified from generate.meta.washer
  if(!is.data.frame(X)) {
    X <- as.data.frame(X)
  }
  #X <- tbl_df(X)

  # Get tuning parameters
  TVParas = JSONforTV(...)

  return(TVParas)
}


JSONforTV <- function(metaVer = 1, startTS = "2015-11-01 00:00:00", endTS ="2015-11-01 23:59:59",
                      computedTime = as.character(Sys.time())){
  ###
  # dummy function
  # copied and modified from washer

  resJSON = list()

  #####
  # META VERSION
  resJSON['meta-version'] = metaVer

  #####
  resJSON$shape_type <- "TV"

  resJSON$generation_info = list(
    data_used = list(
      start = startTS,
      end = endTS
    ),
    samplingRate = 15,
    computed = computedTime)

  #class(resJSON) <- c(class(resJSON), "WashingMachine")
  return(resJSON)

}

meta2tv <- function(data, meta){
  tvDis(data)
}


tvDis <- function(data){

  # simple preprocessing : time stamp
  timeStamp_original <- data[, 1]
  timeStamp <- as.POSIXlt(data[ ,1], "%Y-%m-%d %H:%M:%OS8", tz = "GMT")
  timeStamp <- as.numeric(timeStamp)
  timeStamp <- (timeStamp - timeStamp[1])/3600
  data[ ,1] <- timeStamp
  colnames(data) <- c('timestamp', 'active_power', 'reactive_power')



  # parameters
  blockSize <- 15*60*5 #5minutes
  operationPoint <- 0.8 # tv On threshold, predicion>operationPoint -> On
  tvApValue <- 80 # tvOn : active power
  tvRpValue <- 0 # tvOn : reactive power
  connectThreshold <- 6 # 6*5 = 30min

  tvFeatures <- tvFeatureExtractor(data, blockSize)

  prediction <- predict(xgboost.fit, as.matrix(tvFeatures))
  prediction_connect30 <- connectMethod(prediction, operationPoint, 6)

  tv_disAp <- tvReconstructor(prediction_connect30, data, blockSize, tvApValue)
  tv_disRp <- rep(tvRpValue, length(tv_disAp))
  tvDisaggregated <- data.frame( timestamp = timeStamp_original,
                                 p = tv_disAp,
                                 q = tv_disRp)
  return(tvDisaggregated)
}



######################################################################################################################
#' @title connectMethod
#' @author CH
#' @description upper layer classification algorithm. if the distance between two On-states is shorter than threshold,
#'              connect On-state in a row.
#'              e.g. distanceThreshold 4
#'              0.85 0.5 0.4 0.9 0.5 0.2 -> 1 1 1 1 0.5 0.2
#'
#' @update 2015.12.03
#' @input predictValue : classification result probability
#'        operationPoint : decision threshold, prob>0.8 -> On
#'        distanceThreshold : the distance between two On-states threshold
#' @output modified classification vector





connectMethod <- function(predictValues, operationPoint, distanceThreshold){

  modifiedPredictValues <- predictValues;

  # On-state index list
  onList <- which(predictValues > operationPoint)

  if (length(onList)>0){
    for (i in 1:(length(onList)-1))
    {
      currentOn <- onList[i]
      nextOn <- onList[i+1]

      if ((nextOn-currentOn)<=distanceThreshold){
        modifiedPredictValues[currentOn:nextOn] <- 1
      } else{
        modifiedPredictValues[currentOn] <- 0

        if (i==(length(onList)-1)){ # last two block
          modifiedPredictValues[nextOn] <- 0
        }
      }

    }
  }

  return(modifiedPredictValues)
}

######################################################################################################################

#' @title cutBlocks
#' @author CH
#' @description cut block from sequence
#' @update 2015.12.03
#' @input ap, rp or ts sequence
#' @output block Table(floor(length(ap)/blockSize)*blockSize)
#'


cutBlocks <- function(ap, blockSize)
{

  blockTable = matrix(0, floor(length(ap)/blockSize), blockSize);
  for (i in 1:floor(length(ap)/blockSize)){
    blockTable[i,] = ap[(blockSize*(i-1)+1):(blockSize*i)];
  }

  return(blockTable)
}

######################################################################################################################

#' @title edgeDetector
#' @author (author) prof. / (interpreter) CH
#' @description interpreted version of MATLAB edgeDetector
#' @update 2015.12.01
#' @input the result of getDataFromCSV(ts, rp, ap, idx table)
#' @output transient vectors
#'         type : edge type, e.g. 1X -> rising edge, 2X -> falling edge, 3X-> impulse, -1 -> unknown
#'         startIndex : transient start index
#'         endIndex : transient end index
#'         endStartAp : steadystateHeightAp in old version, (Ap after transient - Ap before transient), for impulse check.
#'         endStartRp : steadystateHeightRp in old version, (Rp after transient - Rp before transient), for impulse check.
#'         maxStartAp : maximum Ap after the start of tranisent
#'         maxStartRp : maximum Rp after the start of tranisent#'
#'         mixStartAp : minimum Ap after the start of tranisent
#'         mixStartRp : minimum Rp after the start of tranisent
#'         numTransAp : the number of points that change Ap compared to (i-1)th point. e.g.) 1 1 3 ->1, 2 0 5 -> 2
#'         numSegmentsAp : the number of segments. e.g. RF -> 2, R->1, RFRRFR ->6
#'         maxSlopeAp : maximum Ap slope
#'         maxSlopeRp : maximum Rp slope
#'         minSlopeAp : minimum Ap slope
#'         minSlopeRp : minimum Rp slope
#'         edgeSequence : the sequence of edges in transient. e.g. RF, R, , RRF
#' @update 2015.12.01 : flatSectionList가 비어?????? 경우 처리.



edgeDetectorTV <- function(data) {
  ##
  # EDGEDETECTOR Summary of this function goes here
  #
  # TODO: Rising edge should be modeled as an object (or strucutre or class
  # or a container depending on language in use)


  ## Handle sign errors of raw data. This can happen due to installation
  #  problem.

  colnames(data) <- c("ts", "ap", "rp")

  if (mean(data$ap) < 0) {
    data$ap <- -data$ap
  }
  if (mean(data$rp) < 0) {
    data$rp <- -data$rp
  }

  ## -------------------------------------------------------------
  # Define parameters
  SAMPLING_RATE <- 15 # 15Hz sampling
  UNIT_TIME <- 3600 # 1.0 in dataTs correspondes to 3600
  # s <- sprintf('SAMPLING_RATE = %d',SAMPLING_RATE)
  # print(s)
  # s <- sprintf('UNIT_TIME = %d seconds corredponds to 1.0 in dataTs',UNIT_TIME)
  # print(s)


  # Slope threshold to declare FLAT
  FLAT_SLOPE_THRESH <- 1 # Per sample (trainslates to 15W/sec)
  # Slope thrshold to declare
  RF_SLOPE_THRESH <- 7 # Per sample (translates to 105W/sec)

  # Min length to be flatSection
  FLAT_LENGTH_THRESH <- 7 # At least this number of samples

  # Max number of missing samples for evaluation
  MAX_MISS_SAMPLES_SECS <- 1 # Up to 1 sec of missing samples are ignored.
  # Max time allowed for 'consecutive missing samples' to be valid. Scale properly
  CONS_MISS_SAMPLES_MAXTIME_THRESH <- MAX_MISS_SAMPLES_SECS/UNIT_TIME


  # Impulse removal related
  MAX_IMP_LENGTH_SECS <- 20 # Upto 20seconds for an impulse
  IMP_TIME_THRESH <- MAX_IMP_LENGTH_SECS/UNIT_TIME

  SIM_POWER_LEVEL_MARGIN <- 3.0 # Similar power level margin
  POWER_DIFF_THRESH <- 10.0 # At least this much needed to declare an appliance on or off event

  # Transient classification related
  MIN_MAX_AP_THRESH <- 7.0 # At least this much of fluctuation needed to be considered as a transient
  GEN_EDGE_LENGTH_THRESH <- 3 # The length of generalized-edge needs to be at least as long as this value

  FR_FLAT_ADD_RATIO_THRESH <- 0.1 # Ratio of difference between left and right peaks of FR compared to impulse's height

  # tag types
  FLAT <- 0
  RISING <- 1
  FALLING <- -1
  UNKNOWN <- -99

  # Diff
  diffVecAp <- diff(data$ap, lag = 1)
  diffVecRp <- diff(data$rp, lag = 1)
  diffVecTs <- diff(data$ts, lag = 1)
  diffVecAp[length(diffVecAp)+1] <- 0 # zero-pad the last sample
  diffVecRp[length(diffVecRp)+1] <- 0 # zero-pad the last sample
  diffVecTs[length(diffVecTs)+1] <- 1/SAMPLING_RATE # Pad with a number
  ## -------------------------------------------------------------
  DATA_LENGTH <- dim(data)[1]
  ## Per sample tags
  tagsAp <- rep(1, DATA_LENGTH) * UNKNOWN
  tagsRp <- rep(1, DATA_LENGTH) * UNKNOWN



  tagsAp[abs(diffVecAp) <= FLAT_SLOPE_THRESH] <- FLAT
  tagsAp[diffVecAp <= -1*RF_SLOPE_THRESH] <- FALLING
  tagsAp[RF_SLOPE_THRESH <= diffVecAp] <- RISING

  tagsRp[abs(diffVecRp) <= FLAT_SLOPE_THRESH] <- FLAT
  tagsRp[diffVecRp <= (-1) * RF_SLOPE_THRESH] <- FALLING
  tagsRp[RF_SLOPE_THRESH <= diffVecRp] <- RISING


  # Deal with missing data samples
  tagsAp[CONS_MISS_SAMPLES_MAXTIME_THRESH < diffVecTs] <- UNKNOWN
  tagsRp[CONS_MISS_SAMPLES_MAXTIME_THRESH < diffVecTs] <- UNKNOWN


  ## -------------------------------------------------------------
  # Detect stable/absolutely FLAT sections

  # Mentality : try to detect short sections of stable FLAT's
  lengthTags <- length(tagsAp)

  # if length(tagsRp ~= lengthTags)
  #   error('Lengths of AP and RP do not match. Pre-processing error.')
  # end

  # Initialize
  NUM_ATTRIBUTES <- 4
  flatSectionList <- matrix(0,lengthTags, NUM_ATTRIBUTES)
  numFlatSections <- 0



  k <- 1
  while(1)
  {
    if (tagsAp[k]==FLAT && tagsRp[k]==FLAT) # Use AND rule: both AP and RP need to be FLAT
    {

      flatStartIndex <- k
      k <- k+1
      while (tagsAp[k] == FLAT && tagsRp[k] == FLAT)
      {
        k <- k+1
        if (k >= lengthTags)
        {
          break
        }
      }
      flatEndIndex <- k
      flatLength <- flatEndIndex - flatStartIndex
      flatStartPowerAp <- data$ap[flatStartIndex]

      # Chunk length only. Ignore power deviation as the flat section becomes large - very low slope everywhere, anyway
      if (flatLength >= FLAT_LENGTH_THRESH)
      {
        numFlatSections <- numFlatSections+1
        flatSectionList[numFlatSections,] <- c(flatStartIndex,flatEndIndex,flatLength,flatStartPowerAp)
      }
    }

    k <- k+1

    # Check if all done
    if (k >= length(data$ap))
    {
      break
    }
  }

  # Truncate
  flatSectionList = flatSectionList[1:numFlatSections,]

  ## -------------------------------------------------------------
  # Flat removal
  # Remove flat sections that are part of impulse

  keepList <- vector() # Keep the first flatSection
  keepList[1] <- 1

  if(numFlatSections>2){# 12.01 update : flatSectionList가 비어?????? 경우 handling.
    for (k in 2:(nrow(flatSectionList)-1))
    {

      # condition 1 : Time-wise, prev and next flatSections are close enough
      # condition 2 : Power levels of prev and next flatSections are similar
      # condition 3 : Middle flatSection needs to be higher power

      if ((data$ts[flatSectionList[k+1,1]] - data$ts[flatSectionList[k-1,2]])< IMP_TIME_THRESH
          && abs(data$ap[flatSectionList[k+1,1]] - data$ap[flatSectionList[k-1,2]]) < SIM_POWER_LEVEL_MARGIN
          && POWER_DIFF_THRESH <= (min(data$ap[flatSectionList[k,1]:flatSectionList[k,2]]) - data$ap[flatSectionList[k-1,2]]))
      {
        # Do nothing. This will cause k'th flat section to be removed
      }
      else
      {
        keepList <- rbind(keepList, k)
      }
    }
  }

  keepList <- rbind(keepList, nrow(flatSectionList))
  # Remove flats inside of impulse shapes
  if (!is.null(nrow(flatSectionList))){
    flatSectionList <- flatSectionList[keepList,]
  }
  ####12.01 update start
  #keepList <- keepList[keepList>0]
  #if (nrow(flatSectionList) > 0){
  #  flatSectionList = unique(flatSectionList[which(keepList>0),])
  #}
  ####12.01 update end


  ## -------------------------------------------------------------
  # Flat addition
  # Add flat sections when between

  newFlat <- vector()
  if (!is.null(nrow(flatSectionList))){
    for (k in 1:(nrow(flatSectionList)-1))
    {
      # A genEdge
      startIndex <- flatSectionList[k,2] # [flatStartIndex flatEndIndex flatLength flatStartPowerAp]
      endIndex <- flatSectionList[k+1,1]
      genEdgeLength <- endIndex - startIndex + 1 # Transient

      # Vector to investigate
      dataApLocal <- data$ap[startIndex:endIndex]
      diffVecTsLocal <- diffVecTs[startIndex:endIndex-1]

      # Features
      maxAp <- max(dataApLocal)
      minAp <- min(dataApLocal)

      # Ignore if FLAT-like, too short, or too many missing data points
      if ((maxAp - minAp) < MIN_MAX_AP_THRESH || genEdgeLength < GEN_EDGE_LENGTH_THRESH)
      {
        next
      }

      # Create a summary table of FALLING and RISING streaks
      streakList <- data.frame()
      j <- startIndex

      while(1)
      {
        type <- tagsAp[j]
        sIndex <- j
        j <- j+1
        while (tagsAp[j] == type)
        {
          j <- j+1
          if(j>=endIndex)
          {
            break
          }
        }
        eIndex <- j

        if (type==FALLING || type == RISING)
        {
          streakList <- rbind(streakList, c(type, sIndex, eIndex, data$ap[sIndex], data$ap[eIndex]))
        }
        if (j >= endIndex)
        {
          break
        }
      }

      # Merge consecutive FALLING's and consecutive Rising's
      streakListMerged <- data.frame()
      if (2 <= nrow(streakList))
      {
        n <- 1
        while (1)
        {
          prevVec <- streakList[n,]
          n <- n+1
          if (nrow(streakList) < n)
          {
            break
          }
          while (prevVec[1] == streakList[n,1])
          {
            # merge
            prevVec[3] <- streakList[n,3]
            prevVec[5] <- streakList[n,5]
            n <- n+1
            if (nrow(streakList) < n)
            {
              break
            }
          }
          streakListMerged <- rbind(streakListMerged, prevVec)
          if (nrow(streakList) < n)
          {
            break
          }
        }
      }
      else
      {
        streakListMerged <- streakList
      }

      # Check if FALLING -> RISING happens
      # If a proper condition is met, add the min value point to the newFlat
      if (nrow(streakListMerged) > 1)
      {
        for (n in 1:(nrow(streakListMerged)-1))
        {
          if ((streakListMerged[n,1] == FALLING) && (streakListMerged[n+1,1] == RISING)
              && ((abs(streakListMerged[n,4] - streakListMerged[n+1,5]) <= SIM_POWER_LEVEL_MARGIN)
                  || (abs(streakListMerged[n,4] - streakListMerged[n+1,5]) <= (streakListMerged[n,4] - streakListMerged[n,5])*FR_FLAT_ADD_RATIO_THRESH)))
          {
            ind <- which.min(data$ap[streakListMerged[n,3]:streakListMerged[n+1,2]])
            newFlat <- rbind(newFlat, streakListMerged[n,3] + ind - 1)
          }
        }
      }
    }
  } else{#### 12.01 update start
    # flat section??? 비어?????? 경우 handling : 처음??? ?????? flat?????? ??????.
    # A genEdge
    startIndex <- 1 # [flatStartIndex flatEndIndex flatLength flatStartPowerAp]
    endIndex <- length(data$ts)
    genEdgeLength <- endIndex - startIndex + 1 # Transient

    # Vector to investigate
    dataApLocal <- data$ap[startIndex:endIndex]
    diffVecTsLocal <- diffVecTs[startIndex:(endIndex-1)]

    # Features
    maxAp <- max(dataApLocal)
    minAp <- min(dataApLocal)

    # Ignore if FLAT-like, too short, or too many missing data points
    if ((maxAp - minAp) < MIN_MAX_AP_THRESH || genEdgeLength < GEN_EDGE_LENGTH_THRESH){}

    # Create a summary table of FALLING and RISING streaks
    streakList <- data.frame()
    j <- startIndex

    while(1)
    {
      type <- tagsAp[j]
      sIndex <- j
      j <- j+1
      while (tagsAp[j] == type)
      {
        j <- j+1
        if(j>=endIndex)
        {
          break
        }
      }
      eIndex <- j

      if (type==FALLING || type == RISING)
      {
        streakList <- rbind(streakList, c(type, sIndex, eIndex, data$ap[sIndex], data$ap[eIndex]))
      }
      if (j >= endIndex)
      {
        break
      }
    }

    # Merge consecutive FALLING's and consecutive Rising's
    streakListMerged <- data.frame()
    if (2 <= nrow(streakList))
    {
      n <- 1
      while (1)
      {
        prevVec <- streakList[n,]
        n <- n+1
        if (nrow(streakList) < n)
        {
          break
        }
        while (prevVec[1] == streakList[n,1])
        {
          # merge
          prevVec[3] <- streakList[n,3]
          prevVec[5] <- streakList[n,5]
          n <- n+1
          if (nrow(streakList) < n)
          {
            break
          }
        }
        streakListMerged <- rbind(streakListMerged, prevVec)
        if (nrow(streakList) < n)
        {
          break
        }
      }
    }
    else
    {
      streakListMerged <- streakList
    }

    # Check if FALLING -> RISING happens
    # If a proper condition is met, add the min value point to the newFlat
    if (nrow(streakListMerged) > 1)
    {
      for (n in 1:(nrow(streakListMerged)-1))
      {
        if ((streakListMerged[n,1] == FALLING) && (streakListMerged[n+1,1] == RISING)
            && ((abs(streakListMerged[n,4] - streakListMerged[n+1,5]) <= SIM_POWER_LEVEL_MARGIN)
                || (abs(streakListMerged[n,4] - streakListMerged[n+1,5]) <= (streakListMerged[n,4] - streakListMerged[n,5])*FR_FLAT_ADD_RATIO_THRESH)))
        {
          ind <- which.min(data$ap[streakListMerged[n,3]:streakListMerged[n+1,2]])
          newFlat <- rbind(newFlat, streakListMerged[n,3] + ind - 1)
        }
      }
    }
  }
  ### 12.01 update end

  # Add the new flats (point-flats) to the list
  flatSectionListNew <- cbind(newFlat, newFlat, rep(1, length(newFlat)), data$ap[newFlat])
  flatSectionList <- rbind(flatSectionList, flatSectionListNew) # Merge
  flatSectionList <- flatSectionList[order(flatSectionList[,1]),]



  ## -------------------------------------------------------------
  # Classify generalized-edges between FLAT's


  # Initialize
  NUM_ATTRIBUTES <- 15


  if (!is.null(nrow(flatSectionList))){ # 12.01 update
    genEdgeList <- matrix(0, nrow(flatSectionList)-1, NUM_ATTRIBUTES )
    numEdges <- 0
    genEdgeStringList <- rep('', nrow(flatSectionList)-1)
    for (k in 1:(nrow(flatSectionList)-1))
    {
      previousFlatSection <- flatSectionList[k,] # [flatStartIndex flatEndIndex flatLength flatStartPowerAp]
      nextFlatSection <- flatSectionList[k+1,]

      startIndex <- previousFlatSection[2]
      endIndex <- nextFlatSection[1]
      genEdgeLength <- endIndex - startIndex + 1

      # Vector to investigate
      dataApLocal <- data$ap[startIndex:endIndex]
      diffVecTsLocal <- diffVecTs[startIndex:(endIndex-1)]



      # Features
      maxAp <- max(dataApLocal)
      minAp <- min(dataApLocal)

      # Ignore if FLAT-like, too short, or too many missing data points
      if ((maxAp - minAp) < MIN_MAX_AP_THRESH || genEdgeLength < GEN_EDGE_LENGTH_THRESH)
      {
        next
      }

      # Vector to investigate: additional
      diffVecApLocal <- diffVecAp[startIndex:(endIndex-1)]

      dataRpLocal <- data$rp[startIndex:endIndex]
      diffVecRpLocal <- diffVecRp[startIndex:(endIndex-1)]

      tagsApLocal <- tagsAp[startIndex:(endIndex-1)]
      tagsRpLocal <- tagsRp[startIndex:(endIndex-1)]

      #browser()

      # Features
      startAp <- dataApLocal[1]
      endAp <- dataApLocal[length(dataApLocal)]

      startRp <- dataRpLocal[1]
      endRp <- dataRpLocal[length(dataRpLocal)]
      maxRp <- max(dataRpLocal)
      minRp <- min(dataRpLocal)

      endStartAp <- endAp - startAp
      maxStartAp <- maxAp - startAp
      minStartAp <- minAp - startAp

      endStartRp <- endRp - startRp
      maxStartRp <- maxRp - startRp
      minStartRp <- minRp - startRp

      maxSlopeAp <- max(diffVecApLocal)
      minSlopeAp <- min(diffVecApLocal)

      maxSlopeRp <- max(diffVecRpLocal)
      minSlopeRp <- min(diffVecRpLocal)

      numTransAp <- sum(tagsApLocal[1:(length(tagsApLocal)-1)] != tagsApLocal[2:length(tagsApLocal)])
      numTransRp <- sum(tagsRpLocal[1:(length(tagsRpLocal)-1)] != tagsRpLocal[2:length(tagsRpLocal)])

      numSegmentsAp <- 0
      shapeAlphabet <- ''
      type <- 'U'
      for (j in 1:length(tagsApLocal))
      {
        if (tagsApLocal[j] == RISING)
        {
          if (type != 'R')
          {
            shapeAlphabet <- paste(shapeAlphabet, 'R', sep = "")
            numSegmentsAp <- numSegmentsAp + 1
            type <- 'R'
          }
        }
        else if (tagsApLocal[j] == FALLING)
        {
          if (type != 'F')
          {
            shapeAlphabet <- paste(shapeAlphabet, 'F', sep = "")
            numSegmentsAp <- numSegmentsAp + 1
            type <- 'F'
          }
        }
        else
        {
          type <- 'U'
        }


      }

      type <- -1



      if(strcmp(shapeAlphabet, 'R') && POWER_DIFF_THRESH < endStartAp)
      {
        # Simple R
        type <- 11
      }
      else if(strcmp(shapeAlphabet, 'RF') && POWER_DIFF_THRESH < endStartAp)
      {
        # 'R' with overshooting
        type <- 12
      }
      else if(strcmp(shapeAlphabet, 'F') && endStartAp < (-1)*POWER_DIFF_THRESH)
      {
        # Simple 'F'
        type <- 21
      }
      else if(strcmp(shapeAlphabet, 'RF') && abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN)
      {
        # Simple impulse
        type <- 31
      }
      else if(strcmp(shapeAlphabet, 'RFF') && abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN)
      {
        # RFF impulse
        type <- 32
      }
      else if(substr(shapeAlphabet,1,1)=='F' && abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN) #12.01 update
      {
        # falling impulse
        type <- 33
      }
      else if(abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN)
      {
        # Other impulses
        type <- 34
      }



      edgeInfo <- c(type, startIndex, endIndex, endStartAp, endStartRp, maxStartAp, maxStartRp, minStartAp, minStartRp,
                    numTransAp, numSegmentsAp, maxSlopeAp, maxSlopeRp, minSlopeAp, minSlopeRp)


      # Now update genEdgeList
      numEdges <- numEdges + 1
      genEdgeList[numEdges,] <- edgeInfo
      genEdgeStringList[numEdges] <- shapeAlphabet
    }
  } else{ #### 12.01 update start
    # flatSectionList가 비어?????? 경우 handling : 처음??? ?????? flat?????? ??????.

    genEdgeList <- data.frame()
    numEdges <- 0
    genEdgeStringList <- data.frame()

    startIndex <- 1
    endIndex <- length(data$ts)
    genEdgeLength <- endIndex - startIndex + 1

    # Vector to investigate
    dataApLocal <- data$ap[startIndex:endIndex]
    diffVecTsLocal <- diffVecTs[startIndex:(endIndex-1)]



    # Features
    maxAp <- max(dataApLocal)
    minAp <- min(dataApLocal)

    # Ignore if FLAT-like, too short, or too many missing data points
    if (!((maxAp - minAp) < MIN_MAX_AP_THRESH)&&!(genEdgeLength < GEN_EDGE_LENGTH_THRESH)){

      # Vector to investigate: additional
      diffVecApLocal <- diffVecAp[startIndex:(endIndex-1)]

      dataRpLocal <- data$rp[startIndex:endIndex]
      diffVecRpLocal <- diffVecRp[startIndex:(endIndex-1)]

      tagsApLocal <- tagsAp[startIndex:(endIndex-1)]
      tagsRpLocal <- tagsRp[startIndex:(endIndex-1)]

      #browser()

      # Features
      startAp <- dataApLocal[1]
      endAp <- dataApLocal[length(dataApLocal)]

      startRp <- dataRpLocal[1]
      endRp <- dataRpLocal[length(dataRpLocal)]
      maxRp <- max(dataRpLocal)
      minRp <- min(dataRpLocal)

      endStartAp <- endAp - startAp
      maxStartAp <- maxAp - startAp
      minStartAp <- minAp - startAp

      endStartRp <- endRp - startRp
      maxStartRp <- maxRp - startRp
      minStartRp <- minRp - startRp

      maxSlopeAp <- max(diffVecApLocal)
      minSlopeAp <- min(diffVecApLocal)

      maxSlopeRp <- max(diffVecRpLocal)
      minSlopeRp <- min(diffVecRpLocal)

      numTransAp <- sum(tagsApLocal[1:(length(tagsApLocal)-1)] != tagsApLocal[2:length(tagsApLocal)])
      numTransRp <- sum(tagsRpLocal[1:(length(tagsRpLocal)-1)] != tagsRpLocal[2:length(tagsRpLocal)])

      numSegmentsAp <- 0
      shapeAlphabet <- ''
      type <- 'U'
      for (j in 1:length(tagsApLocal))
      {
        if (tagsApLocal[j] == RISING)
        {
          if (type != 'R')
          {
            shapeAlphabet <- paste(shapeAlphabet, 'R', sep = "")
            numSegmentsAp <- numSegmentsAp + 1
            type <- 'R'
          }
        }
        else if (tagsApLocal[j] == FALLING)
        {
          if (type != 'F')
          {
            shapeAlphabet <- paste(shapeAlphabet, 'F', sep = "")
            numSegmentsAp <- numSegmentsAp + 1
            type <- 'F'
          }
        }
        else
        {
          type <- 'U'
        }


      }

      type <- -1


      if(strcmp(shapeAlphabet, 'R') && POWER_DIFF_THRESH < endStartAp)
      {
        # Simple R
        type <- 11
      }
      else if(strcmp(shapeAlphabet, 'RF') && POWER_DIFF_THRESH < endStartAp)
      {
        # 'R' with overshooting
        type <- 12
      }
      else if(strcmp(shapeAlphabet, 'F') && endStartAp < (-1)*POWER_DIFF_THRESH)
      {
        # Simple 'F'
        type <- 21
      }
      else if(strcmp(shapeAlphabet, 'RF') && abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN)
      {
        # Simple impulse
        type <- 31
      }
      else if(strcmp(shapeAlphabet, 'RFF') && abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN)
      {
        # RFF impulse
        type <- 32
      }
      else if(substr(shapeAlphabet,1,1)=='F' && abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN) #12.01 update
      {
        # falling impulse
        type <- 33
      }
      else if(abs(endStartAp)<=SIM_POWER_LEVEL_MARGIN)
      {
        # Other impulses
        type <- 34
      }


      edgeInfo <- c(type, startIndex, endIndex, endStartAp, endStartRp, maxStartAp, maxStartRp, minStartAp, minStartRp,
                    numTransAp, numSegmentsAp, maxSlopeAp, maxSlopeRp, minSlopeAp, minSlopeRp)


      # Now update genEdgeList
      numEdges <- numEdges + 1
      genEdgeList <- edgeInfo
      genEdgeStringList <-rbind(genEdgeStringList, shapeAlphabet)
    }
  }#### 12.01 update end

  # Truncate
  if ((!is.null(nrow(genEdgeList)))&&numEdges>0){
    genEdgeList <- as.data.frame(genEdgeList[1:numEdges,])
    if (numEdges==1){
      genEdgeList <- t(genEdgeList)
    }
    colnames(genEdgeList) <- c("type", "startIndex", "endIndex", "endStartAp", "endStartRp", "maxStartAp", "maxStartRp", "minStartAp", "minStartRp",
                               "numTransAp", "numSegmentsAp", "maxSlopeAp", "maxSlopeRp", "minSlopeAp", "minSlopeRp")
    genEdgeStringList <- as.data.frame(genEdgeStringList[1:numEdges])
    colnames(genEdgeStringList) <- "edgeString"
  } else if(numEdges==1){
    genEdgeList <- as.data.frame(t(genEdgeList))
    colnames(genEdgeList) <- c("type", "startIndex", "endIndex", "endStartAp", "endStartRp", "maxStartAp", "maxStartRp", "minStartAp", "minStartRp",
                               "numTransAp", "numSegmentsAp", "maxSlopeAp", "maxSlopeRp", "minSlopeAp", "minSlopeRp")
    genEdgeStringList <- as.data.frame(genEdgeStringList[1:numEdges])
    colnames(genEdgeStringList) <- "edgeString"
  } else{
    genEdgeList <- data.frame()
    genEdgeStringList <- data.frame()
  }


  # 12.01 update : flatSectionList error handling
  if (is.null(nrow(flatSectionList))){
    flatSectionList <- data.frame(t(flatSectionList))
  }

  if (flatSectionList[1,1]==0&&flatSectionList[1,2]==0&&flatSectionList[1,3]==0&&flatSectionList[1,1]==0){
    flatSectionList = flatSectionList[-1,]
  }

  if (is.null(nrow(flatSectionList))){
    flatSectionList <- data.frame(t(flatSectionList))
  }

  colnames(flatSectionList) <- c("flatStartIndex","flatEndIndex", "flatLength", "flatStartPowerAp")



  return(list("edgeList"=genEdgeList, "edgeStringList" = genEdgeStringList, "flatList" = flatSectionList))

  # 12.01 delete
  #edgeInfoTable <- cbind(genEdgeList, genEdgeStringList, flatSectionList)
  #colnames(edgeInfoTable) <- c("type", "startIndex", "endIndex", "endStartAp", "endStartRp", "maxStartAp", "maxStartRp", "minStartAp", "minStartRp",
  #                             "numTransAp", "numSegmentsAp", "maxSlopeAp", "maxSlopeRp", "minSlopeAp", "minSlopeRp", "edgeSequence", "flatStartIndex",
  #                             "flatEndIndex", "flatLength", "flatStartPowerAp")

  #return(edgeInfoTable)

}



######################################################################################################################


#' @title tvFeatureExtractor
#' @author CH
#' @description extract tv block features from total data
#' @update 2015.12.03
#' @input data : c(timestamp, active_power, reactive_power)
#'        blockSize : time block Size
#' @output features of data blocks


# tv feature extractor speed test version without envelope

tvFeatureExtractor = function(data, blockSize){

  ENVELOPE_SIDE_LENGTH = 15*3 # 3 seconds
  QUANTILE_TYPE = 5
  # blockSize = 15*60*5;
  total_ts = data[,1]
  total_ap = data[,2]
  total_rp = data[,3]

  totalApBlock = cutBlocks(total_ap, blockSize)
  totalRpBlock = cutBlocks(total_rp, blockSize)
  totalTsBlock = cutBlocks(total_ts, blockSize)

  # edge features
  type11Count = rep(0,nrow(totalApBlock))
  type12Count = rep(0,nrow(totalApBlock))
  type21Count = rep(0,nrow(totalApBlock))
  type31Count = rep(0,nrow(totalApBlock))
  type32Count = rep(0,nrow(totalApBlock))
  type33Count = rep(0,nrow(totalApBlock))
  type34Count = rep(0,nrow(totalApBlock))

  edgeCount = rep(0,nrow(totalApBlock))
  maxSlopeMean = rep(0,nrow(totalApBlock))
  maxSlopeVar = rep(0,nrow(totalApBlock))
  minSlopeMean = rep(0,nrow(totalApBlock))
  minSlopeVar = rep(0,nrow(totalApBlock))

  stringLengthMax = rep(0,nrow(totalApBlock))
  stringLengthMean = rep(0,nrow(totalApBlock))
  stringLengthVar = rep(0,nrow(totalApBlock))

  # additional edge features
  smallRisingCount = rep(0,nrow(totalApBlock))
  smallFallingCount = rep(0,nrow(totalApBlock))
  largeImpulseCount = rep(0,nrow(totalApBlock))

  # ap features
  Q1MinDiff = rep(0,nrow(totalApBlock))

  # 7 diff(0.5 second diff) features
  diff7Mean = rep(0,nrow(totalApBlock))
  diff7Var = rep(0,nrow(totalApBlock))
  diff7SmallCount = rep(0,nrow(totalApBlock))
  diff7LargeCount = rep(0,nrow(totalApBlock))

  # 30 diff (2 second diff) features
  diff30Mean = rep(0,nrow(totalApBlock))
  diff30Var = rep(0,nrow(totalApBlock))
  diff30SmallCount = rep(0,nrow(totalApBlock))
  diff30LargeCount = rep(0,nrow(totalApBlock))

  # flat features
  flatDiffVar = rep(0,nrow(totalApBlock))
  flatLengthVar = rep(0,nrow(totalApBlock))
  flatCount = rep(0,nrow(totalApBlock))

  #   # 01.07 added features
  #   # flat features
  flatLengthMax = rep(0,nrow(totalApBlock))
  flatLengthMin = rep(0,nrow(totalApBlock))
  flatLengthRange = rep(0,nrow(totalApBlock))
  flatLengthMedian = rep(0,nrow(totalApBlock))

  # edge features
  edgeMaxLength = rep(0,nrow(totalApBlock))
  emptyStringCount = rep(0,nrow(totalApBlock))
  emptyStringRatio = rep(0,nrow(totalApBlock))

  # between edge features
  betweenEdgeMax = rep(0,nrow(totalApBlock))
  betweenEdgeMin = rep(0,nrow(totalApBlock))
  betweenEdgeMean = rep(0,nrow(totalApBlock))
  betweenEdgeVar = rep(0,nrow(totalApBlock))

  # between flat features
  betweenFlatMax = rep(0,nrow(totalApBlock))
  betweenFlatMin = rep(0,nrow(totalApBlock))
  betweenFlatMean = rep(0,nrow(totalApBlock))
  betweenFlatVar = rep(0,nrow(totalApBlock))


  ## 01.11 added features : features related to envelopes
  # upper - lower
  #   envelopeUpperLowerIntegral = rep(0,nrow(totalApBlock))
  #   envelopeUpperLowerMax = rep(0,nrow(totalApBlock))
  #   envelopeUpperLowerMin = rep(0,nrow(totalApBlock))
  #   envelopeUpperLower1Q = rep(0,nrow(totalApBlock))
  #   envelopeUpperLowerMedian = rep(0,nrow(totalApBlock))
  #   envelopeUpperLower3Q = rep(0,nrow(totalApBlock))
  #   envelopeUpperLowerSmallCount = rep(0,nrow(totalApBlock))
  #
  # upper - original
  #   envelopeUpperOriginalIntegral = rep(0,nrow(totalApBlock))
  #   envelopeUpperOriginalMax = rep(0,nrow(totalApBlock))
  #   envelopeUpperOriginalMin = rep(0,nrow(totalApBlock))
  #   envelopeUpperOriginal1Q = rep(0,nrow(totalApBlock))
  #   envelopeUpperOriginalMedian = rep(0,nrow(totalApBlock))
  #   envelopeUpperOriginal3Q = rep(0,nrow(totalApBlock))
  #   envelopeUpperOriginalRainCount = rep(0,nrow(totalApBlock))

  #   # original - lower
  #   envelopeOriginalLowerIntegral = rep(0,nrow(totalApBlock))
  #   envelopeOriginalLowerMax = rep(0,nrow(totalApBlock))
  #   envelopeOriginalLowerMin = rep(0,nrow(totalApBlock))
  #   envelopeOriginalLower1Q = rep(0,nrow(totalApBlock))
  #   envelopeOriginalLowerMedian = rep(0,nrow(totalApBlock))
  #   envelopeOriginalLower3Q = rep(0,nrow(totalApBlock))

  ## 01.12 added features
  # diff7
  diff7VerySmallCount = zeros(size(totalApBlock, 1), 1);
  diff7Small2Count = zeros(size(totalApBlock, 1), 1);

  # diff30
  diff30VerySmallCount = zeros(size(totalApBlock, 1), 1);
  diff30Small2Count = zeros(size(totalApBlock, 1), 1);

  # diff7 bin features
  diff7BinCountMax = zeros(size(totalApBlock, 1), 1);
  diff7BinCountQ3 = zeros(size(totalApBlock, 1), 1);
  diff7BinCountMedian = zeros(size(totalApBlock, 1), 1);
  diff7BinCountQ1 = zeros(size(totalApBlock, 1), 1);
  diff7BinCountEntropy = zeros(size(totalApBlock, 1), 1);

  # diff30 bin features
  diff30BinCountMax = zeros(size(totalApBlock, 1), 1);
  diff30BinCountQ3 = zeros(size(totalApBlock, 1), 1);
  diff30BinCountMedian = zeros(size(totalApBlock, 1), 1);
  diff30BinCountQ1 = zeros(size(totalApBlock, 1), 1);
  diff30BinCountEntropy = zeros(size(totalApBlock, 1), 1);



  ###################################### feature extraction iteration #####################################
  for (i in 1:nrow(totalApBlock)){
    blockTs = totalTsBlock[i,]
    blockAp = totalApBlock[i,]
    blockRp = totalRpBlock[i,]

    dataBlock = data.frame(cbind(blockTs, blockAp, blockRp))
    colnames(dataBlock) = c("ts","ap","rp")

    ## edge features
    edgeDetectorResult = edgeDetectorTV(dataBlock)
    edgeList = edgeDetectorResult[[1]]
    edgeStringList = edgeDetectorResult[[2]]
    flatList = edgeDetectorResult[[3]]

    if (nrow(edgeStringList)>1){
      stringLengthMax[i] = max(apply(edgeStringList, MARGIN = 1, FUN = nchar))
      stringLengthMean[i] = mean(apply(edgeStringList, MARGIN=1, FUN = nchar))
      stringLengthVar[i] = var(apply(edgeStringList, MARGIN=1, FUN = nchar))
    } else if(nrow(edgeStringList)==1){
      stringLengthMax[i] = apply(edgeStringList, MARGIN = 1, FUN = nchar)
      stringLengthMean[i] = apply(edgeStringList, MARGIN = 1, FUN = nchar)
      stringLengthVar[i] = 0
    } else{
      stringLengthMax[i] = 0
      stringLengthMean[i] = 0
      stringLengthVar[i] = 0
    }

    ## edge statistics features
    # type count
    if (nrow(edgeList)>0){
      type11Count[i] = sum(edgeList[,1]==11)
      type12Count[i] = sum(edgeList[,1]==12)
      type21Count[i] = sum(edgeList[,1]==21)
      type31Count[i] = sum(edgeList[,1]==31)
      type32Count[i] = sum(edgeList[,1]==32)
      type33Count[i] = sum(edgeList[,1]==33)
      type34Count[i] = sum(edgeList[,1]==34)
    } else{
      type11Count[i] = 0
      type12Count[i] = 0
      type21Count[i] = 0
      type31Count[i] = 0
      type32Count[i] = 0
      type33Count[i] = 0
      type34Count[i] = 0
    }

    edgeCount[i] = nrow(edgeList)

    if (edgeCount[i] > 1){
      maxSlopeMean[i] = mean(edgeList[,6])
      maxSlopeVar[i] = var(edgeList[,6])
      minSlopeMean[i] = mean(edgeList[,8])
      minSlopeVar[i] = var(edgeList[,8])
      smallRisingCount[i] = sum(edgeList[,6]>0 & edgeList[,6]<40)
      smallFallingCount[i] = sum(edgeList[,8]<0 & edgeList[,8]>-40)
      largeImpulseCount[i] = sum(edgeList[,6]>100 & (edgeList[,1]== 31 | edgeList[,1] == 32 | edgeList[,1] == 33 | edgeList[,1] == 34))

      # 01.07 added
      edgeMaxLength[i] = max(edgeList[,3]-edgeList[,2])
      emptyStringCount[i] = countEmptyString(edgeStringList)
      emptyStringRatio[i] = emptyStringCount[i]/edgeCount[i]

    } else if(edgeCount[i] == 1){
      maxSlopeMean[i] = edgeList[,6]
      maxSlopeVar[i] = 0
      minSlopeMean[i] = edgeList[,8]
      minSlopeVar[i] = 0
      smallRisingCount[i] = sum(edgeList[,6]>0 & edgeList[,6]<40)
      smallFallingCount[i] = sum(edgeList[,8]<0 & edgeList[,8]>-40)
      largeImpulseCount[i] = sum(edgeList[,6]>100 & (edgeList[,1]== 31 | edgeList[,1] == 32 | edgeList[,1] == 33 | edgeList[,1] == 34))

      # 01.07 added
      edgeMaxLength[i] = edgeList[,3]-edgeList[,2]
      emptyStringCount[i] = countEmptyString(edgeStringList)
      emptyStringRatio[i] = emptyStringCount[i]/edgeCount[i]
    } else {
      maxSlopeMean[i] = 0
      maxSlopeVar[i] = 0
      minSlopeMean[i] = 0
      minSlopeVar[i] = 0
      smallRisingCount[i] = 0
      smallFallingCount[i] = 0
      largeImpulseCount[i] = 0

      # 01.07 added
      edgeMaxLength[i] = 0
      emptyStringCount[i] = 0
      emptyStringRatio[i] = 0
    }

    ## 01.07 added : betweenEdgeMax, betweenEdgeMin, betweenEdgeMean
    if (edgeCount[i] > 2){
      endIndex = edgeCount[i]
      frontEdgeEnd = edgeList[1:(endIndex-1), 3]
      backEdgeStart = edgeList[2:endIndex,2]
      betweenEdge = backEdgeStart-frontEdgeEnd
      betweenEdgeMax[i] = max(betweenEdge)
      betweenEdgeMin[i] = min(betweenEdge)
      betweenEdgeMean[i] = mean(betweenEdge)
      betweenEdgeVar[i] = var(betweenEdge)
    } else if(edgeCount[i] == 2){
      betweenEdgeMax[i] = edgeList[2,2]-edgeList[1,3]
      betweenEdgeMin[i] = edgeList[2,2]-edgeList[1,3]
      betweenEdgeMean[i] = edgeList[2,2]-edgeList[1,3]
      betweenEdgeVar[i] = 0
    } else {
      betweenEdgeMax[i] = 0
      betweenEdgeMin[i] = 0
      betweenEdgeMean[i] = 0
      betweenEdgeVar[i] = 0
    }



    ## flat features
    # 01.07 added : flatLengthMax, flatLengthMin, flatLengthRange, flatLengthMedian
    flatDiff = blockAp[flatList[,2]] - blockAp[flatList[,1]]
    if (nrow(flatList)>1){
      flatDiffVar[i] = var(flatDiff)
      flatLengthVar[i] = var(flatList[,3])
      flatCount[i] = nrow(flatList)
      flatLengthMax[i] = max(flatList[,3])
      flatLengthMin[i] = min(flatList[,3])
      flatLengthRange[i] = flatLengthMax[i]-flatLengthMin[i]
      flatLengthMedian[i] = median(flatList[,3])
    } else if(nrow(flatList)==1){
      flatDiffVar[i] = 0
      flatLengthVar[i] = 0
      flatCount[i] = 1
      flatLengthMax[i] = flatList[,3]
      flatLengthMin[i] = flatList[,3]
      flatLengthRange[i] = flatLengthMax[i]-flatLengthMin[i]
      flatLengthMedian[i] = flatList[,3]
    }else{
      flatDiffVar[i] = 0
      flatLengthVar[i] = 0
      flatCount[i] = 0
      flatLengthMax[i] = 0
      flatLengthMin[i] = 0
      flatLengthRange[i] = 0
      flatLengthMedian[i] = 0
    }

    # 01.07 added : betweenFlatMax, betweenFlatMin, betweenFlatMean

    if (nrow(flatList) > 2){
      endIndex = nrow(flatList)
      frontFlatEnd = flatList[1:(endIndex-1), 2]
      backFlatStart = flatList[2:endIndex,1]
      betweenFlat = backFlatStart-frontFlatEnd
      betweenFlatMax[i] = max(betweenFlat)
      betweenFlatMin[i] = min(betweenFlat)
      betweenFlatMean[i] = mean(betweenFlat)
      betweenFlatVar[i] = var(betweenFlat)
    } else if(nrow(flatList) == 2){
      betweenFlatMax[i] = flatList[2,1]-flatList[1,2]
      betweenFlatMin[i] = flatList[2,1]-flatList[1,2]
      betweenFlatMean[i] = flatList[2,1]-flatList[1,2]
      betweenFlatVar[i] = 0
    } else {
      betweenFlatMax[i] = 0
      betweenFlatMin[i] = 0
      betweenFlatMean[i] = 0
      betweenFlatVar[i] = 0
    }


    ## ap features
    Q1MinDiff[i] = quantile(blockAp, 0.25, type = QUANTILE_TYPE) - min(blockAp)

    # 7 samples diff features
    diff7 = blockAp[8:length(blockAp)] - blockAp[1:(length(blockAp)-7)]
    diff7Mean[i] = mean(diff7)
    diff7Var[i] = var(diff7)
    diff7SmallCount[i] = sum(abs(diff7)<10 & abs(diff7)>2)
    diff7LargeCount[i] = sum(abs(diff7)>30 & abs(diff7)<60)

    # 01.12 added features
    diff7VerySmallCount[i] = sum(abs(diff7)>=0 & abs(diff7)<1)
    diff7Small2Count[i] = sum(abs(diff7)>=10 & abs(diff7)<20)

    diff7BinCount = countBin(diff7)
    diff7BinCountMax[i] = max(diff7BinCount)
    diff7BinCountQ3[i] = quantile(diff7BinCount, 0.75, type = QUANTILE_TYPE)
    diff7BinCountMedian[i] = median(diff7BinCount)
    diff7BinCountQ1[i] = quantile(diff7BinCount, 0.25, type = QUANTILE_TYPE)
    diff7BinCountEntropy[i] = getEntropy(diff7BinCount)

    # 30 diff features
    diff30 = blockAp[31:length(blockAp)] - blockAp[1:(length(blockAp)-30)]
    diff30Mean[i] = mean(diff30)
    diff30Var[i] = var(diff30)
    diff30SmallCount[i] = sum(abs(diff30)<10 & abs(diff30)>2)
    diff30LargeCount[i] = sum(abs(diff30)>30 & abs(diff30)<60)

    # 01.12 added features
    diff30VerySmallCount[i] = sum(abs(diff30)>=0 & abs(diff30)<1)
    diff30Small2Count[i] = sum(abs(diff30)>=10 & abs(diff30)<20)

    diff30BinCount = countBin(diff30)
    diff30BinCountMax[i] = max(diff30BinCount)
    diff30BinCountQ3[i] = quantile(diff30BinCount, 0.75, type = QUANTILE_TYPE)
    diff30BinCountMedian[i] = median(diff30BinCount)
    diff30BinCountQ1[i] = quantile(diff30BinCount, 0.25, type = QUANTILE_TYPE)
    diff30BinCountEntropy[i] = getEntropy(diff30BinCount)

    ## envelope features
    #     blockData = data.frame(cbind(blockTs, blockAp, blockRp))
    #     colnames(blockData) = c('timestamp', 'active_power', 'reactive_power')
    #     envelope = envelopeDetector(blockData, ENVELOPE_SIDE_LENGTH, ENVELOPE_SIDE_LENGTH, FALSE, TRUE)
    #     lower = envelope$LowerEnvelope
    #     upper = envelope$UpperEnvelope
    #     upperLower = upper-lower
    #     upperOriginal = upper-blockAp
    #     originalLower = blockAp-lower

    #     # upper - lower
    #     envelopeUpperLowerIntegral[i] = sum(upperLower)
    #     envelopeUpperLowerMax[i] = max(upperLower)
    #     envelopeUpperLowerMin[i] = min(upperLower)
    #     envelopeUpperLower1Q[i] = quantile(upperLower, 0.25, type = QUANTILE_TYPE)
    #     envelopeUpperLowerMedian[i] = median(upperLower)
    #     envelopeUpperLower3Q[i] = quantile(upperLower, 0.75, type = QUANTILE_TYPE)
    #     envelopeUpperLowerSmallCount[i] = sum(upperLower>=3 & upperLower<=15)
    #
    #     # upper - original
    #     envelopeUpperOriginalIntegral[i] = sum(upperOriginal)
    #     envelopeUpperOriginalMax[i] = max(upperOriginal)
    #     envelopeUpperOriginalMin[i] = min(upperOriginal)
    #     envelopeUpperOriginal1Q[i] = quantile(upperOriginal, 0.25, type = QUANTILE_TYPE)
    #     envelopeUpperOriginalMedian[i] = median(upperOriginal)
    #     envelopeUpperOriginal3Q[i] = quantile(upperOriginal, 0.75, type = QUANTILE_TYPE)
    #     envelopeUpperOriginalRainCount[i] = sum(upperOriginal>=10&upperOriginal<=80)
    #
    #     # original - lower
    #     envelopeOriginalLowerIntegral[i] = sum(originalLower)
    #     envelopeOriginalLowerMax[i] = max(originalLower)
    #     envelopeOriginalLowerMin[i] = min(originalLower)
    #     envelopeOriginalLower1Q[i] = quantile(originalLower, 0.25, type = QUANTILE_TYPE)
    #     envelopeOriginalLowerMedian[i] = median(originalLower)
    #     envelopeOriginalLower3Q[i] = quantile(originalLower, 0.75, type = QUANTILE_TYPE)
    #
    sprintf('%d / %d', i, nrow(data))
  } # block loop end



  featureTable = cbind(type11Count, type12Count, type21Count, type31Count, type32Count,
                       type33Count,type34Count, edgeCount, maxSlopeMean, maxSlopeVar, minSlopeMean, minSlopeVar,
                       smallRisingCount, smallFallingCount, largeImpulseCount, Q1MinDiff,
                       stringLengthMax, stringLengthMean, stringLengthVar, diff7Mean, diff7Var, diff7SmallCount, diff7LargeCount,
                       diff30Mean, diff30Var, diff30SmallCount, diff30LargeCount, flatDiffVar, flatLengthVar, flatCount,
                       flatLengthMax, flatLengthMin, flatLengthRange, flatLengthMedian,
                       edgeMaxLength, emptyStringCount, emptyStringRatio,
                       betweenEdgeMax, betweenEdgeMin, betweenEdgeMean, betweenEdgeVar,
                       betweenFlatMax, betweenFlatMin, betweenFlatMean, betweenFlatVar,
                       #                        envelopeUpperLowerIntegral, envelopeUpperLowerMax, envelopeUpperLowerMin, envelopeUpperLower1Q,
                       #                        envelopeUpperLowerMedian, envelopeUpperLower3Q, envelopeUpperLowerSmallCount,
                       #                        envelopeUpperOriginalIntegral, envelopeUpperOriginalMax, envelopeUpperOriginalMin, envelopeUpperOriginal1Q,
                       #                        envelopeUpperOriginalMedian, envelopeUpperOriginal3Q, envelopeUpperOriginalRainCount,
                       #                        envelopeOriginalLowerIntegral, envelopeOriginalLowerMax, envelopeOriginalLowerMin, envelopeOriginalLower1Q,
                       #                        envelopeOriginalLowerMedian, envelopeOriginalLower3Q,
                       diff7VerySmallCount, diff7Small2Count, diff7BinCountMax, diff7BinCountQ3,
                       diff7BinCountMedian, diff7BinCountQ1, diff7BinCountEntropy,
                       diff30VerySmallCount, diff30Small2Count, diff30BinCountMax, diff30BinCountQ3,
                       diff30BinCountMedian, diff30BinCountQ1, diff30BinCountEntropy
  )

  return(featureTable)
}

#############################################################################################
#
# description : count edges that have empty string
countEmptyString = function(edgeStringList){
  return (sum(apply(edgeStringList, MARGIN = 1, FUN = function(x) x=='')))
}

#############################################################################################
#
# description : bin count of diff
# used in feature extraction related to bin counting
countBin = function(numberList){
  # bin Range Setup
  binKnots = c(0,1,5,10,15,20,25,30,40,50,60,70,80,90,100)
  binCount = rep(0, length(binKnots)-1)

  # binning
  for (i in 1:(length(binKnots)-1)){
    lower = binKnots[i]
    upper = binKnots[i+1]
    binCount[i] = sum(numberList>=lower & numberList<upper)
  }

  return(binCount)
}

#############################################################################################
#
# description : calculate entropy from bin count
# used in feature extraction related to bin counting
getEntropy = function(numberList){
  probList = numberList/sum(numberList)
  entropy = 0
  for (i in 1:length(probList)){
    if (probList[i]!=0){
      entropy = entropy - probList[i]*log(probList[i])
    }
  }
  return(entropy)
}




#' @title tvReconstructor
#' @author CH
#' @description reconstruct tv active power sequence using classified block state
#' @update 2015.12.03
#' @input blockStateList : predicted block on-off state
#'        data : total data c(timestamp, active_power, reactive power)
#'        blockSize
#'        tvApValue : active power of tv-on state
#' @output tv active_power vector
#'


######################################################################################################################


tvReconstructor <- function(blockStateList, data, blockSize, tvApValue)
{

  tv_ap = rep(0, length(data[,1]))

  # blockState on -> tvApValue
  for (i in 1:length(blockStateList)){
    if (blockStateList[i]==1){
      tv_ap[((i-1)*blockSize+1):(i*blockSize)] = tvApValue
    }
  }

  # reconstruct last thrown away part
  if (blockStateList[length(blockStateList)]==1){
    tv_ap[((i*blockSize)+1):length(tv_ap)] = tvApValue
  }

  return(tv_ap)
}



######################################################################################################################

#' @title tvStateExtractor
#' @author CH
#' @description extract tv state on total data from tv data and total data
#' @update 2015.12.03
#' @input totalData : c(timestamp, active_power, reactive_power)
#'        tvData : tv c(timestamp, active_ower, reactive_power)
#'        blockSize : classification time size
#' @output tv state vector
#'         e.g. 0 0 0 1 1 1 1 1 1 1 1 1 1 1 0...
#'
#'
#'
#'
#'
tvStateExtractor <- function(totalData, tvData, blockSize)
{
  total_ts <- totalData[,1]
  total_ap <- totalData[,2]
  total_rp <- totalData[,3]

  tv_ts <- tvData[,1]
  tv_ap <- tvData[,2]
  tv_rp <- tvData[,3]
  onThreshold <- 2 # if tv_ap >= (min(tv_ap))+onThreshold, tv_state = 1

  tv_state <- (tv_ap>=(min(tv_ap)+onThreshold))

  # tv on-off point
  tv_start <- rep(0, length(tv_ap))
  tv_end <- rep(0, length(tv_ap))

  # case handling : first time block on-off
  if (tv_state[1]==1){
    tv_start[1] <- 1
  }

  # (i-1)th state : off & ith state : on ->start
  for (i in 2:length(tv_state)){
    tv_start[i] <- ((tv_state[i]==1)&&(tv_state[i-1]==0))
  }


  # (i-1)th state : on & ith state : off -> end
  for (i in 1:(length(tv_state)-1)){
    tv_end[i] <- tv_state[i]==1&&tv_state[i+1]==0
  }

  # case handling : last time block on-off
  if (sum(tv_start)>sum(tv_end)){
    tv_end[length(tv_ap)] <- 1
  }

  tv_OnTime <- cbind(tv_ts[which(tv_start>0)], tv_ts[which(tv_end>0)])

  # label tv-state on total data
  tv_stateOnTotal <- rep(0, length(total_ap))

  if (nrow(tv_OnTime)>0){
    for (i in 1:length(tv_stateOnTotal)){
      for (j in 1:size(tv_OnTime,1)){
        if ((total_ts[i]>=tv_OnTime[j,1])&(total_ts[i]<=tv_OnTime[j,2]))
          tv_stateOnTotal[i] <- 1
      }
    }
  }

  # if 1/5 of block is on, block state is on.
  tv_blockState = rep(0, floor(length(tv_stateOnTotal)/blockSize))
  for (i in 1:length(tv_blockState)){
    blockSum = sum(tv_stateOnTotal[((i-1)*blockSize+1):(i*blockSize)])
    if (blockSum>(blockSize/5)){
      tv_blockState[i]<-1
    }
  }

  return(tv_blockState)
}
