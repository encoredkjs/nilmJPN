#####################################################################
####
#### Energy Dissagregation for washing machine
####
####      Soo-Heang Eo, PhD
####      Encored Techonologies, Inc
####
#####################################################################

#library(lubridate)
#library(quantreg)
#library(zoo)

######
# management JSON parameter
JSONforWashingMachine <- function(metaVer = 1, washingMachineVersion = "0.0.10", epsilon = 3,
    P.washingKernelSize = 180, Q.washingKernelSize = 180, V.washingKernelSize = 180, I.washingKernelSize = 180,
    P.washingEpsilon = 100, Q.washingEpsilon = 25, V.washingEpsilon = .5, I.washingEpsilon = .5,
    P.spindryKernelSize = 300, Q.spindryKernelSize = 300, V.spindryKernelSize = 300, I.spindryKernelSize = 300,
    P.spindryEpsilon = 1, Q.spindryEpsilon = 1, V.spindryEpsilon = 1, I.spindryEpsilon = 1,
    P.minUsage = 50, Q.minUsage = 5,
    P.localWindowSize = 10,
    spindry.startDelay = 1500, spindry.endDelay = 3000, spindry.duration = 600, spindry.difference = 2, spindry.max = 50,
    startTS = "2015-11-01 00:00:00", endTS ="2015-11-01 23:59:59", samplingRate = 15,
    computedTime = as.character(Sys.time()), debug = TRUE, training = FALSE, ...){

    resJSON = list()

    #####
    # META VERSION
    resJSON['meta-version'] = metaVer

    #####
    resJSON$shape_type = "WashingMachine"
    resJSON$washingMachineVersion = washingMachineVersion
    
    resJSON$ActivePower <- list(
    washingKernelSize = P.washingKernelSize,
    washingEpsilon = P.washingEpsilon,
    spindryKernelSize = P.spindryKernelSize,
    spindryEpsilon = P.spindryEpsilon,
    minUsage = P.minUsage,
    localWindowSize = P.localWindowSize)

    resJSON$ReactivePower <- list(
        washingKernelSize = Q.washingKernelSize,
        washingEpsilon = Q.washingEpsilon,
        spindryKernelSize = Q.spindryKernelSize,
        spindryEpsilon = Q.spindryEpsilon,
        minUsage = Q.minUsage)
  
    resJSON$Voltage <- list(
        washingKernelSize = V.washingKernelSize,
        washingEpsilon = V.washingEpsilon,
        spindryKernelSize = V.spindryKernelSize,
        spindryEpsilon = V.spindryEpsilon
        )

    resJSON$Current <-list(
        washingKernelSize = I.washingKernelSize,
        washingEpsilon = I.washingEpsilon,
        spindryKernelSize = I.spindryKernelSize,
        spindryEpsilon = I.spindryEpsilon
        )

    resJSON$epsilon = epsilon
    
    resJSON$SpinDry = list(
            startDelay = spindry.startDelay,
            endDelay = spindry.endDelay,
            duration = spindry.duration,
            difference = spindry.difference,
            max = spindry.max,
            fix = TRUE)

    resJSON$generation_info = list(
        data_used = list(
            start = startTS,
            end = endTS
            ),
        samplingRate = samplingRate,
        computed = computedTime,
        trainingMode = training,
        debug = debug)

    #class(resJSON) <- c(class(resJSON), "WashingMachine")
    return(resJSON)
}

######################
# Get washing mode
.zeroOneScale <- function(x){
    (x-min(x))/(max(x)-min(x))
}

.segmentWashingMachine <- function(ts, kernelSize){
    #if(!is.vector(ts)) stop("STOP in .segmentWashingMachine: ts must be a numeric vector.")
    g = length(ts) %/% kernelSize
    segmentedGroup = rep(1:g, each = kernelSize)
    remains = rep(g+1, times = (length(ts) - length(segmentedGroup)))
    segmentedGroup <- c(segmentedGroup, remains)
    #segmentedGroup = ((ts - ts[1]) %% kernelSize) == 0
    #segmentedGroup <- cumsum(segmentedGroup)
    return(segmentedGroup)
}

getWashingMode <- function(X, washerParas){
    cat("NILM1Hz now tries to find wahsing mode...")
    #########
    ## Step 1-1. Fit linear quantile regression with tau = 0.2 and 0.8 for Active Power
    timestamp = X$timestamp
    X$timestamp <- as.numeric(X$timestamp)
    Xap = as.zoo(X)
    #merge(y = as.zoo(X$active_power), x = as.zoo(X$timestamp), all = FALSE)

    xCoef = rollapply(Xap, 
        width = 30 * 15 ,
        FUN = function(z) {
            z <- as.data.frame(z)
            fit1 = quantile(z$active_power, 0.2)
            #coef(rq.fit(x = cbind(Int = 1, z$timestamp), y = z$active_power, tau = 0.15))
            fit2 = quantile(z$active_power, 0.8)
            fit3 = quantile(z$reactive_power, 0.2)
            fit4 = quantile(z$reactive_power, 0.8)
            fit5 = length(unique(round(z$reactive_power, -1)))
            fit6 = min(z$reactive_power)
            fit7 = max(z$reactive_power)
            #fit2 = coef(rq.fit(x = cbind(Int = 1, z$timestamp), y = z$active_power, tau = 0.85))
            #fit3 = coef(rq.fit(x = cbind(Int = 1, z$timestamp), y = z$reactive_power, tau = 0.15))
            #fit4 = coef(rq.fit(x = cbind(Int = 1, z$timestamp), y = z$reactive_power, tau = 0.85))
            pVal = tryCatch(
                dip.test(rescale(z$reactive_power))$p.value,
                error = function(e) 0
                )
            return(data.frame(
                apDiff = fit2 - fit1,
                rpDiff = fit4 - fit3,
                modes = fit5,
                rpMin = fit6,
                rpMax = fit7,
                pVal = pVal
                ))
        },
        by = 30 * 15 ,
        by.column = FALSE, 
        align = "right")

    ###########
    ## Detect washer based on active power
    b0APIndex = xCoef$apDiff > 100 & 
                xCoef$rpDiff > 10 & # 11
                xCoef$apDiff < 750 & # 1000
                xCoef$rpDiff < 150 & # 200
                xCoef$modes %in% c(2:20) & # 2:15
                xCoef$pVal <= 0.2
                
    # 18, 6
    xWhich = rollapply(b0APIndex, width = 20,
            FUN = function(z) sum(z) > 6,
        by = 10,
        by.column = FALSE, align = "right")

    if(any(xWhich)){
        xCandid = index(xWhich)[xWhich]
        checkCandid = apply(as.matrix(xCandid), 1, function(s, S) which(s == S), S = index(xWhich))
        checkIndex = index(xWhich)[checkCandid]
        checkIndex <- checkIndex - 3 * washerParas$ActivePower[['washingKernelSize']] * washerParas$generation_info$samplingRate
        checkIndex <- ifelse(checkIndex < index(b0APIndex)[1], index(b0APIndex)[1], checkIndex)
        if(length(checkIndex) > 2){
            checkIndex <- checkIndex[c(TRUE, diff(checkIndex) > 30 * 60 * washerParas$generation_info$samplingRate)]
        }
    } else{
        checkIndex  = NA
    }
    cat("Finished \n")

    if(washerParas$generation_info$trainingMode){
        estimatedTimeStamp = timestamp[checkIndex]
        resLearningMeta = data.frame(
            washingTimeStamp = estimatedTimeStamp,
            washingHour = hour(estimatedTimeStamp),
            washingDay = day(estimatedTimeStamp),
            washingMonth = month(estimatedTimeStamp),
            washingDayOfWeek = wday(estimatedTimeStamp)
            )

        # Calculate Box Duration Period (BDP) for every day of the week
        #bdp = rollapply(b0APIndex, width = 50,
        #    FUN = function(z) sum(z),
        #    by = 1,
        #    by.column = FALSE, align = "right")        
        #resLearningMeta$bdp <- bdp[which(index(b0APIndex) %in% checkIndex)]

        return(resLearningMeta)
    } else if(all(is.na(checkIndex))){
        return(list(
                startIndex = NA,
                startTimeStamp = NA
                ))
    } else{
        return(list(
            startIndex = checkIndex,
            startTimeStamp = timestamp[checkIndex]
            ))
    }
}
######################

getSpinDryMode <- function(X, resWashing, washerParas){
    ##############################################################################
    ### Get SpinDry Mode
    segmentStart = washerParas$SpinDry$startDelay * washerParas$generation_info$samplingRate
    segmentEnd = washerParas$SpinDry$endDelay * washerParas$generation_info$samplingRate
    endTimeStamp = endIndex = rep(0, length(resWashing$startTimeStamp))

    endIndex <- resWashing$startIndex + segmentEnd
    endTimeStamp <- X$timestamp[endIndex]
    #for(i in 1:length(checkTimeStamp)){
        #candidateIndex = (checkIndex[i] ):(checkIndex[i] + segmentStart + segmentEnd - 1)

        #Pnew = .localCentering(x = X$active_power[candidateIndex], 
        #    localWindow = jsonParas$ActivePower$localWindowSize, 
        #    type = "mean")$Xnew

        #PsegWashing = .segmentWashingMachine(
        #    ts = X$timestamp[candidateIndex], 
        #    kernelSize = jsonParas$ActivePower$spindryKernelSize)

        #rowValue = max(PsegWashing) - 1
        #checkMat = checkTestStat = matrix(0, nrow = rowValue, ncol = 2)

        #for(j in 1:rowValue){
        #    PtmpTable = .buildContingencyTable(C1 = Pnew[PsegWashing == j ], 
        #                                       C2 = Pnew[PsegWashing == (j+1)])
        #    checkMat[i,1] <- .checkPeakNumber(PtmpTable, 10)
        #    if(checkMat[i,1]){
        #        checkTestStat[i,1] <- max(chisq.test(PtmpTable)$statistic, 0.01)
        #    }
        #}

        #if( any(candidateSet > 10000) ){
            #endIndex[i] <- index(candidateSet[which.max(candidateSet)])
            #endIndex[i] <- checkIndex[i] + segmentEnd
            #endTimeStamp[i] <- X$timestamp[endIndex[i]]
        #} else {
            #endTimeStamp[i] <- NA
            #endIndex[i] <- NA
        #}                
    #}

    finalIndexCheck = cbind(resWashing$startIndex, endIndex)
    finalIndexCheck <- apply(finalIndexCheck, 1, function(x) all(!is.na(x)) )

    return(list(startIndex   = resWashing$startIndex[finalIndexCheck],
                startTimeStamp = resWashing$startTimeStamp[finalIndexCheck],
                endIndex     = endIndex[finalIndexCheck],
                endTimeStamp = as.POSIXct(endTimeStamp[finalIndexCheck], origin = "1970-01-01 00:00:00")
                )
    )
}

######################
# Result manipulation
.checkPeakNumber <- function(conTable, numPeaks){
    booleanIndex = apply(conTable, 1, function(x, numPeaks) all(x > numPeaks), numPeaks = numPeaks)
    booleanIndexSum = sum(booleanIndex)
    if(booleanIndexSum == 2L){
        return(TRUE)
    #} else if (booleanIndexSum == 1L){
    #    return("hold") # maybe TV or Vacuum Cleaner?
    } else{
        return(FALSE)
    }
}

.buildContingencyTable <- function(C1, C2){

    if( sum(C1 == 0) > 500 | sum(C2 == 0) > 500) {
        return(matrix(0, nrow = 2, ncol = 2))
    } else{
        C1 = ifelse(C1 > 0, 1, -1)
        C2 = ifelse(C2 > 0, 1, -1)
        conTable <- matrix(c(table(C1), table(C2)), nrow = 2, ncol = 2)
        return(conTable)
    }
}

.localCentering <- function(x, localWindow, type = c("mean", "var", "LB")) {
    if(!is.vector(x)) stop("STOP in .localCentering: x must be a numeric vector.")

    type <- match.arg(type)

    Xobs = length(x)
    Xdim = c(localWindow, ceiling( Xobs / localWindow))
    if(Xobs != prod(Xdim)){
        numNA = prod(Xdim) - Xobs
        x <- c(x, rep(NA, numNA))
    }
    dim(x) = Xdim
    if(type == "mean"){
        Xcenter = colMeans(x, na.rm = TRUE)
        res = x - rep(Xcenter, rep.int(nrow(x), ncol(x)))
        resList = data.frame(Xnew = as.vector(res), Xcenter = rep(Xcenter,localWindow) )
        if(Xobs != prod(Xdim)) resList <- resList[1:(nrow(resList)-numNA),,drop = FALSE]
        #class(resList) <- c("CommonModel", "WashingMachine")
    } else if (type == "var"){
        Xvar = apply(x, 2, var, na.rm = TRUE)
        resList = data.frame(Xvar = Xvar)
    } else if (type == "LB"){
        # find lower bound of the local window
        Xcenter = apply(x, 2, function(z) quantile(z, 0.01, na.rm = TRUE))
        res = x - rep(Xcenter, rep.int(nrow(x), ncol(x)))
        resList = data.frame(Xnew = as.vector(res), Xcenter = rep(Xcenter,localWindow) )
    }
    return(resList) 
}

calculatePowerUsageWashingMachine <- function(startIndex, endIndex, resObject, P, Q, washerParas){
    ###############################
    ### Under Construction

    for(i in 1:length(startIndex)){
        apHat = .localCentering(P[startIndex[i]:endIndex[i]], localWindow = 10 *  washerParas$generation_info$samplingRate, type = "LB")
        if(length(P[startIndex[i]:endIndex[i]]) !=  length(apHat$Xnew)){
            Phat = na.omit(apHat$Xnew)[1:length(P[startIndex[i]:endIndex[i]])]
        }
        resObject$p[startIndex[i]:endIndex[i]] <- pmax(Phat, 0 )
        resObject$q[startIndex[i]:endIndex[i]] <- Q[startIndex[i]:endIndex[i]]
    }
    if(washerParas$generation_info$debug){
        plot(resObject$p[startIndex:endIndex], type = "l", ylab = "Estimated Active Power")
        plot(resObject$q[startIndex:endIndex], type = "l", ylab = "Estimated Reactive Power")
    } 
    return(resObject)
}



######################
# Training Mode (Fit learning model)
generate.meta.washer <- function(X = NULL, training = TRUE, debug = FALSE, ...){

    ###
    # Step 0. Check integrity 
    if(!is.data.frame(X)) {
        X <- as.data.frame(X)
    }
    #X <- tbl_df(X)

    # Get tuning parameters
    washerParas = JSONforWashingMachine(...)
    washerParas$generation_info$debug <- debug
    washerParas$generation_info$trainingMode <- training

    ###
    # Step 1. find washing mode for every day of the week
    #
    # A. divide dataset into the day of the week
    # B. calculate box duration period (BDP) for every day of the week
    # C. 
    #resLearning = getWashingMode(X, washerParas)


    washerParas$generation_info$trainingMode <- FALSE
    return(washerParas)
} 


######################
# predict forecasting time and consumed load for washing machine
predict.washer <- function(meta, data, debug = FALSE){
    cat("Search washing machine appliances based on the unimodality test...  \n")
    ###
    # set meta
    washerParas = meta
    washerParas$generation_info$debug <- debug
    washerParas$generation_info$samplingRate <- 15
        
    ###
    # storage for the result object
    resObject = data.frame(
        timestamp = data$timestamp,
        p = rep(0, length = nrow(data)),
        q = rep(0, length = nrow(data))
    )
    
    ###################
    # Detect washing mode
    resWashing = getWashingMode(X = data, washerParas)
    
    if(is.na(resWashing$startIndex)){
        return(resObject)
    }


    ###################
    # Detect spindry mode
    resSpinDry = getSpinDryMode(X = data, resWashing, washerParas)

    ###############################################################################
    # NOTICE ME: need to improve by vectorized calculcation ussing apply()-related functiosn
    ###############################################################################
    if(length(resSpinDry$endIndex) >= 1){
        resObject = calculatePowerUsageWashingMachine(
                startIndex = resSpinDry$startIndex, 
                endIndex = resSpinDry$endIndex, 
                resObject = resObject, 
                P = data$active_power, 
                Q = data$reactive_power,
                washerParas = washerParas)
        return(resObject)
    } else{
        return(resObject)
    }

}
# END
################################

