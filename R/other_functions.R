## Potential evapotranspiration  ###############################################
# Modified Blaney-Criddle method

potET <- function(Temp){
  #potET calculates mean daily potential evapotranspiration
  #Temp   time series object of air temperatures of class 'xts'/'zoo' or coercable to xts/zoo.
  #       Resolution should be at least daily.
  if (!is.xts(Temp) || !is.zoo(Temp)) {
    try(Temp <- xts(Temp[,2], order.by = Temp[,1]))
  }
  # Solar indices (Virtual Water-Science Laboratory)
  # Potential daily sunshine duration from (Meszaroš et al. 2002, Parajka et al. 2003)
  SI <- c(0.1966, 0.2341, 0.2715, 0.3090, 0.3464, 0.3652,
          0.3558, 0.3277, 0.2856, 0.2481, 0.2060, 0.1873)

  Temp <- to_timestep(Temp,by = "days")
  Temp$month <- month(Temp)
  PET <- -1.55 + 0.96*(8.128 + 0.457 * Temp[,1])*SI[Temp[,2]]
  names(PET) <- "PET(mm/day)"
  return(PET)
}



## rbind at specific datetime ##################################################
rbind_at <- function(x,y, at){
  #binds 'x' and 'y' together only keeping only data from 'x' before 'at' and
  #only data from 'y' after 'at'

  #'x' a xts object. Single or multivariate
  #'y' second xts object. Single or multivariate (same as 'x')
  #'at' datetime value or vector when binding should a occur.

  #check number of variables in 'x' and 'y'
  if (ncol(x)!=ncol(y)){
    stop("'x' and 'y' must have same number of columns!")
  }

  #check 'at'
  if (!(length(at) == 1 || length(at) == ncol(x))){
    stop("'at' should be of length 1 or ncol(x)!")
  }

  #make sure 'at' is a POSIXct object
  at = as.POSIXct(at, tz="Etc/GMT-1")
  if (length(at) == 1){
    at = rep(at,ncol(x))
  }
  #create empty xts object
  by = as.numeric(median(diff(index(x))), "secs")
  xy <- xts(order.by = seq(start(x),end(y), by = by))

  for (i in 1:ncol(x)){
    temp <- rbind(x[paste0("/",at[i]),i],y[paste0(at[i]+by,"/"),i])
    xy <- merge(xy,temp)
  }
  names(xy)<-names(x)
  return(xy)
}

## make index unique ###########################################################
make.ind.unique <- function(x, FUN = "mean"){
  #'FUN' function to use on dupplicate index enties (mean, first,last)
  if(FUN == "mean"){
    ind_dupl <- index(x)[which(duplicated(index(x)))]
    if(length(ind_dupl)!=0){
      x_mean <- xts(aggregate(x[ind_dupl], index(x[ind_dupl]), mean,na.rm=T))
      x <- make.index.unique(x, drop = T, fromLast = F)
      x[ind_dupl,] <- x_mean[,]
    }
    return(x)
  }
  if(FUN == "first"){
    return(make.index.unique(x, drop = T, fromLast = F))
  }
  if(FUN == "last"){
    return(make.index.unique(x, drop = T, fromLast = T))
  }
}

## remove outliers #############################################################
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

## Mean monthly values  ########################################################
MonMean <- function(data, plot = F){
  #MonMean calculates average value for each month (e.g. mean of all Aprils in dataset)
  #INPUT:
  #   'data'  xts object. could be uni- or multi-variate
  #VALUE;
  #   data.frame of size 12 x ncol(data)

  # change timestep to monthly if not already
  data <- to_timestep(data, "month")

  #array of logical values of size nrow(data)x12
  #each column indicates which rows in data corespond to certain month
  y <- mapply(function(x){.indexmon(data) == x}, 0:11)

  # calculates mean for each month for each variable
  monmean <- apply(data,2, function(x){apply(y,2,function(a){mean(x[a],na.rm=T)})})
  months <- c("Jan", "Feb", "Mar", "Apr","May", "Jun",
              "Jul", "Aug", "Sep","Oct","Nov","Dec")
  rownames(monmean)<- months
  if(plot){
    x11()
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    matplot(monmean,type = c("b"),pch=1, main="Monthly average",
            xlab = "Month", ylab = "Monthly average", col = 1:ncol(monmean),
            lwd = 1.5, lty = 1:ncol(monmean))
    legend("right", colnames(monmean), inset=c(-0.25,0), pch = 1,
           lty = 1:ncol(monmean), col = 1:ncol(monmean),bty='n',
           xpd = T, title = "Stations")

    # normalized graph
    x11()
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    nmonmean <- apply(monmean, 2, function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x, na.rm=T))})
    matplot(nmonmean, type = c("b"),pch=1, col = 1:ncol(monmean),
            xlab = "Month", ylab = "Normalized Monthly average",
            lwd = 1.5, lty = 1:ncol(monmean), main="Normalized Monthly average")
    legend("right", colnames(monmean), inset=c(-0.25,0), pch = 1,
           lty = 1:ncol(monmean), col = 1:ncol(monmean),bty='n',
           xpd = T, title = "Stations")
  }#END IF

  return(monmean)
}

#==============================================================================#
# Mode ####
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#==============================================================================#
# Data normalization ####

## normalize data to the interval [0,1]
dfNorm <- function(x){
  as.data.frame(lapply(x, function(y) {(y - min(y, na.rm=TRUE))/(max(y,na.rm=TRUE) -
                                                                   min(y, na.rm=TRUE))}))
}
# 'stdx' function from hydroTSM
stdxts <- function(x){
  return(xts(apply(x, 2, stdx), order.by = index(x)))
}

## subtract min or mean value of the time-series

XtsChange <- function(x, FUN = min){
  xts(vapply(x, function(col) col - FUN(col, na.rm = T),
             FUN.VALUE = numeric(nrow(x))),
      order.by = index(x))
}

#==============================================================================#
# Date of min, max ####
date.min <- function(x){
  #date.min finds the date of the minimum in the provided 'xts' object 'x'

  #INPUTS:
  # 'x' xts object (could be multivariate)

  #Value:
  # vector of minima dates with the length of ncol(x)

  d <- c()
  for(i in 1:ncol(x)){
    x1 <- x[!is.na(x[,i]),i]
    d <- append(d, index(x1[which.min(x1)]))
  }

  return(d)
}

date.max <- function(x){
  #date.max finds the date of the minimum in the provided 'xts' object 'x'

  #INPUTS:
  # 'x' xts object (could be multivariate)

  #Value:
  # vector of minima dates with the length of ncol(x)

  d <- c()
  for(i in 1:ncol(x)){
    x1 <- x[!is.na(x[,i]),i]
    d <- append(d, index(x1[which.max(x1)]))
  }

  return(d)
}

#==============================================================================#
# rollmean with NA ####
rollmean2 <- function(x,k=3, sides=2, method="convolution", shrink = F){
  #rollmean2 calculates rolling mean of timeseries with NAs
  if(sides == 2 & k%%2 == 0) warning("if sides = 2, 'k' should be odd!")
  if (shrink == F){
    return(as.xts(apply(x, 2, function(z) stats::filter(z, rep(1/k, k), sides=sides, method="convolution")),
                  order.by = index(x)))
  }
  if(shrink == T){
    return(as.xts(apply(x, 2, function(z) stats::filter(z, rep(1/k, k), sides=sides, method="convolution")),
                  order.by = index(x))[ceiling((k-1)/2+1):(nrow(x)-floor((k-1)/2)),])
  }
}

#==============================================================================#
#  Precipitation Events ####

PrecEvent <- function(prec, t, sigprec, plot = F){
  #PrecEvent finds the episodes of significant precipitation in the time series

  #INPUT:
  #'prec' - precipitation time series in xts format
  #'t' - duration of precipitation effect. numeric value in the units of 'prec'
  #'sigprec' - numeric value of minimum precipitation to consider
  #'plot' - should a plot be drawn

  #OUTPUT:
  #data frame with start time, end time and string of the period ("2017-01-03/2017-01-06")

  # if 'sigprec' is missing set it based on the mean precipitation value
  if(missing(sigprec)){
    sigprec = mean(prec[prec > 0],na.rm=T)*0.8
  }

  # Time step in seconds
  dt <- as.numeric(median(diff(index(prec))), "secs")

  # Significant precipitation times
  sigPrec <- index(prec[which(prec > sigprec)])
  # Significant precipitation start and end times
  sigPrec <- data.frame(Start = sigPrec, End = sigPrec + t*dt)

  ii <- 1
  jj <- 1
  #Data frame with stored episode times
  Events <- data.frame(Start = sigPrec[1,1], End = sigPrec[1,2])
  while (ii < nrow(sigPrec)){
    if(sigPrec[ii,2] >= sigPrec[ii+1,1]){
      ii <- ii + 1
    }else{
      Events[jj,2] <- sigPrec[ii,2]
      Events[jj+1,1] <- sigPrec[ii+1,1]
      ii <- ii + 1
      jj <- jj + 1
    }
  }
  Events[jj,2] <- sigPrec[ii,2]
  #Episode period as string
  Events$Period <- paste0(Events[,1], "/", Events[,2])
  cat("Number of precipitation epsiodes found: ",nrow(Events),"\n")
  #Total rainfall
  Events$P_total <- sapply(Events$Period, function(z) sum(prec[z],na.rm = T))
  #Max intensity
  Events$i_max.mm.h <- sapply(Events$Period, function(z) max(prec[z],na.rm = T)/dt*3600)
  #Plot
  if(plot){
    PrecCum <- cumsum(prec)
    dy <-dygraph(merge(PrecCum,PrecCum[Events[,3]]))%>%
      dyOptions(useDataTimezone = T, fillGraph = T)
    print(dy)
  }

  return(Events)
}

#------------------------------------------------------------------------------#

PrecipitationEvent <- function(prec, duration.h, sigprec, plot = F){
  #PrecEvent finds the episodes of significant precipitation in the time series

  #INPUT:
  #'prec' - precipitation time series in xts format
  #'duration.h' - duration of precipitation effect. numeric value in the hours
  #'sigprec' - numeric value of minimum precipitation to consider
  #'plot' - should a plot be drawn

  #OUTPUT:
  #data frame with start time, end time and string of the period ("2017-01-03/2017-01-06")

  # if 'sigprec' is missing set it based on the mean precipitation value
  if(missing(sigprec)){
    sigprec = mean(prec[prec > 0],na.rm=T)*0.8
  }

  # Time step in seconds
  dt <- as.numeric(median(diff(index(prec))), "secs")

  # Event duration in seconds
  t <- duration.h * 3600

  # minimum total precipitation
  min.P.tot <- 4 #[mm]
  P.tot <- rollapply(prec, width = round(t/dt), FUN = sum, na.rm = T, align = "right")

  #
  P.logic <- prec
  P.logic[P.logic < sigprec] <- 0
  P.rle <- rle(as.double(P.logic))
  P.dry <- rep(P.rle$lengths >= t/dt, times = P.rle$lengths)
  P.dry.rle <- rle(P.dry)
  for (i in 1:length(P.dry.rle$lengths)){
    if(P.dry.rle$values[i] == F){
      ind.end <- sum(P.dry.rle$lengths[1:i])
      ind.start <- ind.end - P.dry.rle$lengths[i]
      P.sum <- sum(prec[ind.start:ind.end])
      if(P.sum < min.P.tot) P.dry.rle$values[i] = T
    }
  }
  P.dry <- inverse.rle(P.dry.rle)

  Events <- data.frame(Start = index(prec)[cumsum(P.dry.rle$lengths) - P.dry.rle$lengths + 1],
                       End = index(prec)[cumsum(P.dry.rle$lengths) + t/dt])[!P.dry.rle$values,]

  #Episode period as string
  Events$Period <- paste0(Events[,1], "/", Events[,2])
  cat("Number of precipitation epsiodes found: ",nrow(Events),"\n")
  #Total rainfall
  Events$P_total <- sapply(Events$Period, function(z) sum(prec[z],na.rm = T))
  #Max intensity
  Events$i_max.mm.h <- sapply(Events$Period, function(z) max(prec[z],na.rm = T)/dt*3600)
  #Plot
  if(plot){
    PrecCum <- cumsum(prec)
    dy <-dygraph(merge(PrecCum,PrecCum[Events[,3]]))%>%
      dyOptions(useDataTimezone = T, fillGraph = T)
    print(dy)
  }

  return(Events)
}

#------------------------------------------------------------------------------#



#==============================================================================#
#  Recharge Events ####

RechargeEvents <- function(GWL,prec,minD,PrecEvents, plot=F){

  #changes in GWL
  dGWL <- diff.xts(GWL)
  colnames(dGWL) <- paste0("d",colnames(GWL))
  # significant Recharge/change in GWL
  sigR <- index(dGWL[dGWL > minD])

  #time step
  dt <- as.numeric(median(diff(index(GWL))), "secs")

  ## find episodes of significant recharge
  REvents <- data.frame(Start=sigR[1]-dt, End = sigR[1])
  ii <- 1
  jj <- 1
  kk <- 0
  while (ii < length(sigR)){
    #check if next sig. recharge is within two time steps
    if(as.numeric(sigR[ii+1]-sigR[ii],units = "secs") <= 2*dt){
      ii <- ii + 1
      #check if recharge in next step is positive
    }else if(dGWL[sigR[ii]+(kk+1)*dt] > 0){
      kk <- kk + 1
      #save the end of recharge event and start of the next one
    }else{
      REvents[jj,2] <- sigR[ii] + kk*dt
      if(REvents[jj,2] >= sigR[ii+1]-dt){
        #find the first sig. recharge after the end of last episode
        ii <- which((sigR - dt) > REvents[jj,2])[1]
        REvents[jj+1,1] <- sigR[ii]-dt
      }else{
        REvents[jj+1,1] <- sigR[ii+1]-dt
        ii <- ii + 1
      }
      jj <- jj + 1
      kk <- 0
    }
  }
  # save the end of last recharge event
  REvents[jj,2] <- sigR[ii]

  #recharge period as string
  REvents$Period <- paste0(REvents[,1], "/", REvents[,2])
  cat("Number of recharge epsiodes found: ",nrow(REvents),"\n")

  # only when precipitation data is missing
  if(!(missing(PrecEvents) & missing(prec))){
    ## Precipitation events
    if(missing(PrecEvents)){
      dtp <- as.numeric(median(diff(index(prec))), "secs")
      PrecEvents <- PrecEvent(prec, t=2*24*3600/dtp)
    }
    PrecEvents$Interval <- interval(PrecEvents[,1],PrecEvents[,2])

    ## Relate Recharge and Precipitation events
    temp <- REvents[0,]
    jj <- 1
    for(ii in 1:nrow(REvents)){
      if(any(REvents[ii,1] %within% PrecEvents$Interval)){
        temp[jj,]<-REvents[ii,]
        jj <- jj+1
      }
    }
    cat("Number of recharge episodes not attributable to precipitation:",
        nrow(REvents)-nrow(temp),"\n")
    print(REvents[!REvents[,1]%in%temp[,1],])
    REvents <- temp
  } #end if missing

  #Plot
  if(plot){
    dy <-dygraph(merge(GWL,dGWL,GWL[REvents$Period]))%>%
      dyOptions(useDataTimezone = T)%>%
      dySeries(paste0(colnames(GWL),".1"), axis="y2",fillGraph = T)%>%
      dySeries(colnames(GWL),axis="y2")%>%
      dyLimit(0, color="red")%>%
      dyLimit(minD, color="black")
    print(dy)
    if(!missing(prec)){
      PrecCum <- cumsum(prec)
      dy2 <-dygraph(merge(PrecCum,PrecCum[PrecEvents[,3]],GWL,GWL[REvents$Period]))%>%
        dyOptions(useDataTimezone = T, fillGraph = T)%>%
        dySeries(colnames(PrecCum), axis="y2")%>%
        dySeries(paste0(colnames(PrecCum),".1"), axis="y2")
      print(dy2)
    }

  }
  return(REvents)
}

#==============================================================================#
#  Start of rise of hydrograph ####
# Last Update: 14.012019
# Notes: Works bad - not to be used

HydRise <- function(x,sig.rise = 0.001, k = 20){
  #HydRise determines the sart time of the rise in water level response to a
  #rainfall event by finding the first point of three consecutive increases in
  #water level of at least 'sig.rise'

  #INPUTS:
  #'x' xts of water level data
  #'sig.rise' significant rise of water level in units of 'x'

  #Value:
  #Vector of start times

  if(!is.xts(x))stop("'x' has to be an xts objects")
  x <- rollmean2(x, k = k)

  #Compute the lengths and values of runs of equal values in a vector
  runs <- rle(as.vector(diff(x) >= sig.rise))
  # Indices of the runs with length of at least 3
  runs.3 <- which(runs$values == TRUE & runs$lengths >= 3)

  if(!any(runs.3)) stop("No events found")

  runs.lengths.cumsum <- cumsum(runs$lengths)
  newindex = ifelse(runs.3 > 1, runs.3 - 1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)

  return(index(x)[starts])
}

#==============================================================================#
#  RISE Method ####
RISE <- function(GWL,RechargeEvents, Sy=0.036, pos=T, prec){

  dGWL <- diff.xts(GWL)
  if(!missing(RechargeEvents)){
    dGWL <- dGWL[RechargeEvents[,3]]
  }
  dGWLpos <- dGWL[dGWL>0]
  R <- sum(Sy*dGWL, na.rm=T)*1000
  Rpos <- sum(Sy*dGWLpos, na.rm=T)*1000
  cat("Recharge:",R,"\n")
  cat("Recharge (only positive):",Rpos,"\n")
  if(!missing(prec)){
    P <- sum(prec)
    RPR <- R/P
    RPRpos <- Rpos/P
    cat("Precipitation:", P, "\n")
    cat("Recharge-Precipitation Ratio:",RPR,"\n")
    cat("Recharge-Precipitation Ratio (only positive):",RPRpos,"\n")
  }

  if(pos){
    return(Rpos)
  }else{
    return(R)
  }
}

#==============================================================================#
####### Drought periods
#Needs improvements - better use 'dryPer'
drought_periods <- function(prec, min_dur = 10){
  #function that finds longest periods of drought in precipitation time series
  start = index(prec)[which(prec>0)[which.max(diff(which(prec>0)))]]
  end = index(prec)[which(prec>0)[which.max(diff(which(prec>0)))+1]]

  if(end-start < 10) stop(paste("No periods of droughts of minimum duration 'min_dur'=",min_dur ,"were found."))

  return (c(start, end))
}

#==============================================================================#
#  Dry period ####
dryPer <- function(prec,t, mindrytime=1, sigprec,plot=T){
  #dryPer finds the dry periods in the precipitation time-series

  #INPUT:
  #'prec' - precipitation time series in xts format
  #'t' - duration of precipitation effect. numeric value in time-step of 'prec'
  #'mindrytime' - minimum length of dry period in time-step of 'prec'
  #'sigprec' - numeric value of minimum precipitation to consider
  #'plot' - should a plot be drawn

  #OUTPUT:
  #data frame with start time, end time and string of the dry period ("2017-01-03/2017-01-06")

  # if 'sigprec' is missing set it based on the mean precipitation value
  if(missing(sigprec)){
    sigprec = mean(prec[prec > 0],na.rm=T)*0.8
    cat("'sigprec' set to",sigprec,"\n")
  }

  #time step
  dt <- as.numeric(median(diff(index(prec))), "secs")

  # precipitation events
  rain <- PrecEvent(prec,t=t,sigprec=sigprec, plot=F)

  # dry periods
  dry <- data.frame(Start=rain[1,1],End=rain[1,2])
  ii <- 1 #'rain' indices
  jj <- 1 #'dry' indices
  # first period
  if (rain[1,1] - dt != index(prec[1,])){
    dry[1,1] <- index(prec[1,])
    dry[1,2] <- rain[1,1] - dt
    ii <- ii + 1
    jj <- jj + 1
  }#if

  while(ii < nrow(rain)){
    dry[jj,1] <- rain[ii,2]
    dry[jj,2] <- rain[ii+1,1] - dt
    ii <- ii + 1
    if(abs(as.numeric(dry[jj,2]-dry[jj,1],"secs")) > mindrytime * dt){
      jj <- jj + 1
    }#if
  }#while

  #last period
  if(rain[nrow(rain),2] != tail(index(prec),1)){
    dry[jj,1] <- rain[ii,2]
    dry[jj,2] <- tail(index(prec),1)
  }

  #Episode period as string
  dry$Period <- paste0(dry[,1], "/", dry[,2])
  cat("Number of dry epsiodes found: ",nrow(dry),"\n")

  #Episode duration
  dry$Duration <- as.numeric(dry$End-dry$Start,units="days")

  #Plot
  if(plot){
    PrecCum <- cumsum(prec)
    data <- merge(PrecCum,PrecCum[dry[,3]],PrecCum[rain[,3]])
    colnames(data) <- c("CumPrecipitation","DryPeriod","PrecEvent")
    dy <-dygraph(data)%>%
      dyOptions(useDataTimezone = T,colors=c("black","red","blue"))%>%
      dySeries("DryPeriod", fillGraph = T)%>%
      dySeries("PrecEvent", fillGraph = T)
    print(dy)
  }

  return(dry)
}#function



#==============================================================================#
## HOAL zones ####
HOALzone <- PiezoInstall()[,c("Station","Group")]
HOALzone$Group <- as.factor(HOALzone$Group)

#==============================================================================#
## Hysteresis paper ####

# Lag times
lag_times <- function(data,max.lag=1){
  #INPUT:
  #'data' multivariate xts object
  #'max.lag' max lag time in days to pass to ccfxts function
  #
  #Value:
  #data frame with one row for each station pair (long format)

  data <- rollmean2(data,k=12)
  data_diff <- diff.xts(data) # differentiented timeseries
  lags <- data.frame(st1=NA,st2=NA, lag.h=NA, lag.h.peak=NA, cor=NA) #data.frame to store results
  k=1 # row index
  for (i in 1:ncol(data_diff)) {
    for (j in 1:ncol(data_diff)) {
      temp <- ccfxts(data_diff[,i],data_diff[,j], lag = max.lag)
      lags[k,"st1"] <- colnames(data)[i]
      lags[k,"st2"] <- colnames(data)[j]
      lags[k,"lag.h"] <- ifelse(temp$maxcor< 0.4,NA,temp$lagmax)
      lags[k,"lag.h.peak"] <- head(which(stdxts(data[,i])==1)-which(stdxts(data[,j])==1),1)/12
      lags[k,"cor"] <- temp$maxcor
      k=k+1
    } #end for
  } #end for
  return(lags)
} #end function

## Hysteresis index
#Last update: 9.1.2019

# Literature:
#  Lloyd et al. (2016) Testing an improved index for analysing storm
#   discharge-concentration hysteresis

HI <- function(Q, C, sections = 5){
  #HI function calculates Lloyd Hysteresis index

  #INPUT:
  #'Q' discharge data for an event, data.frame or xts
  #'C' concentration or other dependant data for the same event as 'Q'
  #'sections' numeric value of data section size in percent

  #VALUe:
  # A data.frame with HI values for each section

  # NA check
  if(any(is.na(Q),is.na(C))) {
    #warning("Input data contains NAs - NA returned")
    Q.quant <- seq(from = sections, to = 100 - sections, by = sections)/100
    return(data.frame(Q.quant, HI = NA))
  }

  if(is.data.frame(Q) & is.data.frame(C)) {
    QC <- cbind(Q, data.frame(approx(C, xout = Q[,1]))[,2])

  }else if(is.xts(Q) & is.xts(C)){
    QC <- fortify.zoo(na.approx(merge(Q,C), xout = index(Q)))

  }else {
    stop("'Q' and 'C' must both be either data.frame or xts")
  }

  # normalize
  QC.norm <- data.frame(Index = QC[,1], dfNorm(QC[,2:3]))

  # Rising and Falling limb
  ind.Qmax <- base::which.max(QC[,2])
  QC.RL <- QC.norm[1:ind.Qmax,]
  QC.FL <- QC.norm[(ind.Qmax + 1):nrow(QC.norm),]

  # Q quantiles
  Q.quant <- seq(from = sections, to = 100 - sections, by = sections)/100

  #browser()
  if(nrow(QC.RL) < length(Q.quant) | nrow(QC.FL) < length(Q.quant)){
    return(data.frame(Q.quant, HI = NA))
  }

  C.RL <- data.frame(approx(QC.RL[,2:3], xout = Q.quant, ties = mean))
  C.FL <- data.frame(approx(QC.FL[,2:3], xout = Q.quant, ties = mean))

  # Hysteresis index
  HI <- data.frame(Q = Q.quant, HI = C.RL[,2] - C.FL[,2])

  return(HI)
}

#------------------------------------------------------------------------------#
## Hysteresis index Zuecco
# Reference: Zuecco G., Penna, D., Borga M., van Meerveld H.J. (2016) A versatile
#               index to characterize hysteresis between hydrological variables
#               at the runoff event timescale

#' Area of the Hysteresis loop
#'
#' \code{HystArea2} returns the area of the Hysteresis loop between -1 and 1.
#'
#' @details
#' shape.check was added to discarde or indicate when data is not appropriate to
#' calculate the hysteresis area. This can happen when C only increases in the time of the event
#'
#' @param Q vector, data.frame or xts of event discharge data (and second variable data)
#' @param C vector, data.frame or xts of second event response variable, optional if
#' \code{Q} is appropriate object.
#' @param min.Q Numeric. Minimum normalized \code{Q} from where the area will be
#'   calculated.
#' @param shape.check Logical. Should the the data be checked for loop shape?
#' @param min.dC Numeric. Minimum change in \code{C} to calculate area.
#'   Only when \code{shape.check == T}
#' @param min.narrow Minimum ration of last point from 0 or 1 to still
#'   calculate area
#'
#' @return
#' Numeric value of area between -1 and 1.
#' Returns Inf if data does not form a loop (only if \code{shape.check == T}).
#' Returns NA if change of C is less or equal to min.dC
#' (only if \code{shape.check == T})
#'
#' @examples\dontrun{
#' HystArea2(Q,C)
#' HystArea2(Q, C, min.Q = 0.1)
#' HystArea2(Q, C, min.Q = 0.1, shape.check = T)
#' }
#'
#' @export
#'
HystArea2 <- function(Q, C = NULL, min.Q = 0.1,
                      shape.check = F, min.dC = 0.02, min.narrow = 0.05){
  if(missing(C)){
    if(length(Q) == 0) stop("'C' is not provided and 'Q' is of length 0")
    if(is.xts(Q)){
      if(ncol(Q)==2){
        QC <- coredata(Q)
      }else{
        stop("If no 'C' is provided, 'Q' must have 2 columns!")
      }
    }else if(is.data.frame(Q)){
      if(ncol(Q)==2){
        QC <- Q
      }else if(ncol(Q) == 3 & class(Q[,1])){
        QC <- Q[,2:3]
      }else{
        stop("'C' is missing and 'Q' is not of correct format!")
      }
    }else{
      stop("'C' is missing and 'Q' is not of correct format!")
    }
    # C is provided
  }else{
    if(length(Q) == 0 | length(C) == 0) stop("Provided 'C' or 'Q' is of length 0")
    if(is.vector(Q)) Q <- as.data.frame(Q)
    if(is.vector(C)) C <- as.data.frame(C)

    if(is.xts(Q) & is.xts(C)){
      QC <- coredata(merge(Q,C))
    }
    if(is.data.frame(Q) & is.data.frame(C)){
      if(ncol(Q) == 1 & ncol(C) == 1){
        QC <- cbind(Q,C)
      }else if(ncol(Q) == 2 & ncol(C) == 2){
        if(is.POSIXct(Q[,1]) & is.POSIXct(C[,1])){
          colnames(Q) <- c("Index","Q")
          colnames(C) <- c("Index","C")
          QC <- merge(Q,C)
        }else{
          stop("'Q' is not of correct format!")
        }
      }else{
        stop("'Q' and 'C' are not of recognisable foramt.")
      }
    }else{
      stop("'Q' and 'C' are not of recognisable foramt.")
    }
  }

  QC <- QC[complete.cases(QC),]
  if(nrow(QC) < 4) {
    warning("Provided data is too short or does not overlap in time! NA returned.")
    return(NA)
  }
  QC.norm <- apply(QC,2,scales::rescale)

  # if(is.vector(Q)) Q <- as.data.frame(Q)
  # if(is.vector(C)) C <- as.data.frame(C)
  #
  # if(missing(C)){
  #   if(is.xts(Q) & ncol(Q) > 1){
  #     C <- as.data.frame(Q[,2])
  #     Q <- as.data.frame(Q[,1])
  #   }else if (is.data.frame(Q) & ncol(Q) > 2){
  #     if(ncol(Q)>2)
  #     C <- as.data.frame(Q[,3])
  #     Q <- as.data.frame(Q[,2])
  #   }else{
  #     stop("'C' missing with no default")
  #   }
  # }
  #
  # if(nrow(Q) < 9 | nrow(C) < 9) return(NA)
  #
  # # create data.frame
  # if(is.data.frame(Q) & is.data.frame(C)) {
  #   if(length(which(is.na(Q[,2]))) > 0.2*nrow(Q)) return(NA)
  #   if(length(which(is.na(C[,2]))) > 0.2*nrow(C)) return(NA)
  #   QC <- cbind(Q, data.frame(approx(C, xout = Q[,1]))[,2])
  #
  # }else if(is.xts(Q) & is.xts(C)){
  #   if(length(which(is.na(Q))) > 0.2*nrow(Q)) return(NA)
  #   if(length(which(is.na(C))) > 0.2*nrow(C)) return(NA)
  #   QC <- fortify.zoo(na.approx(merge(Q,C), xout = index(Q)))
  #
  # }else {
  #   stop("'Q' and 'C' must both be either data.frame or xts")
  # }
  #
  # # normalize
  # QC.norm <- data.frame(Index = QC[,1], dfNorm(QC[,2:3]))
  # QC.norm <- QC.norm[complete.cases(QC.norm),]
  # if(nrow(QC) < 9) return (NA)
  # #browser()

  # check hysteresis shape
  if(shape.check){
    # min change in C
    if(abs(max(QC[,2],na.rm = T) - min(QC[,2], na.rm = T)) <= min.dC) return(NA)

    # is it approx a loop
    ind.last <- last(which(!is.na(QC.norm[,2])))
    if( QC.norm[1,2] < 0.5){
      if (QC.norm[ind.last, 2] > (1 - min.narrow)) return(Inf)
    }else if (QC.norm[ind.last, 2] > min.narrow) return(Inf)
  }

  QC.norm <- rbind(QC.norm[QC.norm[,1] >= min.Q,], head(QC.norm[QC.norm[,1] >= min.Q,],1))
  A <-sum(diff(QC.norm[,1]) * (QC.norm[-nrow(QC.norm),2] + 0.5 * diff(QC.norm[,2])))
  return(A/(1-min.Q))
}



## Hysteresis area
# Last update: 12.02.2019
# Notes: slow
HystArea <- function(Q,C){
  #HystArea calculates the normalized area of the hysteresis loop

  #INPUT:
  #'Q' discharge data for an event, data.frame or xts
  #'C' concentration or other dependant data for the same event as 'Q'

  #VALUe:
  # numeric value of area between 0 and 1

  # create data.frame
  if(is.data.frame(Q) & is.data.frame(C)) {
    if(length(which(is.na(Q[,2]))) > 0.2*nrow(Q)) return(NA)
    if(length(which(is.na(C[,2]))) > 0.2*nrow(C)) return(NA)
    QC <- cbind(Q, data.frame(approx(C, xout = Q[,1]))[,2])

  }else if(is.xts(Q) & is.xts(C)){
    if(length(which(is.na(Q))) > 0.2*nrow(Q)) return(NA)
    if(length(which(is.na(C))) > 0.2*nrow(C)) return(NA)
    QC <- fortify.zoo(na.approx(merge(Q,C), xout = index(Q)))

  }else {
    stop("'Q' and 'C' must both be either data.frame or xts")
  }

  # normalize
  QC.norm <- data.frame(Index = QC[,1], dfNorm(QC[,2:3]))

  # shoelace formula
  n <- nrow(QC.norm)
  QC.norm <- rbind(QC.norm,QC.norm[1,])
  A = 0.5*sum(sapply(1:n, function(ii)(QC.norm[ii+1,2] + QC.norm[ii,2]) *
                       (QC.norm[ii+1,3] - QC.norm[ii,3])))

  return(A)
}


## Hysteresis lag time
# Last update: 9.1.2019

EventLag <- function(x,units = "hours", k=NA){
  # EventLag calculates lag time based on time difference between maxima of two
  # variables for an event

  #INPUT:
  #'x' bi- or multi -variate data for an event, data.frame or xts
  #'units' units of lag
  #'k' numeric smoothing window - works only for xts objects

  #VALUe:
  # numeric vector of length (nrow(x) -1) of lag time in units specified.
  # Positive lag denotes that second variable is lagging behind the first,
  # and negative the reverse.

  if(is.xts(x)){
    if(!is.na(k)) x <- rollmean2(x, k = k, shrink = T)
    ind.maxima <- apply(x,2,which.max)
    sapply(2:ncol(x),function(k) as.double(index(x)[ind.maxima[k]] -
                                             index(x)[ind.maxima[1]], units = units))

  }else if(is.data.frame(x)){
    ind.maxima <- apply(x[,-1],2,which.max)
    sapply(2:length(ind.maxima),function(k) as.double(x[ind.maxima[k],1] -
                                                        x[ind.maxima[1],1], units = units))

  }else stop("'x' is not in the right format!")

}

##----------------------------------------------------------------------------##

#' Area between two hydrograms
#'
#' @description Calculates the area between two hydrographs
#'
#' @param x,y xts or data.frame of hydrogram data
#' @param tstep numeric. Time step in seconds on which to calculate the area.
#'
#' @return An xts object of the same length as overlap of \code{x} and \code{y}.
#' @details
#'   The returned area is positive when \code{(x[t]+x[t+1])/2 > (y[t] + y[t+1])/2}
#'
#' @export
#'
HystArea3 <- function(x, y = NULL, tstep = 300){
  if(missing(y)){
    if(is.xts(x)){
      if(ncol(x) < 2) stop("Missing second time-series.")
      XY <- x[,1:2]
    }
    if(is.data.frame(x)){
      if(ncol(x) < 3)  stop("Missing second time-series.")
      if(!is.POSIXct(x[,1])) stop("Index of 'x' is not POSIXct!")
      XY <- xts(x[,2:3], order.by = x[,1], tzone = tz(x[1,1]))
    }
  }else{
    if(is.xts(x) & is.xts(y)){
      XY <- merge(x,y)
    }else if(is.data.frame(x) & is.data.frame(y)){
      x <- xts(x[,2], order.by = x[,1], tzone = tz(x[1,1]))
      y <- xts(y[,2], order.by = y[,1], tzone = tz(y[1,1]))
      XY <- merge(x,y)
    }else{
      stop("'x' and 'y' must be of the same class")
    }
  }
  #browser()

  XY.norm <- stdxts(XY)
  XY.norm <- na.approx(to_timestep(XY.norm, by = tstep))

  XY.mean <- rollmean(XY.norm, k=2, align = "left")
  XY.mean$Area <- -apply(XY.mean, 1, diff)
  return(XY.mean$Area)
}

#------------------------------------------------------------------------------#
#' Classification of Station response
#'
#' Classifies the response of a hydrological variable at a station in relation
#'   the streamflow.
#'
#' @param stats vector or a data.frame of response statistics. See Details for more.
#'
#' @details
#'   Input vector or data.frame needs to contain following statistics:
#'   \describe{
#'     \item{cor}{Pearson's correlation factor}
#'     \item{dcor}{Distance correlation factor}
#'     \item{HystIndex}{A hysteresis index (HystaAre2)}
#'     \item{ccf.value}{Maximum cross-correlation value}
#'     \item{ccf.lag}{Lag time determined by ccf}
#'   }
#'   Each row in a data.frame is considered another event.
#'
#' Responses are classified into following classes:
#' \describe{
#'   \item{1}{Perfect or near perfect matching of streamflow and the station
#'   response - line or arc shape (cor > 0.8; dcor > 0.8;-0.1 < HystIndex < 0.1).}
#'   \item{2.1}{The station response laggs after the streamflow but still close -
#'   loop or figure-of-eight shape (0.4 < cor < 0.8; 0.5 < dcor < 0.8;
#'   abs(HystIndex) > 0.1; ccf.lag < 0).}
#'   \item{2.2}{The steamflow laggs after the station response but still close -
#'   loop or figure-of-eight shape (abs(HystIndex) > 0.1; 0.5 < dcor < 0.8;
#'   0.4 < cor < 0.8; ccf.lag > 0).}
#'   \item{3}{Disconnected response - L, F, E shapes (abs(HystIndex) < 0.2;
#'   0.5 < dcor < 0.7; -0.6 < cor < 0.6; abs(ccf.value) < 0.3).}
#'   \item{4}{Indistinct response - does not fall in other categories}
#' }
#' If any statistics is NA response class is 4 and if all are NAs the returned
#' value is also NA.
#'
#' @return a factor (vector).
#'
#' @export
#'
#' @examples
#' stat1 <- c(0.7, 0.5, -0.2, 0.7, -2.5)
#' RespClassification(stats = stat1)
#'
#' stat.df <- as.data.frame(rbind(stat1, stat1 + 0.1))
#' RespClassification(stats = stat.df)
#'
RespClassification <- function(stats){
  if(!is.vector(stats)){
    return(apply(stats, 1, FUN = RespClassification))
  }else{
    if(all(is.na(stats))) {fct <- NA
    } else if(any(is.na(stats))) {fct <- 4
    } else  if(stats[1] > .8 & stats[2] > .8 & abs(stats[3]) < 0.1){ fct <- 1
    } else if(stats[1] > .4 & stats[1] < .8 & stats[2] > .5 & stats[2] < .8 &
              abs(stats[3]) > 0.1 & stats[5] < 0) {fct <- 2.1
    } else if(stats[1] > .4 & stats[1] < .8 & stats[2] > .5 & stats[2] < .8 &
              abs(stats[3]) > 0.1 & stats[5] < 0) {fct <- 2.1
    } else if(abs(stats[1]) < 0.6 & stats[2] > .5 & stats[2] < .7 &
              (abs(stats[3]) < 0.2 | is.infinite(stats[3])) & abs(stats[4]) < .3){
      fct <- 3
    } else fct <- 4

    return(factor(fct,levels = c(1, 2.1, 2.2, 3, 4)))
  }
}
#END RespClassification
#------------------------------------------------------------------------------#

#' determine the
GWLPeak <- function(x, Q.peak, npeaks = 3, nups = 24, ndowns = 4){
  #browser()
  ind.max <- which.max(x)
  time.max <- index(x)[ind.max]
  if(!is.numeric(ind.max) | length(ind.max) < 1){
    return(as.POSIXct(NA, tz = tzone(x)))
  } else if(length(which(is.na(x))) > 0.1*nrow(x)){
    return(as.POSIXct(NA, tz = tzone(x)))
  }else if(time.max %within% lubridate::interval(first(index(x)), first(index(x)) + 1800)){
    return(as.POSIXct(NA, tz = tzone(x)))
  }else if(time.max %within% lubridate::interval(last(index(x)) - 1800, last(index(x)))){
    return(as.POSIXct(NA, tz = tzone(x)))
  }else if(as.numeric(x[first(which(!is.na(x))),]) - min(x, na.rm = T) > 0.5 * diff(range(x, na.rm = T))){
    return(as.POSIXct(NA, tz = tzone(x)))
  }else{
    peaks <- pracma::findpeaks(as.numeric(x), nups = nups, ndowns = ndowns,
                               minpeakdistance = 1, npeaks = npeaks, sortstr = T,
                               minpeakheight = min(x)+.6*diff(range(x)))
    peaks.time <- index(x)[peaks[,2]]

    if(length(peaks) == 0 | sum(!is.na(Q.peak)) == 0){
      return(time.max)
    }else{
      return(c(time.max, peaks.time)[which.min(c(abs(time.max - sort(Q.peak)[1,1]),
                                                 abs(peaks.time - sort(Q.peak)[1,1])))])
    }
  }
}

#------------------------------------------------------------------------------#

PeakLag <- function(St.peak, Q.peak){
  if(is.na(St.peak) | sum(!is.na(Q.peak)) < 1) return(NA)

  as.double(sort(as.data.frame(Q.peak))[1,1] - St.peak, units = "hours")
}


#------------------------------------------------------------------------------#

InitialResp <- function(x, peak.time, nups = 18, ndowns = 6){
  if(any(is.na(x))) return(as.POSIXct(NA, tz = tzone(x)))
  tstep <- median((diff(index(x))))
  if(missing(peak.time)) peak.time <- index(x)[which.max(x)]
  if(is.na(peak.time)) peak.time <- index(x)[which.max(x)]

  ind.min <- which.min(window(x, end = peak.time))
  if(ind.min >= which.max(x)) return(as.POSIXct(NA, tz = tzone(x)))


  dales <- pracma::findpeaks(-as.numeric(window(x, end = peak.time)),
                             nups = ndowns, ndowns = nups, npeaks = 0)
  if (length(dales) > 0){
    if(length(dales) > 4)
      dales <- dales[order(dales[,4] - dales[,3], decreasing = T),]
    if (any(dales[,2] == ind.min)){
      return(index(x)[ind.min])
    }else{
      return(index(x)[dales[1,2]])
    }
  }else{
    return(index(x)[ind.min])
  }
}


#------------------------------------------------------------------------------#

#==============================================================================#
# Tracer tests ####

#' Reads procesed .ppb file from fluorometer software
#'
#' @param file Path to the .ppb file.
#'
#' @export
readPPB <- function(file){
  df <- read.csv(file = file, header = F, skip = 3, sep ="",
                 stringsAsFactors = F, skipNul =T,
                 col.names = c("ID", "Date", "Tracer1","Tracer2","Tracer3",
                               "Turbidity", "Temperature"),
                 colClasses = c("numeric", "character",rep("numeric",5)))
  df.xts <- xts(df[,3:7], order.by = as.POSIXct(df$Date, tz = "Etc/GMT-1",
                                                format = "%d/%m/%y-%H:%M:%S"))

  return(df.xts)
}

#' reads raw fluorometer file .mv
#'
#' @param file Path to the .mw file
#'
#' @export
readMV <- function(file){
  df <- read.csv(file = file, header = F, skip = 3, sep ="",
                 stringsAsFactors = F, skipNul =T,
                 col.names = c("ID", "Date", "R","Tracer1","Tracer2","Tracer3",
                               "Turbidity", "Baseline", "BatteryV",
                               "Temperature", "conductivity"),
                 colClasses = c("numeric", "character",rep("numeric",9)))
  df.xts <- xts(df[,4:11], order.by = as.POSIXct(df$Date, tz = "Etc/GMT-1",
                                                 format = "%d/%m/%y-%H:%M:%S"))

  return(df.xts)
}

#==============================================================================#
# General R functions ####

# Stops the execturion of a function without showing the error
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()
