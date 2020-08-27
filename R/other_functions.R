## Potential evapotranspiration  ###############################################

#' Potential Evapotranspiration
#'
#' Calculates the mean daily potential evapotranspiration by modified Blaney-Criddle method
#'
#' @param Temp xts or coercable to xts. At least daily air temperature.
#' @param SI vector of solar indices for each month of the year (\code{length(SI) == 12}).
#'   If NULL, the defaults for HOAL will be used.
#' @details This method uses Solar indices
#'   (\href{http://www.water-switch-on.eu/index.html}{defined by Switch On project}).
#'   Potential daily sunshine duration from
#'   Meszaroš et al. 2002, Parajka et al. 2003.
#'
#' @return xts of daily PET
#' @export
potET <- function(Temp, SI = NULL){
  if (!is.xts(Temp) || !is.zoo(Temp)) {
    try(Temp <- xts(Temp[,2], order.by = Temp[,1]))
  }

  if(is.null(SI)){
    SI <- c(0.1966, 0.2341, 0.2715, 0.3090, 0.3464, 0.3652,
            0.3558, 0.3277, 0.2856, 0.2481, 0.2060, 0.1873)
  }
  Temp <- to_timestep(Temp,by = "days")
  Temp$month <- lubridate::month(Temp)
  PET <- -1.55 + 0.96*(8.128 + 0.457 * Temp[,1])*SI[Temp[,2]]
  names(PET) <- "PET(mm/day)"
  return(PET)
}


#==============================================================================#
#' Remove outliers
#'
#' Removes outlies based on quantiles
#'
#' @param x numeric vector.
#' @param na.rm logical. Should the NAs be removed?
#' @param ... additional parameters passed to \link[stats]{quantile}.
#'
#' @return numeric vector.
#' @export
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- stats::quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * stats::IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

#==============================================================================#
#' Mode (Modus)
#'
#' @param v numeric vector.
#'
#' @return numeric
#' @export
#'
#' @examples
#' {
#' a <- c(1:10,3:15,6:2)
#' mode(a)
#' }
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#==============================================================================#
## Events Definition ####

#' Precipitation Events
#'
#' Finds the episodes of significant precipitation in the time series.
#'
#' @param prec xts of precipitation data
#' @param t numeric. Duration of precipitation effect in the units of \code{prec}.
#' @param sigprec numeric value of minimum precipitation to consider.
#' @param plot logical. Should plot be shown.
#'
#' @details   If 'sigprec' is missing set it based on the mean precipitation value.
#' @return A data frame with start time, end time and string of the period ("2017-01-03/2017-01-06").
#' @export
PrecEvent <- function(prec, t, sigprec, plot = F){
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
    dy <- dygraph(merge(PrecCum,PrecCum[Events[,3]]))%>%
      dyOptions(useDataTimezone = T, fillGraph = T)
    print(dy)
  }

  return(Events)
}

#------------------------------------------------------------------------------#
#' Precipitation Events (ver 2)
#'
#' Finds the episodes of significant precipitation in the time series.
#'
#' @param prec xts of precipitation data
#' @param duration.h numeric. Duration of precipitation effect in hours.
#' @param sigprec numeric value of minimum precipitation to consider.
#' @param plot logical. Should plot be shown.
#'
#' @details   If 'sigprec' is missing set it based on the mean precipitation value.
#' @return A data frame with start time, end time and string of the period ("2017-01-03/2017-01-06").
#' @export
PrecipitationEvent <- function(prec, duration.h, sigprec, plot = F){
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

#==============================================================================#
#' Recharge Events
#'
#' Finds the recharge events based on groundwater and precipitation data.
#'
#' @param GWL xts of groundwater data.
#' @param prec xts of rainfall data
#' @param minD numeric. Minimum change in groundwater table [m]
#' @param PrecEvents data.frame of rainfall events returned by \code{PrecEvent}
#'   or \code{PrecipitationEvent}.
#' @param plot logical. Should the plot be shown?
#'
#' @return A data frame with start time, end time and string of
#'   the period ("2017-01-03/2017-01-06").
#'
#' @export
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
# Rainfall events - taken from Rainmaker package
# https://github.com/USGS-R/Rainmaker

#' RMevents
#'
#' Rainfall event determination
#'
#' @description
#' Compute rainfall event variables based on time series of rain data with only one rain
#' gage or one mean radar rain column.
#' This function is taken from Rainmaker package (\url{https://github.com/USGS-R/Rainmaker})
#'
#' @param df dataframe with rainfall
#' @param ieHr numeric Interevent period in hours, defaults to 6,
#' @param rainthresh numeric Minimum event depth in units of the rain column, default is given as 5.1 assuming millimeters (0.2")
#' @param rain string Column name of rainfall unit values, defaults to "rain"
#' @param time string column with as.POSIXctdate, defaults to "pdate"
#' @return list of all rain events that surpass rainthresh (storms2) and all rain events (storms). Also returns all
#' a data frame of all rain observations > 0 with the associated date/time and assigned event number (tipsbystorm) and
#' the minimum time difference between observations (timeInterval)
#' @import dplyr
#' @importFrom rlang sym
#' @export
RMevents <- function (df, ieHr = 12, rainthresh = 4, rain = "rain", time = "pdate")
{
  if (!time %in% names(df)) {
    stop("Supplied 'time' column name not in df")
  }
  if (all(is.na(df[[time]]))) {
    stop("All time values are NA")
  }
  ieMin <- ieHr * 60
  df <- df[order(df[, time]), ]
  df <- unique(df)
  df <- df[df[rain] != 0, ]
  df <- df[df[rain] > 1e-05, ]
  df["event"] <- NA
  df[1, "event"] <- 1
  dif_time <- diff(df[[time]])
  timeInterval <- min(dif_time)
  df$dif_time[2:nrow(df)] <- dif_time
  for (i in 2:nrow(df)) {
    if (dif_time[[i - 1]] >= ieMin) {
      df$event[i] <- df$event[i - 1] + 1
    }
    else {
      df$event[i] <- df$event[i - 1]
    }
  }
  rain.events <- aggregate(x = df[[rain]], by = list(df$event),
                           sum)
  time_quo <- sym(time)
  start.dates <- group_by(df, event) %>% summarize(start_date = min(!!time_quo))
  start.dates <- start.dates$start_date - timeInterval
  end.dates <- group_by(df, event) %>% summarize(end_date = max(!!time_quo))
  end.dates <- end.dates$end_date
  out <- data.frame(stormnum = rain.events[, 1], StartDate = start.dates,
                    EndDate = end.dates, rain = rain.events[, 2])
  out2 <- subset(out, rain >= rainthresh, row.names = FALSE)
  return(list(storms2 = out2, storms = out, tipsbystorm = df[,
                                                             c(rain, time, "dif_time", "event")], timeInterval = timeInterval))
}

#==============================================================================#

#' Precipitation Events
#'
#' Finds the episodes of significant precipitation in the time series.
#'
#' @param prec xts of precipitation data
#' @param Q xts of discharge data
#' @param t numeric. Duration of precipitation effect in hours.
#' @param sigprec minimum event depth in units of rainfall data.
#' @param t_rec recession time in hours
#' @param plot logical. Should plot be shown.
#' @param sigQ numeric. Minimum change of discharge in units of discharge. Defaults to 0.0016 mm.(2.5 l/s)
#'
#' @return A data frame with start time, end time, total rainfall,
#' end time of discharge event and string of the period ("2017-01-03/2017-01-06").
#' @export
RainRunoffEvents <- function(prec, Q, t = 6, sigprec = 4, t_rec = 48,plot = F, sigQ = 0.0016){
  prec.df <- fortify.zoo(prec)
  t_step <- as.double(median(diff(index(prec))),units="secs")
  colnames(prec.df) = c("Index", "Rain")
  rain.events <- RMevents(prec.df, ieHr = t, rainthresh = sigprec, rain = "Rain", time = "Index")
  ev <- rain.events$storms2
  ev$EndDateQ <- ev$EndDate
  storms <- rain.events$storms
  #extending recession period
  for(ii in 1:nrow(ev)){
    ind_ev <- which(storms$StartDate == ev$StartDate[ii]) + 1
    if(ind_ev <= nrow(storms)){
      ev$EndDateQ[ii] <- min(ev$EndDate[ii] + t_rec*3600,
                             storms$StartDate[ind_ev] - t_step)
    }else {
      ev$EndDateQ[ii] <- min(ev$EndDate[ii] + t_rec*3600,
                             last(index(prec)))
    }
  }

  ev$Period <- paste0(ev$StartDate, "/",ev$EndDateQ)

  ev$Qdata <-
    sapply(ev$Period,
           function(ii){
             !any(is.na(Q[ii]))
           })
  ev$Qresponse <-
    sapply(1:nrow(ev), function(ii){
      if(ev$Qdata[ii] == T){
        if(diff(range(Q[ev$Period[ii]])) > sigQ){
          return(T)
        }else{
          return(F)
        }
      }else{
        return(F)
      }
    })

  if(plot){
    dy <- dygraphs::dygraph(merge(prec,Q))%>%
      dygraphs::dyOptions(useDataTimezone = T, digitsAfterDecimal = 5)%>%
      dygraphs::dySeries(colnames(prec), axis = "y2", stepPlot = T, fillGraph = T)%>%
      HOAL::dyShadings(ev$StartDate, ev$EndDateQ)%>%
      HOAL::dyShadings(ev$StartDate[!ev$Qresponse], ev$EndDateQ[!ev$Qresponse], color = "red")%>%
      dygraphs::dyEvent(ev$StartDate, label = ev$StartDate)
    print(dy)
  }
  return(ev)
}

#==============================================================================#
#  Start of rise of hydrograph ####
# Last Update: 14.012019
# Notes: Works bad - not to be used

#' Hydrograph rise
#'
#' @description
#'   Determines the sart time of the rise in water level response to a
#'   rainfall event by finding the first point of three consecutive increases in
#'   water level of at least \code{sig.rise}.
#'
#'
#' @param x xts of groundwater level data.
#' @param sig.rise significant rise of water level in units of \code{x}.
#' @param k numeric. Rolling mean window size.
#'
#' @return Vector of start of rise times.
#' @export
HydRise <- function(x,sig.rise = 0.001, k = 20){
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
#' RISE Method - GW recharge
#'
#' @param GWL xts of groundwater table data.
#' @param RechargeEvents data.frame returned by \code{RechargeEvents}.
#' @param Sy numeric. Specific yield [-].
#' @param pos logical. Should only positive recharge be considered.
#' @param prec (oprional) xts of precipitation. Should be on the same period as \code{GWL}.
#'
#' @return Numeric of total recharge.
#' @export
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
#' Dry period
#'
#' @param prec xts of precipitation
#' @param t numeric. Duration of precipitation effect units of \code{prec}.
#' @param mindrytime numeric. Minimum length of dry period in time-step of \code{prec}.
#' @param sigprec numeric value of minimum precipitation to consider.
#' @param plot logical. Should the plot be shown.
#'
#' @return A data.frame with start time, end time and string of the dry period ("2017-01-03/2017-01-06").
#' @export
dryPer <- function(prec,t, mindrytime=1, sigprec,plot=T){
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
#' USDA Texture calculator
#'
#' Classifies soil texture based on clay, silt and sand contens.
#'
#' @param clay,silt,sand single numeric value or one numeric vector of length 3.
#'   Clay, silt and sand content in percent.
#'
#' @details
#'   This function classifies soils based on their clay, silt and sand content.
#'   as defined by the USDA
#'
#' @return Soil texture class as factor.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   TextureCalc(40, 20, 40)
#'   TextureCalc(c(10,40,50))
#' }
TextureCalc <- function(clay,silt =NULL,sand=NULL){
  if(missing(silt) & missing(sand)){
    silt <- clay[2]
    sand <- clay[3]
    clay <- clay[1]
  }

  if(any(is.na(clay), is.na(silt), is.na(sand))) return(NA)

  if((silt + 1.5*clay) < 15){
    text <- "sand"
  }else if((silt + 1.5*clay) >= 15 & ((silt + 2*clay) < 30)){
    text <- "loamy sand"
  }else if((clay >= 7 & clay < 20) & (sand > 52) & ((silt + 2*clay) >= 30) | (clay < 7 & silt < 50 & (silt+2*clay)>=30)){
    text <- "sandy loam"
  }else if((clay >= 7 & clay < 27) & (silt >= 28 & silt < 50) & (sand <= 52)){
    text <- "loam"
  }else if((silt >= 50 & (clay >= 12 & clay < 27)) | ((silt >= 50 & silt < 80) & clay < 12)){
    text <- "silt loam"
  }else if(silt >= 80 & clay < 12){
    text <- "silt"
  }else if((clay >= 20 & clay < 35) & (silt < 28) & (sand > 45)){
    text <- "sandy clay loam"
  }else if((clay >= 27 & clay < 40) & (sand > 20 & sand <= 45)){
    text <- "clay loam"
  }else if((clay >= 27 & clay < 40) & (sand <= 20)){
    text <- "silty clay loam"
  }else if(clay >= 35 & sand > 45){
    text <- "sandy clay"
  }else if(clay >= 40 & silt >= 40){
    text <- "silty clay"
  }else if(clay >= 40 & sand <= 45 & silt < 40){
    text <- "clay"
  }else{
    text <- NA
  }
  text <- factor(text, levels = c("sand","loamy sand","sandy loam","loam",
                                  "silt loam","silt","sandy clay loam","clay loam",
                                  "silty clay loam","sandy clay","silty clay","clay"))
  return(text)
}

#==============================================================================#
# General R functions ####

# Stops the execturion of a function without showing the error
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()
