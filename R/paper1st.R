#==============================================================================#
## Hysteresis paper ##
#==============================================================================#

#' Lag times
#'
#' Determins peak-to-peak lag times between all columns in \code{data} based on cross-correlation.
#'
#' @param data multivariate xts object
#' @param max.lag max lag time in days to pass to \code{\link{ccfxts}}.
#'
#' @details \code{data} should not include NA values.
#'
#' @return A data.frame with one row for each station pair (long format).
#' @note Use of this function is not recommended!
#' @examples
#' \dontrun{
#' data(sample_xts)
#' lag_times(sample_xts, max.lag = 3)
#' }
lag_times <- function(data,max.lag=1){
  data <- rollmean2(data,k=9)
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

##============================================================================##
## Hysteresis index ####
#Last update: 9.1.2019

# Literature:
#  Lloyd et al. (2016) Testing an improved index for analysing storm
#   discharge-concentration hysteresis

#' Hysteresis index - Lloyd
#'
#' @description Lloyd et al. (2016) Testing an improved index for analysing
#'   storm discharge-concentration hysteresis
#'
#' @param Q discharge data for an event, data.frame or xts.
#' @param C concentration or other dependant data for the same event as 'code{Q}.
#' @param sections numeric value of data partitioning size in percent.
#' @param Q.min numeric [0,1]. fraction of lowest \code{Q} values to exclude
#'   from the calculation. Defaults to 0.15.
#' @param plt logical. Plot results?
#'
#' @return A data.frame with HI values for each section
#' @export
#'
HI <- function(Q, C, sections = 5, Q.min = 15, plt = F){
  # NA check
  if(any(is.na(Q),is.na(C))) {
    #warning("Input data contains NAs - NA returned")
    Q.quant <- seq(from = Q.min, to = 100 - sections, by = sections)/100
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
  Q.quant <- seq(from = Q.min, to = 100 - sections, by = sections)/100
  QC.RL.quant <- data.frame(approx(QC.RL[,2:3], xout = sort(c(QC.RL[,2],Q.quant)), ties = mean))
  QC.FL.quant <- data.frame(approx(QC.FL[,2:3], xout = sort(c(QC.FL[,2],Q.quant)), ties = mean))

  if(nrow(QC.RL) < length(Q.quant) | nrow(QC.FL) < length(Q.quant)){
    return(data.frame(Q.quant, HI = NA))
  }

  Sections <- seq(Q.min, 100, by = sections)/100

  A.RL <-  sapply(2:length(Sections),
                  function(ii){
                    QC.section <- QC.RL.quant[QC.RL.quant[,1]>= Sections[ii-1] &
                                                QC.RL.quant[,1] <= Sections[ii],]
                    if(nrow(QC.section) < 2 | any(is.na(QC.section))) return(NA)
                    #QC.section
                    sum(sapply(2:nrow(QC.section),
                               function(jj){
                                 0.5 * (QC.section[jj,1] - QC.section[jj-1,1]) *
                                   (QC.section[jj,2] + QC.section[jj-1,2])

                               }))
                  }
  )

  A.FL <-  sapply(2:length(Sections),
                  function(ii){
                    QC.section <- QC.FL.quant[QC.FL.quant[,1]>= Sections[ii-1] &
                                                QC.FL.quant[,1] <= Sections[ii],]
                    if(nrow(QC.section) < 2 | any(is.na(QC.section))) return(NA)

                    sum(sapply(2:nrow(QC.section),
                               function(jj){
                                 0.5 * (QC.section[jj,1] - QC.section[jj-1,1]) *
                                   (QC.section[jj,2] + QC.section[jj-1,2])

                               }))
                  }
  )

  dA <- A.RL-A.FL

  # Hysteresis index
  HI <- data.frame(section = Q.quant, A.RL, A.FL, dA)

  #plot
  if(plt == T){
    par(mfrow= c(1,3))
    plot(QC.norm[,2:3], type = "l")
    points(QC.RL.quant[QC.RL.quant[,1] %in% Sections,], col = "blue")
    points(QC.FL.quant[QC.FL.quant[,1] %in% Sections,], col = "red")

    plot(Sections[-length(Sections)], A.RL, type = "s",col = "blue",
         ylim = c(min(A.RL, A.FL), max(A.RL, A.FL)), xlab = "Q [-]", ylab="A [-]")
    lines(Sections[-length(Sections)], A.FL, type = "s",col = "red")

    plot(Sections[-length(Sections)], dA, type="s", col = "green", ylim = c(min(dA,0), max(dA,0)))
    abline(h=0, col = "black", xlab="Q [-]", ylab="dA [-]")
  }

  return(HI)
}


#------------------------------------------------------------------------------#
## Hysteresis index Zuecco
# Reference: Zuecco G., Penna, D., Borga M., van Meerveld H.J. (2016) A versatile
#               index to characterize hysteresis between hydrological variables
#               at the runoff event timescale

#' Hysteresis index - Zuecco
#'
#' \code{HystArea2} returns the area of the Hysteresis loop between -1 and 1.
#'
#' @details
#'   shape.check was added to discarde or indicate when data is not appropriate
#'   to calculate the hysteresis area. This can happen when \code{C} only
#'   increases in the time of the event.
#'
#' Reference: Zuecco G., Penna, D., Borga M., van Meerveld H.J. (2016) A versatile
#'              index to characterize hysteresis between hydrological variables
#'              at the runoff event timescale.
#'
#'  The method is modified so that the area is calculated based on all data points.
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
#' @examples
#' \dontrun{
#' data(sample_xts)
#'
#' # two xts objects
#' Q <- sample_xts[,1]
#' C <- sample_xts[,2]
#' HystArea2(Q,C)
#'
#' # two data.frames
#' Q <- data.frame(Index = xts::index(sample_xts), sample_xts[,1])
#' C <- data.frame(Index = xts::index(sample_xts), sample_xts[,2])
#' HystArea2(Q, C, min.Q = 0.1)
#'
#' # one xts object
#' Q <-sample_xts[,1:2]
#' HystArea2(Q, min.Q = 0.1, shape.check = T)
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

#' Hysteresis index - shoelace formula
#'
#' Calculates the normalized area of the hysteresis loop using the
#'   \href{https://en.wikipedia.org/wiki/Shoelace_formula}{shoelace formula}.
#'
#' @param Q vector, data.frame or xts of event discharge data (and second variable data)
#' @param C vector, data.frame or xts of second event response variable, optional if
#' \code{Q} is appropriate object.
#'
#' @return
#'   Numeric value of area between -1 and 1.
#'   Returns Inf if data does not form a loop (only if \code{shape.check == T}).
#'   Returns NA if change of C is less or equal to min.dC
#'   (only if \code{shape.check == T})
HystArea <- function(Q,C){

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


## ===========================================================================##
# Last update: 9.1.2019

#' Hysteresis lag time
#'
#' @description Calculates lag time based on time difference between maxima of two
#'   variables for an event
#'
#' @param x bi- or multi -variate data for an event, data.frame or xts.
#' @param units character. units of lag ("secs", "hours", "mins","days")
#' @param k numeric size of smoothing window - works only for xts objects.
#'
#' @return
#'   numeric vector of \code{length = (nrow(x) -1)} of lag time in units specified.
#'   Positive lag denotes that second variable is lagging behind the first,
#'   and negative the reverse.
#' @export
#'
EventLag <- function(x,units = c("hours", "secs", "mins","days"), k=NA){
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
#' \dontrun{
#' stat1 <- c(0.7, 0.5, -0.2, 0.7, -2.5)
#' RespClassification(stats = stat1)
#'
#' stat.df <- as.data.frame(rbind(stat1, stat1 + 0.1))
#' RespClassification(stats = stat.df)
#' }
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

#' GWL peak time
#'
#' @description Determins the time of GWL peak.
#'
#' @param x single variate xts data of groundwater levels.
#' @param Q.peak POSIXct or coercable to POSIXct. Time of discharge peak.
#' @param npeaks intiger. Maximum number of peaks to return.
#' @param nups intiger. Minimum number time steps before peak when water table increases.
#' @param ndowns intiger. Minimum number of time steps after peak when water table decreases.
#'
#' @return POSIXct of peak time.
#'
#' @export
#'
#' @seealso pracma::findpeaks
GWLPeak <- function(x, Q.peak, npeaks = 3L, nups = 24, ndowns = 4){
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

#' Peak-to-peak lag time
#'
#' Returns the lag time between station and discahrge peak.
#'
#' @param St.peak,Q.peak POSIXct vector or vector coercable to POSIXct.
#'   Peak time of a station and discharge.
#'
#' @details Returned lag time is > 0 if station peaks before the dischare and < 0
#'   if the station peaks after discharge.
#'
#' @return numeric vector of \code{length(St.peak)}
#' @export
PeakLag <- function(St.peak, Q.peak){
  if(is.na(St.peak) | sum(!is.na(Q.peak)) < 1) return(NA)

  as.double(sort(as.data.frame(Q.peak))[1,1] - St.peak, units = "hours")
}


#------------------------------------------------------------------------------#

#' Initial response time
#'
#' Returns the time when station first responds to the precipitation input.
#'
#' @param x xts. Station data (e.g. groundwater table)
#' @param peak.time POSIXct or coercable to POSIXct. Time of station peak.
#' @param nups intiger. Minimum number of time steps after initial response when water table increases.
#' @param ndowns intiger. Minimum number of time steps before initial response when water table decreases.
#'
#' @details This function finds local minimums and select the one closes to the
#'   start of the event.
#'
#' @return POSIXct.
#' @export
#' @seealso pracma::findpeaks, GWLPeak
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

