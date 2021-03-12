#==============================================================================#
## xts functions ####
#==============================================================================#

# Author: Lovrenc Pavlin
# Created: 11.05.2019
# Last Update: 11.05.2019


#' Change xts' time step
#'
#' Interpolate or aggregate an xts object to desired time step.
#'
#' @param x An xts object.
#' @param by Time step as a character (e.g. "30 secs" "mins", "hours", "6 hours")
#'   or numeric number of seconds.
#' @param FUN Aggregation function (e.g. mean, sum).
#' @param force Logical. Set to TRUE to ignore \code{x} periodicity.
#' @param interpolation Type of interpoletion method: "linear" or "spline".
#' @param maxgap Intiger. How many missing time steps are still interpolated.
#'
#' @return An xts object.
#' @export
#'
#' @examples
#' \dontrun{
#' # sample data
#' data(sample_xts)
#' # Aggregation to monthly means
#' x.monthly <- to_timestep(sample_xts, by = "months")
#'
#' # Interpolation to 12 hours values
#' x.hourly <- to_timestep(sample_xts, by = "12 hours")
#' }
#'
to_timestep <- function(x, by, FUN = mean, force = F, interpolation = c("linear", "spline"), maxgap = 20L){

  #args cleanup
  interpolation <- match.arg(interpolation)
  first = xts::first
  last = xts::last

  if(!is.xts(x) | nrow(x) < 2){
    stop("'x' must be an xts object with minimum 2 observation.")
  }

  # periodicity of x as difftime object
  p <- as.difftime(median(diff(index(x))), units = "secs")

  valid = 0L
  # 'by' is not supplied or too long
  if (length(by) != 1L) {
    stop("'by' must be of length 1")

    # 'by' is numeric
  }else if (is.numeric(by)){
    #less than 60 seconds
    if (by/60 < 1){
      by2 <- c(by, "secs")
      valid = 1
      # more than 60 and less than 3600 sedonds -> minutes
    }else if(by/3600 < 1){
      by2 <- c(by%/%60, "mins")
      valid = 2

      # more than 3600 and less than 86400 sedonds <- hours
    }else if (by/86400 < 1){
      by2 <- c(by%/%3600, "hours")
      valid = 3

      # more than 86400 seconds <- days
    } else{
      by2 <- c(by %/% 86400, "days")
      valid = 4
    }
    by = as.difftime(by, units = "secs")

    # 'by' is character
  }else if(is.character(by)){
    by2 = strsplit(by, " ", fixed = TRUE)[[1L]]
    # character must consist of one or two parts
    if (length(by2) > 2L || length(by2) < 1L) {
      stop("invalid 'by' string")
    }
    # get index by matching last part of by2 to vector of time step names
    valid <- pmatch(by2[length(by2)], c("secs", "mins",
                                        "hours", "days", "weeks", "months", "years",
                                        "quarters"))

    # if units could not be matched
    if (is.na(valid)){
      stop("invalid string for 'by'")
    }

    if (length(by2)==1){
      by2 = c(1, by2[1])
    }
    # converting units
    # units are secs, mins, hours, days
    if (valid <=4 ){
      by = as.difftime(as.numeric(by2[1]), units = by2[2])

      #units are weeks, months, years, quarters
    }else{
      h = c (NA,NA,NA,1,7, 30, 365, 91)
      by = as.difftime(as.numeric(by2[1])*h[valid], units = "days")
    }


  }else {
    stop("invalid 'by'. Must be character")
  }

  ## Chose what to do based on comparrison of 'p' and 'by'
  step <- paste(by2[1], by2[2])

  ## rounding function
  rounding <- lubridate::round_date
  if (valid > 1){
    rounding <- lubridate::floor_date
  }

  # 'p' and 'by' are same
  if (p == by){
    if (force == F){
      # return unchanged 'x'
      return(x)

    }else if (force == T){
      # rounded start and end time stamps
      st <- rounding(zoo::index(xts::first(x)), step)
      ed <- rounding(zoo::index(xts::last(x)) + p - by, step)
      # endpoints
      ep <- lubridate::floor_date(zoo::index(x), step)
      # averaging to time stamp
      x_new <- xts(do.call(stats::aggregate, args = list(x = x, by = ep, FUN = FUN, na.rm = TRUE)))
    }

    # time step of 'x' is longer as specified by 'by' - interpolation
  }else if (p > by){
    # rounded start and end time stamps
    st <- rounding(index(xts::first(x)), step)
    ed <- rounding(index(xts::last(x)), step)
    g <- seq(st, ed, by = step)
    # interpolation to rounded grid
    if (interpolation == "linear") x_new <- zoo::na.approx(x, xout = g, maxgap = maxgap)
    if (interpolation == "spline") x_new <- zoo::na.spline(x, xout = g, maxgap = maxgap)


    # time step of 'x' is shorter as specified by 'by' - averaging
  }else if (p < by){
    # rounded start and end time stamps
    st <- rounding(index(xts::first(x)), step)
    ed <- rounding(index(xts::last(x)), step)
    # endpoints
    ep <- lubridate::floor_date(zoo::index(x), step)
    # averaging to time stamp
    x_new <- xts(do.call(stats::aggregate, args = list(x = x, by = ep, FUN = FUN, na.rm = TRUE)))
  }

  names(x_new) <- names(x)
  return(x_new)
}# End function to_timestep

#------------------------------------------------------------------------------#

#' intersect_date
#'
#' Returns dates that occur in two time series or date(time) vectors
#'
#' @param ts1,ts2 An xts object or a POSIXct vector
#'
#' @return Vector of class POSIXct.
#' @export
#'
#' @examples
#' \dontrun{
#' # sample data
#' data(sample_xts)
#'
#' x1 <- sample_xts["/2007-04", ]
#' x2 <- sample_xts["2007-03/", ]
#' intersect_date(x1, x2)
#' }
intersect_date <- function(ts1,ts2){
  if(xts::is.xts(ts1) & xts::is.xts(ts2)){
    lubridate::as_datetime(intersect(index(stats::na.omit(ts1)), index(na.omit(ts2))),
                           tz = xts::tzone(ts1))

  }else if(lubridate::is.POSIXct(ts1) & lubridate::is.POSIXct(ts2)){
    lubridate::as_datetime(intersect(stats::na.omit(ts1), stats::na.omit(ts2)),
                           tz = xts::tzone(ts1))
  }else stop("'ts1' and 'ts2' should be a 'xts' object or 'POSIXct vectors'.")
}


#==============================================================================#

#' Cross-correlation method for \code{xts}
#'
#' @description
#' Method of \code{ccf} function from \code{base} for \code{xts}.
#'   Returns a cross-correlation of two \code{xts} objects.
#'
#' @param x,y An object of class xts. \code{x} could be bivariate and \code{y}
#'   omitted.
#' @param type Character string giving the type of ccf to be computed. Allowed
#'   values are \code{"correlation"} (the default) or \code{"covariance"}.
#'   Will be partially matched.
#' @param plot logical. If TRUE the ccf is plotted.
#' @param lag  maximum lag  in days at which to calculate the ccf. Default is 3.
#' @param diff Logical. If TRUE the time-series will be differentiated before
#'   ccf calculation
#'
#' @return
#' An object of class \code{"acf"} (see \link[stats]{acf}) with additional elements:
#' lag - lag in hours
#' maxcor - absolutely maximal value of correlation
#' lagmax - lag where ccf is 'maxcor'
#' snames - names of input time-series generated from column names
#'
#' The lag k value returned by ccf(x, y) estimates the correlation between x[t+k] and y[t]
#' If lag < 0 y is lagging after x
#' The result is returned invisibly if plot is TRUE.
#'
#' @examples
#' t <- seq(as.POSIXct("2017-01-01"), as.POSIXct("2017-06-01"), by = "days")
#' x <- xts(runif(length(t)), order.by = t)
#' y <- xts(runif(length(t)), order.by = t)
#'
#' ccfxts(x,y, type = "correlation", plot = TRUE, lag = 3, diff = FALSE)
#'
#' @export
#'
ccfxts <- function(x,y = NULL, type = c("correlation", "covariance"),
                   plot = F, lag = 3, diff = FALSE){


  if(missing(type)) type = "correlation"
  if(!xts::is.xts(x)) stop("'x' is not an xts object!")
  if(!missing(y)){
    if(!xts::is.xts(y)) stop("'y' is not an xts object!")
    ts <- merge(x,y)
  }else{
    if(ncol(x) != 2) stop("If 'y' is missing, 'x' must have two columns!")
    ts <- x
  }

  xname <- colnames(ts)[1]
  yname <- colnames(ts)[1]

  # Combine both time series
  if (diff == T) ts <- diff(ts)

  # Time step
  tstep <- as.double(median(diff(index(ts))), units = "hours")

  # Maksimum lag to consider in index number
  max_lag <- lag*24/tstep

  # check if the provided time series are long enough
  if(all(is.na(ts[,1])) | all(is.na(ts[,2]))){
    warning("Provided time-series contains only NA values! Returned 'NA'!")
    crosscorr <- list()
    crosscorr$lag <- NA
    crosscorr$acf <- NA
    crosscorr$acf <- NA
    crosscorr$n.used <- nrow(ts)
    crosscorr$snames <- c(xname, yname)
    crosscorr$maxcor <- NA
    crosscorr$lagmax <- NA
  }else if(sum(!is.na(base::rowMeans(ts))) < max_lag){
    warning("Provided time series contain to few numeric values! Returned 'NA'!")
    crosscorr <- list()
    crosscorr$lag <- NA
    crosscorr$acf <- NA
    crosscorr$acf <- NA
    crosscorr$n.used <- nrow(ts)
    crosscorr$snames <- c(xname, yname)
    crosscorr$maxcor <- NA
    crosscorr$lagmax <- NA
  }else{
    # approximation of missing values
    ts <- zoo::na.approx(ts, na.rm = T, maxgap = 5)
    # carry on only non NA values
    ts <- ts[stats::complete.cases(ts),]

    x <- as.numeric(ts[,1])
    y <- as.numeric(ts[,2])

    # Cross correlation
    crosscorr <- stats::ccf(x, y, type = type, plot = F,
                            na.action = na.pass, lag.max = max_lag)
    # Scale lag to hours
    crosscorr$lag <- crosscorr$lag*tstep # lag in hours

    # Plotting
    if (plot) {
      layout(matrix(c(1,2), ncol=1, byrow=TRUE))
      par(mar=c(3,3,0.5,1)+0.03, mgp=c(1.5,0.3,0), tcl=-.2, xaxs="i", yaxs="r")
      plot(crosscorr, xlab = "Lag [h]")

      par(mar=c(3,3,0.5,2)+0.03)
      plot(index(ts), coredata(ts[,1]), col = "red", type = "l",
           ylab = "Groundwater level [m]", xlab = "Date", yaxt='n', xaxt='n')
      axis(2, at=pretty(range(ts[,1], na.rm=T)) ,col = "red", col.axis= "red", cex.axis = .7)
      dates <- lubridate::pretty_dates(index(ts),10)
      axis(1, dates, format(dates, "%d.%m.%Y"), cex.axis = .7)
      par(new = T, mar=c(3,3,0.5,2)+0.03)
      plot(index(ts), coredata(ts[,2]), col = "blue" , type = "l",
           axes=F, xlab = "", ylab ="")
      axis(4, at=pretty(range(ts[,2], na.rm=T)) ,col = "blue", cex.axis = .7, col.axis= "blue")
      legend("topleft", legend = c(paste0("x = ",xname),paste0("y = ",yname)), col = c("red", "blue"), lty = c(1,1))
    }

    ind.max <- which.max(abs(crosscorr$acf))
    # maximum correlation value
    crosscorr$maxcor <- crosscorr$acf[ind.max]
    # lag at max correlation value
    crosscorr$lagmax <- crosscorr$lag[ind.max]
    # change names
    crosscorr$snames <- paste0("x = ",xname," & y = ",yname)
  }
  return(crosscorr)
}

#------------------------------------------------------------------------------#
# Last update: 11.12.2018

#' IndexShift
#'
#' Shifts the index of an xts object.
#'
#' @param x An xts object.
#' @param shift Numeric value of \code{units} to shift for.
#' @param units character of the name of unit of time ("secs","mins","hours", "days",
#'   weeks", "months", "years"). Default is "secs". "months" and "years"
#'   work only with monthly and yearly data respectivly.
#'
#' @return returns an xts object of same dimensions as \code{x}
#'
#' @note This function is esspecially useful for adjusting the time zone.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_matrix, package = "xts")
#' x <- xts::as.xts(sample_matrix)
#'
#' IndexShift(x, 2, "days")
#' }
IndexShift <- function(x, shift, units = "secs"){

  if (is.xts(x) == F) stop("'x' is not an xts object.")
  if (is.numeric(shift) == F) stop("'shift' must be numeric")
  if (!units %in% c("secs", "mins", "hours", "days", "weeks", "months", "years"))
    stop("'units' is not a valid time unit name - supported units:
         'secs', 'mins', 'hours', 'days', 'weeks', 'months', 'years'")

  if (units == "secs"){ index(x) <- index(x) + shift * 1
  } else if (units == "mins"){ index(x) <- index(x) + shift * 60
  } else if (units == "hours"){ index(x) <- index(x) + shift * 3600
  } else if (units == "days"){ index(x) <- index(x) + shift * 3600 * 24
  } else if (units == "weeks"){index(x) <- index(x) + shift * 3600 * 24 * 7
  } else if (units == "months"){
    if (median(diff(index(x)), na.rm = T) < 30) stop ("'x' must contain monthly data")
    index(x) <- make_datetime(year = year(index(x)) + (month(index(x)) + shift - 1) %/% 12,
                              month = (month(index(x)) + shift - 1) %% 12 + 1,
                              day = day(index(x)),
                              tz = tzone(x))
  } else if (units == "years"){ year(index(x)) <- year(index(x)) + shift
  } else {
    stop("'units' is not a valid time unit name - supported units:
         'secs', 'mins', 'hours', 'days', 'weeks', 'months', 'years'")
  }

  return (x)
}

#------------------------------------------------------------------------------#

#' Make index unique
#'
#' Function to use on dupplicate index enties (mean, first,last)
#'
#' @param x xts object.
#' @param FUN aggregation function.
#'
#' @return xts with unique index values
#' @export
#'
#'
make.ind.unique <- function(x, FUN = "mean"){
  if(FUN == "mean"){
    ind_dupl <- unique(index(x)[duplicated(index(x))])
    if(length(ind_dupl)!=0){
      x_mean <- xts::xts(aggregate(x[ind_dupl], index(x[ind_dupl]), mean,na.rm=T))
      x <- xts::make.index.unique(x, drop = T, fromLast = F)
      #browser()
      x[ind_dupl,] <- x_mean[,]
    }
    return(x)
  }
  if(FUN == "first"){
    return(xts::make.index.unique(x, drop = T, fromLast = F))
  }
  if(FUN == "last"){
    return(xts::make.index.unique(x, drop = T, fromLast = T))
  }
}

#------------------------------------------------------------------------------#

#' Regularization of an xts object
#'
#' @description
#'   Creates a series that is ordered and has strictly regular time steps by
#'   filling missing timesteps with NAs or aggregating duplicated timesteps.
#'
#' @param x An xts object.
#' @param FUN Function for data aggregation at dupplicate timesteps.
#' @param by optional character of a format "k unit", where k is an integer and unit a valid time unit (e.g. hours)
#'
#' @return A regular xts object
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_matrix, package = "xts")
#' x <- xts::as.xts(sample_matrix)[1:20, ]
#' x.mis <- x[c(-2,-5,-7,-15,-16,-18),]
#' xtsreg(x.mis)
#'
#' x.dupl <-x[c(1,2,3,3,4:20),]
#' xtsreg(x.dupl)
#' }
xtsreg <- function(x,FUN="mean", by){
  if(missing(by)) by <-  as.numeric(median(diff(zoo::index(x))),units="secs")

  reg <- xts::xts(order.by = seq(start(x),end(x), by = by))
  x_unique <- HOAL::make.ind.unique(x, FUN = FUN)
  x_new <- xts::merge.xts(x_unique,reg)
  names(x_new) <- names(x)

  return (x_new)

}#End function xtsreg


#------------------------------------------------------------------------------#

#' Subtract series minimum/maximum
#'
#' @param x An xts object.
#'
#' @return An xts object of same dimensions as input.
#'
#' @details Works on each column seperately.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_matrix, package = "xts")
#' x <- xts::as.xts(sample_matrix)
#'
#' subMin(x)
#' subMax(x)
#' }
subMin <-function(x){
  xts(apply(coredata(x), 2, function(y) y - min(y, na.rm = T)),
      order.by = index(x), tz = tzone(x))
}

#' @describeIn subMin Subtract series maximum from the series
subMax <-function(x){
  xts(apply(coredata(x), 2, function(y) y - max(y, na.rm = T)),
      order.by = index(x), tz = tzone(x))
}

#------------------------------------------------------------------------------#
#' Create empty xts
#'
#' @description Creates an empty xts object or xts filled with NA spanning from
#'   'start' to 'end' and with timestep 'by'.
#'
#' @param start datetime of the first timestep as object coercable to POSIXct.
#' @param end datetime of the last timestep as object coercable to POSIXct.
#' @param by timestep given as number of seconds or appropriate character time descriptor.
#' @param dim number of columns to create.
#' @param tz Time zone.
#'
#' @return xts object of dimensions (\code{start}-\code{end})/\code{by} x \code{dim}
#'
#' @export
#'
#' @examples
#' emptyxts("2017-01-01", "2018-02-01", "2 days")
emptyxts <- function(start,end, by, dim = 1, tz = "Etc/GMT-1"){
  start = as.POSIXct(start, tz = tz)
  end = as.POSIXct(end, tz= tz)
  g <- seq(start, end, by = by)
  x <- xts::xts(order.by = g)
  if(dim == 0){
    return(x)
  }
  else{
    for(i in 1:dim){
      x <- merge.xts(x,NA)
    }
    return(x)
  }
}#End function emptyxts

#------------------------------------------------------------------------------#
#' rowMeans.xts
#'
#' Calculates the mean of each row of a multivariate xts object.
#'
#' @param x xts object.
#' @param na.rm logical. Should missing values (including NaN) be removed?
#'
#' @return The arithmetic mean of each row in 'x'
#'  as a xts object of the same length as 'x'
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_matrix, package = "xts")
#' x <- xts::as.xts(sample_matrix)
#'
#' rowMeans.xts(x)
#' }
rowMeans.xts <- function(x, na.rm = T){
  x_new <- xts::xts(rowMeans(x,na.rm = na.rm),
                    order.by = index(x),
                    tzone = xts::tzone(x))
  colnames(x_new) <- paste0(colnames(x)[1],"_mean")
  return(x_new)
}#End function rowMeans.xts

#------------------------------------------------------------------------------#
#' Loewss smoothing for xts
#'
#' Smoothes an xts time-series using LOWESS smoothing (\code{stats::lowess})
#'
#' @param x An (multivariate) xts object
#' @param span Optional. the smoother span. This gives the proportion of points
#' in the plot which influence the smooth at each value.
#' Larger values give more smoothness.
#' @param max.span optional numeric. Maximum allowed span. Used to prevet from over smoothing.
#' @param degree the degree of the local polynomials to be used. It can ben 0, 1 or 2.
#' @param ... Additional parameters to \link[fANCOVA]{loess.as}.
#'
#' @details If \code{span} is not given, \code{stats:optim} is used to determine,
#' the optimal value that minimizes SSE.
#'
#' \code{...} can contain additional arguments to \link[fANCOVA]{loess.as}
#' (see examples).
#'
#' @return An xts object of the same dimensions as \code{x}
#'
#' @examples
#' \dontrun{
#' x <- xts::xts(rnorm(231)+10, order.by = as.Date(13514:13744, origin="1970-01-01"))
#'
#' # user-defined span
#' x.smth.05 <- loess.xts(x, span = 0.5)
#'
#' # Automatic span selection
#' x.smth <- loess.xts(x)
#'
#' # Multivariate xts object
#' x2 <- merge(x, x+rnorm(231))
#' x2.smth <- loess.xts(x2)
#'
#' # Additional parameters
#' x2.smth <- loess.xts(x2, criterion = "gcv", degree = 2, plot = T )
#'}
#' @export
loess.xts <- function(x, span, max.span, degree = 2, ...){
  if(!xts::is.xts(x)) stop("'x' must be an xts object.")
  if(nrow(x) < 10 | length(which(!is.na(x))) < 10) stop("'x' is too short.")

  x_new <- xts::xts(order.by = index(x))
  if(missing(span)) span <- NULL

  for(ii in 1:ncol(x)){
    if(any(is.na(x[ ,ii]))){
      warning(paste("'x' contains missing values! (Name =", colnames(x[ ,ii]),
                    ", at", xts::first(index(x)),". NA returned."))
      x_new <- merge.xts(x_new, NA)
    }else{
      temp <- fANCOVA::loess.as(index(x[ ,ii]), x[ ,ii],
                                user.span = span, degree = degree, ...)
      if(!missing(max.span)){
        if(temp$par$span > max.span & is.null(span)){
          temp <- loess.as(index(x[ ,ii]), x[ ,ii],
                           user.span = max.span, degree = degree, ...)
          warning("Calculated span > max.span. Recalculating with 'max.span'.")
        }
      }
      x_new <- merge.xts(x_new, temp$fitted)
    }
  }
  colnames(x_new) <- colnames(x)

  return(x_new)
}

#------------------------------------------------------------------------------#
#' Expand or shrink subset period
#'
#' @description
#' Modifies provided subset period (e.g. "2017-01-15/2018-03-20") by extending or
#' shrinking it by 'by' number of seconds.
#'
#' @param period A character subseting period (see
#'   href{https://cran.r-project.org/web/packages/xts/xts.pdf#Rfn..parseISO8601.1}
#'   {xts::.parseISO8601}).
#' @param by A numeric value. Number of seconds to add or subtract from
#'   \code{period}. If \code{by} is a vector with length 2, first value is applied
#'   to left side and second value to the right side of \code{period}.
#'
#' @return A character vector of \code{period} format.
#'
#' @examples
#' # Expand for one day
#' expandPeriod(period = "2017-01-15/2018-03-20", by = 3600*24)
#'
#' # Shrink for 2 days
#' expandPeriod(period = "2017-01-01/2017-02-01", by = 2 *3600* 24)
#'
#' # Expand for 2 hours at the begining
#' expandPeriod(period = "2017-03-04 03:20/2017-03-06 05:00", by = c(7200,0))
#'
#' @export
#'
#' @note end time might change for one timestep
expandPeriod <- function(period, by){
  per <- xts::.parseISO8601(period)
  if(all(is.na(per))){
    stop("Provided 'period' is not in a recognisable format!")
  }

  if(!length(by) %in% c(1,2)){
    stop(paste0("'by' must be of length 1 or 2, insted ", length(by),"!"))
  }
  period_new <- ""
  if(!is.na(per$first.time)){
    per$first.time <- per$first.time + by[1]
    period_new <- paste0(per$first.time, "/")
  }else{
    period_new <- "/"
  }

  if(!is.na(per$last.time)){
    per$last.time <- per$last.time + tail(by,1)
    period_new <- paste0(period_new,per$last.time)
  }
  return(period_new)
}

################## Time series decomposition ###################################

#' Seasonal Decomposition of Time Series by Loess for xts
#'
#' @inherit stlplus::stlplus description
#'
#' @param x univariate xts object.
#' @inheritParams stlplus::stlplus
#' @param ... additional parameters
#'
#' @inherit stlplus::stlplus details
#' @return An xts object containing columns:
#' \describe{
#'   \item{raw}{Input series.}
#'   \item{seasonal}{Seasonal component of the series,}
#'   \item{trend}{Trend component of the series,}
#'   \item{remainder}{Remaining signal.}
#' }
#' @export
#'
#' @importFrom stlplus stlplus
#'
#' @seealso stlplus::stlplus
#'
#' @examples
#' \dontrun{
#' data(sample_xts)
#'
#' stlplus.xts(sample_xts, n.p = 10, s.window = 100)
#' }
#'
stlplus.xts <- function(x, n.p, s.window, s.degree = 1, t.window = NULL,
                        t.degree = 1, fc.window = NULL, fc.degree = NULL, fc.name = NULL,
                        l.window = NULL, l.degree = t.degree, s.jump = ceiling(s.window/10),
                        t.jump = ceiling(t.window/10), l.jump = ceiling(l.window/10),
                        fc.jump = NULL, critfreq = 0.05, s.blend = 0, t.blend = 0,
                        l.blend = t.blend, fc.blend = NULL, inner = 2, outer = 1,
                        sub.labels = NULL, sub.start = 1, zero.weight = 0.000001,
                        details = FALSE,...){


  # Check for univariaty
  if(ncol(x)!=1){
    stop("'x' is not univariate!")
  }

  # Periodicity
  per <- median(diff(index(x)))
  freq <- 1/as.double(per, units = "secs")

  # 'ts' object
  ts_x <- ts(as.double(x), frequency = freq)

  # decomposition
  stl <- stlplus(ts_x, n.p =n.p, s.window = s.window, s.degree = s.degree,
                 t.window = t.window, t.degree = t.degree, fc.window = fc.window,
                 fc.degree = fc.degree, fc.name = fc.name, l.window = l.window,
                 l.degree = l.degree, s.jump = s.jump, t.jump = t.jump,
                 l.jump = l.jump, fc.jump = fc.jump, critfreq = critfreq,
                 s.blend = s.blend, t.blend = t.blend, l.blend = l.blend,
                 fc.blend = fc.blend,inner = inner, outer = outer,
                 sub.labels = sub.labels, sub.start = sub.start,
                 zero.weight = zero.weight, details = details)

  stl$data <- xts(stl$data[,1:4], order.by = index(x))
  return(stl)
}


################.

#' Seasonal decomposition by Moving Averages for xts
#'
#' @description
#'   Decompose a time series into seasonal, trend and irregular components using
#'   moving averages. Deals with additive or multiplicative seasonal component.
#'
#' @param x An univariate xts object.
#' @param type The type of seasonal component. Can be abbreviated.
#' @param units character time descriptor.
#' @param plot logica. Should the result be ploted?
#'
#' @return An xts with raw data column and trend, seasonal and random components
#'   column.
#'
#'
#' @inherit stats::decompose details
#'
#' @examples
#' \dontrun{
#' data(sample_matrix, package = "xts")
#' x <- xts::as.xts(sample_matrix)[,1]
#'
#' decompose.xts(x)
#'}
decompose.xts <- function(x, type = c("additive", "multiplicative"), units = "days", plot = T){

  decomp_xts <- x
  if (missing(type)) type = "additive"

  ## Time series preparation
  # chechk for NA
  if(length(which(is.na(x))) > 0){
    x <- zoo::na.approx(x, maxgap = 13)
    if(length(which(is.na(x))) > 0) {stop("'x' contains data gaps longer than 13 timesteps.")}
    warning("Data gaps in 'x' were removed by linear interpolation")
  }
  # Periodicity
  per <- xts::periodicity(x)
  freq <- 1/as.double(per$difftime, units = units)


  # 'ts' object
  ts_x <- stats::ts(x, frequency = freq)

  ## decomposition
  decomp_x <- stats::decompose(ts_x,type = type)
  # back to xts object
  df <- data.frame(trend = decomp_x$trend,
                   seasonal = decomp_x$seasonal,
                   rand = decomp_x$random)
  decomp_xts <- merge(decomp_xts, xts(df, order.by = index(x), tzone = tzone(x)))
  names(decomp_xts) <- c(names(x), "trend", "seasonal", "random")
  if(plot == T) plot(decomp_x)

  return(decomp_xts)
}
#------------------------------------------------------------------------------#

#' Plot decomposed xts
#'
#' @param x decomposed xts object.
#'   Must contain 4 columns (observed data, trend, seasonal, random)
#' @param type Plot type: zoo or dygraph
#' @param main plot title.
#' @param yax.flip logical. flip y axis?
#' @param xlab x axis label.
#' @param ... optional parameters for plot.zoo
#'
#' @return None
#'
#' @importFrom  htmltools browsable tagList
#'
plot.decompose.xts <- function(x, type = c("zoo", "dygraph"),
                               main = "Decomposed time series",
                               yax.flip = T, xlab = "Date",...){
  #x  deomposed 'xts' time series. Must contain 4 columns (observed data, trend, seasonal, random)
  #...  optional parameters for plot function

  if(ncol(x) != 4) stop(paste("Invalid number of columns in 'x'. Expected 4, got", ncol(x)))
  if(missing(type)) type = "zoo"

  if(type == "zoo") plot.zoo(x,main = main, xlab = "Date",yax.flip = T, ...)

  if(type == "dygraph"){
    # obs <- dygraph(x, main = main, ylab = "GWL (m a.s.l.)" ) %>%
    #     dyOptions(useDataTimezone = T) %>%
    #     dySeries(names(x)[3], axis = "y2") %>%
    #     dySeries(names(x)[4], axis = "y2") %>%
    #     dyAxis("y2", label = "Amplitude (m)")
    #
    # create a list of dygraphs objects
    dy_graph <- list(
      dygraphs::dygraph(x[,1:2], group="decomp", main=main, ylab = "GWL (m a.s.l.)", height = "300px")%>%
        dyOptions(useDataTimezone = T)%>%dyCrosshair(direction="vertical"),
      dygraphs::dygraph(x[,3], group="decomp", main="seasonal", ylab = "Amplitude (m)", height = "100px")%>%
        dyOptions(useDataTimezone = T)%>%dyCrosshair(direction="vertical"),
      dygraphs::dygraph(x[,4], group="decomp", main="random", ylab = "Amplitude (m)", height = "200px")%>%
        dyOptions(useDataTimezone = T)%>%dyCrosshair(direction="vertical")
    )  # end list

    # render the dygraphs objects using htmltools
    browsable(tagList(dy_graph))

  }#END if
}

#==============================================================================#
#' Detrend an xts series
#'
#' @inheritParams decompose.xts
#'
#' @details uses classical seasonal decomposition by moving averages to remove
#'  the trend from the series.
#' @seealso \code{\link{decompose.xts}}, \code{\link[stats]{decompose}},
#'  \code{\link[HOAL]{stlplus.xts}}
#' @return detrended xts object.
#'
#'
detrend.xts <- function(x, type = c("additive", "multiplicative"), units = "days", plot = T){
  x1 <- decompose.xts(x, type = type, units, plot)
  return(x1[,1]-x1[,2])
}

#==============================================================================#
# Date of min, max ##
#' Date of minimum/maximum
#'
#' Returns the date of the minimum/maximum in the provided 'xts' object 'x'.
#'
#' @param x (multivariate) xts object.
#'
#' @return vector of minima/maxima dates (POSIXct) with the length of \code{ncol(x)}
#'
#' @export
#' @examples
#' \dontrun{
#' x <- data(sample_matrix, package = "xts")
#'
#' #univariate xts
#' date.min(x[,1])
#'
#' # multivariate xts
#' date.max(x[,])
#' }
date.min <- function(x){
  d <- c()
  for(i in 1:ncol(x)){
    x1 <- x[!is.na(x[,i]),i]
    d <- append(d, index(x1[which.min(x1)]))
  }

  return(d)
}

#' @describeIn date.min Returns the date of the maxima.
#' @export
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
#' Rollmean
#'
#' Calculates rolling mean allowing NAs in the data.
#'
#' @param x (multivariate) xts object.
#' @param k numeric. window size.
#' @inheritParams stats::filter
#' @param shrink logical. Should the leading and trailing NAs be removed?
#'
#' @return xts of dimensions as \code{x} if \code{shrink == F}, and of
#'   \code{nrow(x) - (k-1) x ncol(x)} if \code{shrink == T}.
#' @export
#' @seealso stats::filter
#' @examples
#' \dontrun{
#' data(sample_matrix, package = "xts")
#' x <- rollmean2(sample_matrix["2017",], k = 7)
#' }
rollmean2 <-
  function(x,
           k = 3,
           sides = 2,
           method = "convolution",
           shrink = F) {
    #rollmean2 calculates rolling mean of timeseries with NAs
    if (sides == 2 &
        k %% 2 == 0)
      warning("if sides = 2, 'k' should be odd!")
    if (shrink == F) {
      return(as.xts(apply(x, 2, function(z)
        stats::filter(z, rep(1 / k, k), sides = sides, method = "convolution")),
        order.by = index(x)))
    }
    if (shrink == T) {
      return(as.xts(apply(x, 2, function(z)
        stats::filter(z, rep(1 / k, k), sides = sides, method = "convolution")),
        order.by = index(x))[ceiling((k - 1) / 2 + 1):(nrow(x) -
                                                         floor((k - 1) / 2)), ])
    }
  }

#==============================================================================#
# Data normalization ####

#' normalize data.frame to the interval [0,1]
#'
#' @param x data.frame of numeric data.
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_xts)
#' df <- fortify.zoo(sample_xts["2017",c(1:3)])
#' dfNorm(df)
#' }
dfNorm <- function(x){
  as.data.frame(lapply(x, function(y) {(y - min(y, na.rm=TRUE))/(max(y,na.rm=TRUE) -
                                                                   min(y, na.rm=TRUE))}))
}

# 'stdx' function from hydroTSM
#' Normalize xts to the [0,1]
#'
#' Applies \link[hydroTSM]{stdx} to all columns of the xts provided.
#'
#' @param x xts object
#'
#' @return xts object
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_xts)
#' stdxts(sample_xts["2017-01-01",1:2])
#' }
stdxts <- function(x){
  return(xts(apply(x, 2, hydroTSM::stdx), order.by = index(x)))
}

# ## subtract min or mean value of the time-series
#
# XtsChange <- function(x, FUN = min){
#   xts(vapply(x, function(col) col - FUN(col, na.rm = T),
#              FUN.VALUE = numeric(nrow(x))),
#       order.by = index(x))
# }

##============================================================================##
#' rbind at specific time
#'
#' @param x,y (multivariate) xts objects.
#' @param at datetime value or vector coercable to POSIXct when binding should a occur.
#'
#' @return xts object with values from \code{x} on time interval
#'   \code{[start(x),at) and values from y on time interval [at,end(y)]}.
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_xts)
#' x <- sample_xts["2017-01/2017-06",1]
#' y <- sample_xts["2017-03/2017-08",1]
#' rbind_at(x,y, at="2017-04-23")
#' }
rbind_at <- function(x,y, at){
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
    temp <- rbind(x[paste0("/",at[i]-by),i],y[paste0(at[i],"/"),i])
    xy <- merge(xy,temp)
  }
  names(xy)<-names(x)
  return(xy)
}

##============================================================================##
#' Mean monthly values
#'
#' @description Calculates average value for each month (e.g. mean of all Aprils in dataset).
#'
#' @param data (multivariate) xts object.
#' @param plot logical. Should the plot be shown.
#'
#' @return A data.frame of dimensions  \code{12 x ncol(data)}.
#' @export
#'
#' @examples
#' \dontrun{
#' data(sample_xts)
#' MonMean(sample_xts)
#' }
MonMean <- function(data, plot = F){
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

  return(data.frame(monmean))
}
