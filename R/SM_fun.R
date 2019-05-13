#==============================================================================#
#                   Soil Moisture                                           ####
#==============================================================================#


#' Total soil moisture
#'
#' @param sm xts object of soil moisture content data.
#'
#' @return xts object of total soil moisture.
#'
#' @export
#'
#'
SMtot <- function(sm){
  # if necessary convert from % to fraction
  if(mean(sm, na.rm = T) > 1){
    sm <- sm/100
  }#if
  # drop excess measurements where there are more values for single time
  sm <- make.index.unique(sm, drop = T)

  z = c(0.05,0.1, 0.2, 0.5, 1)

  # calculation for each time step
  sm_tot <- xts::xts(order.by = zoo::index(sm))
  sm_tot <- merge.xts(sm_tot, sm_tot=NA)
  for (t in 1:nrow(sm)){
    ii <- which(!is.na(sm[t,])) # working sensor
    if(length(ii)==0){next} # skip time step if all sensors are down
    # if at least one sensor works
    dz = 1 # whole profile = 1m
    for (j in 2:length(ii)) {# determine representation depth range for each additional working sensor
      dz2 <- 1 - (z[ii[j]]+z[ii[j-1]])/2
      dz[j-1] <- dz[j-1] - dz2
      dz <- c(dz,dz2)
    }#for
    sm_tot[t,] <- sum(dz*as.double(sm[t,ii]))
  }#for

  return(sm_tot)
}#function

#==============================================================================#
# LastEdit: 18.09.2018

#' Antecedente soil moisture index
#'
#' @description ASI calculates antecedent soil moisture index as total soil moisture
#'
#' @param sm soil moisture data.
#' @param depth depth to which to scale.
#' @param min_st minimum number of stations neccessary for calculation.
#'
#' @return xts object.
#' @export
#'
ASI = function(sm, depth = 0.6, min_st = 1){

  n = ncol(sm)
  z = c(0.05,0.1, 0.2, 0.5, 1)[1:n]

  #if(z[n] > depth){stop("'depth' is less than the depth of the deepest sensor!")}

  asi <- sm[,1]
  asi[,] <- NA

  ASI.row <- function(k, z,depth, min_st){
    k <- as.numeric(k)
    i <- !is.na(k) # which value is not NA
    n <- sum(i,na.rm = T) #number of values
    if(n < min_st){
      return(NA)
    }
    z2 <-z[i] #z where value
    k2 <- k[i]#k where value
    S <- 0 #sum
    S <- S + (z2[1]-0)*k2[1] #first available sensor below ground
    if(n>1){
      for(j in 2:n){
        S <- S + (k2[j]+k2[j-1])/2*(z2[j]-z2[j-1])
      }
    }
    S <- S + k2[n]*(depth-z2[n])
    #in case when deepest sensor is below 'depth'
    if(z2[n] > depth){
      k3 <- k2[n-1]+(depth-z2[n-1])/(z2[n]-z2[n-1])*(k2[n]-k2[n-1])
      S = S - (z2[n] - depth) * (k2[n] - k3) / 2
    }

    return(S)
  }

  asi[,] <- apply(sm, 1, ASI.row, z=z, depth=depth, min_st=min_st)
  colnames(asi) <- "ASI.m"

  return(asi)
}

#==============================================================================#
#' Mean Soil moisture
#'
#' @description calculates the mean non NA value of soil moisture content in the profile.
#' @param sm soil moisture data.
#'
#' @return an xts object of mean soil moisture
#' @export
#'
SMmean <- function(sm){
  sm_mean <- xts::xts(base::rowMeans(sm, na.rm = T), order.by = zoo::index(sm))
  names(sm_mean) <- "sm_mean"
  return(sm_mean)
}#function

#==============================================================================#
#' Change in soil moisture
#'
#' @description Differenciate the soil moisture time series and optionally applies smoothing.
#' @param sm soil moisture data.
#' @param smooth logical. Should the series be smoothed?
#' @param na.pad logical. Should the series be padded by NAs to original length?
#'
#' @return An xts object.
#' @export
#'
SMchange <- function(sm, smooth = F, na.pad=F){
  dsm <- diff(sm, na.pad = na.pad)
  if (is.numeric(smooth)){
    dsm <- zoo::rollapply(dsm, width = smooth, function(x){mean(x,na.rm = t)})
  }#if
  names(dsm) <- paste0("d",names(sm))
  return(dsm)
}#function
