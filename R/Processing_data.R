#==============================================================================#
## GWL processing ####
#==============================================================================#

#' Add new barometric data
#'
#' Adds new baromteric data to the end omitting the overlapping period.
#'
#' @param pOld,pNew xts objects with barometric data.
#'
#' @return xts object.
#'
#' @export
#'
pBaroUpdate <- function(pOld,pNew){
  ## Shift baro3 accordingly
  t_over <- intersect_date(pOld, pNew)
  #t_over <- index(pNew)[index(pNew[!is.na(pNew),]) %in% index(pOld)] #overlap period
  if(length(t_over) == 0){stop("'pNew' does not overlap with 'pOld'!")}
  print(paste("Overlap period:",  head(t_over,1), "-", tail(t_over,1) ))
  t_ed <- tail(t_over,1) # last index value in the overlap
  t_ed2 <- t_ed - 3600 # end of comparison period
  t_st <- t_ed2 - 3600*24 # start of comparison period
  t_per <- paste0(t_st,"/",t_ed2) # comparison period as string
  # timestep of new data
  t_new <- index(pNew)[(!index(pNew) %in% index(pOld)) &
                         index(pNew) > index(last(pOld)) ]
  d <- mean(pOld[t_per,1]-pNew[t_per, 1], na.rm=T) # mean difference

  # updated dataset
  pBaroNew <- rbind(pOld,pNew[t_new,1] + d)
  colnames(pBaroNew) <- "pBaro"
  return(pBaroNew)

}

#' Barometric compensation of Diver data
#'
#' BaroComp barometricaly compensates pressure data from GW divers
#'
#' @param diver xts object of total pressure from GW diver.
#' @param baro xts object of air (barometric) pressure time series of the same
#' or broader time period as \code{diver}.
#' @param dry logical. Should timesteps when diver is dry or almost dry be printed out
#' @param maxdiff numeric. Maximum difference of consecutive compensated data.
#'
#' @details
#'
#' @export
#'
BaroComp <- function(diver, baro, dry = F, maxdiff = 50){
  # Check periodicity
  perd <- xts::periodicity(diver)
  perb <- xts::periodicity(baro)

  if (!zoo::is.regular(diver, strict = T)) diver <- HOAL::xtsreg(diver)

  baro <- zoo::na.approx(baro, xout = zoo::index(diver))



  #
  # # Add frequency in seconds
  # perd$freqs <- as.double(perd$difftime, units = "secs")
  # perb$freqs <- as.double(perb$difftime, units = "secs")
  #
  # # Check that baro covers whole diver time period
  # if(perb$start > perd$start | perb$end < perd$end){
  #   warning("'baro' does not cover the whole time period of 'diver'/n
  #           Partial time series will be returned!")
  #   diver <- diver[index(baro)]
  # }
  #
  # # linear interpolation of barometric data to the time stamps from diver
  # if (perb$freqs != perd$freqs){
  #   baro <- na.approx(baro[paste0(perd$start, "/", perd$end)], xout = index(diver))
  # }

  # Calculates water column above diver in cmH20 by substracting diver pressure
  # by air pressure
  comp <- diver[, 1] - baro[, 1]
  comp[zoo::index(diver[is.na(diver[,1])])] <- NA
  comp[zoo::index(baro[is.na(baro[,1])])] <- NA
  comp[which(diff(comp) > maxdiff)] <- NA
  comp <- comp[paste0(perb$start, "/", perb$end)]
  colnames(comp) <- "WaterColumn(cmH2O)"

  # check for potential problems
  if(dry){
    dry.per <- which(comp < 10 & !is.na(comp))
    if (length(dry.per) > 0){
      warning(paste("diver si dry or almoste dry at following times:\n"), immediate. = T)
      cat(as.character(index(comp[dry.per])), sep="\n")
    }
  }

  return(comp)
}

#==============================================================================#

#' Reads manual measurements from Excel
#'
#' @param st Optional character vector. If not \code{NA}, it should be a GW
#'    station name.
#' @param man.meas.file character. Full path to the Excel file with manual
#'    measurements at GW stations ("D:/PhD/HOAL/raw_data/piezometer/GWstations.xlsx")
#' @param FTP logical. Should the latest Excel file be fetched from FTP server.
#'    Not yet implemented.
#'
#' @return A tibble with columns \code{Station}, \code{Date},
#' \code{h} and \code{Comment}.
#'
#' @export
#'
ManualMeasurments <- function(st = NA,
                              man.meas.file,
                              FTP = FALSE){
  # Local files
  file <- man.meas.file

  # all piezometers
  df <- data.frame(readxl::read_excel(file, sheet = "Manual_measurements",
                              na = c("-", "x", "?", "X"), skip = 3,
                              col_names = c("Date","Time","Station","h", "comment"),
                              col_types = c("date", "date", "text", "numeric", "text")))

  df$Date <- as.POSIXct(paste(df[,1], substr(df[,2],12,19)),
                        format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT-1")
  df <- df[,c(3,1,4,5)]

  # removing rows with no real data
  df <- df[!is.na(df$h) & !is.na(df$Date),]
  # ordering
  df <- df[order(df$Date),]
  # station subset
  if (!is.na(st)) df <- subset(df, Station == st)
  return(df)
}#End function ManualMeasurements

#------------------------------------------------------------------------------#

#' Find diver events
#'
#' @description
#' Extracts
#'
#' @details
#' This function that reads the Excel file with procesed manual measurements of
#' GWL from the field and prepares them for plotting with
#' \code{\link[dygraphs]{dyEvent}}.
#'
#' @param st optional character vector of a GW station names.
#' @param man.meas.file character. Full path to the Excel file with manual
#' measurements at GW stations
#' @param  diver.data.path character. Path to the directory of raw diver data
#' (MON files).
#'
#' @return Tibble with three columns:
#' \item{Station}{Name of the GW station.}
#' \item{Date}{POSIXct of the event.}
#' \item{Events}{Event description.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' DiverEvents("H01")
#'}
#'
DiverEvents <- function(st,
                        man.meas.file = "D:/PhD/HOAL/raw_data/piezometer/GWstations.xlsx",
                        diver.data.path = "D:/PhD/HOAL/raw_data/piezometer/"){


  # Manual measurements data
  df <- data.frame(readxl::read_excel(man.meas.file, sheet = "Manual_measurements",
                              na = c("-", "x", "?", "X"), skip = 3,
                              col_names = c("Date","Time","Station","h", "Event"),
                              col_types = c("date", "date", "text", "numeric", "text")))

  df$Date <- as.POSIXct(paste(df[,1], substr(df[,2],12,19)),
                        format = "%Y-%m-%d %H:%M:%S", tz = "Etc/GMT-1")
  #rearranging
  df <- df[,c(3,1,4,5)]

  ## add read-out times

  if(missing(st)){
    st <- unique(df$Station)
  }#end if

  for(ii in st){
    MONstart <- MONfiles(diver.data.path, keywords = paste0("_",ii,"_"),
                         start = "2011-01-01", end = "2080-01-01")[,2]
    if(length(MONstart) > 0){
      df <- rbind(df, data.frame(Station = ii, Date = MONstart,
                                 h = 0, Event = "MON start"),
                  make.row.names = F)
    }#end if
  }#end for


  df <- df[!is.na(df$Date),]

  # ordering
  df <- df[order(df$Date),]

  # add comments
  df[is.na(df[,4]),4] <- "read-out/man. meas."

  if(missing(st)) return(df[,c(1,2,4)])
  else return(df[df[,1] %in% st,c(1,2,4)])
}

#------------------------------------------------------------------------------#

#' Piezometer installation information
#'
#' @description installation information of all piezometers currently in use in the HOAL
#'   (ID,  East_GK_East, North_GK_East, Elevation,  Elevation_old, h_pipe)
#'
#' @param st optional. One or more station names to be included in the output.
#'    If NULL (default), all available data will be exported.
#' @param man.meas.file character. Full path to the Excel file with manual
#'    measurements at GW stations
#'
#' @return A dataframe with piezometer installation information
#'
#' @export
#'
#' @importFrom readxl read_xlsx
#'
PiezoInstall <- function(st = NULL,
                         man.meas.file = "D:/PhD/HOAL/raw_data/piezometer/GWstations.xlsx"){

  df <- readxl::read_xlsx(man.meas.file,
                          sheet = "Installation_overview",
                          skip = 3,
                          na = c("x","?","-","X","NA"),
                          col_names = c("Station", "X","Y", "Z", "PipeHeight","TotalDepth",
                                        "ScreenSt", "ScreenEnd", "DiverCurrent", "DiverType", "DiverDepth",
                                        "PipeColor","Status", "InstallationDate","Group","Comments"))



  if(!is.null(st)) df <- subset(df, Station == st)

  return(as.data.frame(df))
}

#------------------------------------------------------------------------------#

#' Save piezometer installation information to SHP
#'
#' @export
#' @param file character path where SHP should be saved (without ".shp")
#'
PiezoInstallGIS <- function(file = "D:/PhD/HOAL/GIS/piezometers/GW_stations"){
  #'piezo_install_GIS" writes piezometer installation information as shp file

  #INPUT:
  #'file' name of desired output file

  #VALUE
  #ESRI shapefile file

  data <- PiezoInstall()
  spp <- sp::SpatialPointsDataFrame(coords = data[,2:3],
                                    data = data,
                                    coords.nrs = 2:3,
                                    proj4string = sp::CRS("+init=epsg:31256"))

  suppressWarnings(rgdal::writeOGR(spp, dirname(file), basename(file),
                            driver="ESRI Shapefile", overwrite_layer = T))
}

#------------------------------------------------------------------------------#

#' Euclidean distance between stations
#'
#' @param piezo.install data.frame given by the \code{piezo_install}
#' @param filename character filepath to save data to. Not saved if NULL.
#'
#' @return data.frame
#' @export
PiezoDist <- function(piezo.install, filename = NULL){
  piezo.dist <- data.frame()
  for(i in 1:nrow(piezo.install)){
    for(j in 1:nrow(piezo.install)){
      piezo.dist[i,j] <- sqrt((piezo.install[i,2]-piezo.install[j,2])^2+
                                (piezo.install[i,3]-piezo.install[j,3])^2)
    }
  }
  rownames(piezo.dist) <- piezo.install[,1]
  names(piezo.dist) <- piezo.install[,1]

  # write to file
  if(!is.null(filename)){
    try(write.table(piezo.dist, file="spatial_data/piezo_distances.csv"))
    try(saveRDS(piezo.dist, file="Groundwater/spatial_data/piezo_distances.RDS"))
  }

  return(piezo.dist)
}

#------------------------------------------------------------------------------#

#' Diver elevation
#'
#' Function that calculates absolute diver hight from manual measurements
#'
#' @param manual data frame of manual measurements for station
#' @param diver xts object with pressure data from diver
#' @param baro xts object with pressure data from barometric sensor
#' @param g_elev numeric. ground elevation at station.
#' @param h_pipe numeric. hight of piezometer pipe above ground surface.
#'
#' @return data.frame
#' @export
#'
DiverElev <- function(manual, diver, baro, g_elev, h_pipe){

  # create xts object from manual measurements data
  df <- xts::xts(manual[3]/100, order.by = manual[[2]], tz = "Etc/GMT-1")
  # add a column with diver and baro pressure data
  maxgap <- 10
  df$diver <- na.locf(diver, xout = index(df), maxgap = maxgap)
  if(any(is.na(df$diver))){
    df$diver[is.na(df$diver), ] <- sapply(index(df[is.na(df$diver), ]),
                                          function(x){value <- stats::na.omit(head(diver[index(diver) >= x],maxgap))
                                          ifelse(length(value)>0,value[1],NA)})
  }
  df$baro  <- na.approx(baro , xout = index(df))
  # pipe height
  df$h_pipe <- h_pipe
  # ground surface elevation
  df$g_elev <- g_elev
  # ROK - top of piezometer without the cap
  df$ROK <- g_elev + h_pipe
  # add a column with water column in mH2O
  df$WC <- (df$diver - df$baro)/100
  # calculate diver elevation
  df$diver_elev <- df$ROK - df$WC - df$h
  # calculate cable length
  df$CL <- df$WC + df$h
  # convert back to data.frame
  df <- data.frame(df)
  df$comment <- manual[[4]]
  return(df)
}#end diverElev

#------------------------------------------------------------------------------#

#' Diver elevation list
#'
#' Saves diver elevation table to Excel
#'
#' @param Stations character vector. GW stations to use.
#' @param start Coercable to POSIXct. Start of investigation period.
#' @param end Coercable to POSIXct. End of investigation period.
#' @param filepath character path where Excel file should be saved.
#' @param man.measure data.frame of manual measurements. By default read by \link{ManualMeasurments}.
#' @param st.install data.frame of station installation data. By default read by \link{PiezoInstall}.
#'
#' @export
#'
DiverElevList <- function(Stations, start, end, filepath = NULL,
                          man.measure = ManualMeasurments(),
                          st.install = PiezoInstall()){

  # Review manual measurements with added diver elevations
  elev.list <- list()
  for (st in Stations){
    cat(st, "  ")
    # extract manual measurements for current station and wanted time window
    manual.temp <- man.measure[man.measure$Station == st &
                                 man.measure$Date >= start &
                                 man.measure$Date <= end, ]
    st.install.temp <- subset(st.install, Station == st)

    if(nrow(manual.temp) < 1){
      elev.list[[st]] <- NA
      next()
    }#END if

    # calculate diver elevation for manual measurements for current station
    elev.list[[st]] <- DiverElev(manual = manual.temp,
                                 diver = get(st)[,1],
                                 baro = baro,
                                 g_elev = st.install.temp$Z,
                                 h_pipe = st.install.temp$PipeHeight)
    elev.list[[st]]$median <- median(elev.list[[st]]$diver_elev, na.rm = T)
    elev.list[[st]]$sd <- sd(elev.list[[st]]$diver_elev, na.rm = T)

    outlier <- boxplot.stats(elev.list[[st]]$diver_elev)$out
    elev.list[[st]]$outlier <- elev.list[[st]]$diver_elev %in% outlier
  }#END for

  if(!is.null(filepath)){
    #save to excel
    info <- c("Divers absolute elevation calculated for all available manual measurements",
              paste("Date:", today()),
              "Created by: Lovrenc Pavlin")
    write.xlsx(info, file = paste0(filepath,"/03 manual measurements.xlsx"),
               sheetName = "info",row.names = F)
    for (st in Stations){
      write.xlsx(elev.list[[st]],
                 file = paste0(filepath,"03 manual measurements.xlsx"),
                 sheetName = st, append = T, row.names = T)
    }#end for
  }
  return(elev.list)

}#END function

#------------------------------------------------------------------------------#

#'DiverElevTable
#'
#' @description Write a table with diver elevation
#'
#' @param diver.elev.list list of diver elevations
#' @param filepath character filepath.
#' @param start start date as character or POSIXct.
#' @param end end date as character or POSIXct.
#' @param type type of output.
#'
#' @details Type of output:
#'   "all" - all available manual measurements
#'   "median" - median from all available manual measurements
#'   "one" - one elevation value for the whole period
#'
#' @export
#' @return None. Excel file.
#'
#'
DiverElevTable <- function(diver.elev.list, filepath, start, end, type = c("all","median","one")){
  Stations <- names(diver.elev.list)

  #'one' one value for the whole period
  if(type == "one"){
    diver.elev.table <- data.frame(Station = Stations,
                                   diver_elev = as.double(sapply(diver.elev.list, function(x) x$median[1])),
                                   valid_from = as.POSIXct(start, tz = "Etc/GMT-1"),
                                   valid_to = as.POSIXct(end, tz = "Etc/GMT-1"))

    #'median' median value in each of sections of the period
  }else if(type == "median"){
    diver.elev.table <- data.frame(matrix(ncol = 4, nrow = 0))
    for(st in Stations){
      valid.from <- as.POSIXct(start, tzone = "Etc/GMT-1")
      ii <- 1
      while(ii <= nrow(diver.elev.list[[st]])){
        valid.to <- as.POSIXct(rownames(diver.elev.list[[st]])[ii], tzone = "Etc/GMT-1")
        diver.elev.table <- rbind(diver.elev.table,
                                  data.frame(station = st,
                                             diver_elev = diver.elev.list[[st]]$median,
                                             valid_from = valid.from,
                                             valid_to= valid.to))
        valid.from <- valid.to
        ii <- ii + 1
      }#end while
      diver.elev.table <- rbind(diver.elev.table,
                                data.frame(station = st,
                                           diver_elev = diver.elev.table[nrow(diver.elev.table),2],
                                           valid_from = valid.from,
                                           valid_to = end))
    }#end for
    #'all' different value for each section of the period
  }else{
    diver.elev.table <- data.frame(matrix(ncol = 4, nrow = 0))
    for(st in Stations){
      valid.from <- as.POSIXct(start, tzone = "Etc/GMT-1")
      ii <- 1
      while(ii <= nrow(diver.elev.list[[st]])){
        if(!is.na(diver.elev.list[[st]][ii,"diver_elev"]) | !diver.elev.list[[st]][ii,"outlier"]){
          valid.to <- as.POSIXct(rownames(diver.elev.list[[st]])[ii], tzone = "Etc/GMT-1")
          diver.elev.table <- rbind(diver.elev.table,
                                    data.frame(station = st,
                                               diver_elev = diver.elev.list[[st]][ii,"diver_elev"],
                                               valid_from = valid.from,
                                               valid_to= valid.to))
          valid.from <- valid.to
        }#end if
        ii <- ii + 1
      }#end while
      diver.elev.table <- rbind(diver.elev.table,
                                data.frame(station = st,
                                           diver_elev = diver.elev.table[nrow(diver.elev.table),2],
                                           valid_from = valid.from,
                                           valid_to = as.POSIXct(end, tz = "Etc/GMT-1")))
    }#end for
  }#end else

  diver.elev.table$valid_from <- as.character(diver.elev.table$valid_from)
  diver.elev.table$valid_to <- as.character(diver.elev.table$valid_to)

  info2 <- c("Table of absolute diver elevations",
             paste("Date:", today()),
             "Created by: Lovrenc Pavlin",
             "Comments:",
             "Absolute diver elevations and periods of validity were taken from '03 manual measurements.xlsx' found in the same folder as this one",
             "Diver elevations are given in m above mean sea level")
  xlsx::write.xlsx(info2, file = paste0(filepath, "04 diver elevation.xlsx"),
             sheetName = "info2", col.names = F, row.names = F,append = F)
  xlsx::write.xlsx(diver.elev.table, file = paste0(filepath, "04 diver elevation.xlsx"),
             sheetName = "diver_elev",append = T, col.names = T, row.names = F)

  return(diver.elev.table)
}#end function

#============================ absGWL =====================================#

#' Absolute Groundwater level
#'
#' @description Function, that converts water column time series to absolute groundwater levels.
#'
#' @param WC xts object containing water column data of a groundwater gauge.
#'   This is equal to baro compensated diver presure. It should be in cmH2O.
#' @param diver_elev table of diver absolute elevation and its validity
#'   interval - from manually created excel file
#'
#' @return xts object
#' @export
#'
AbsGWL <- function(WC, diver_elev){

  if (ncol(WC) > 1) stop("'WC' should be univariate xts object of Water column above the diver")
  if (ncol(diver_elev) != 3) stop("'diver_elev' should have 3 columns (diver elevation, valid_from, valid_to")

  GWL <- c()
  for ( i in 1:nrow(diver_elev)){
    if( i == 1 | length(GWL)==0){
      # first interval is take as a whole
      GWL <- WC[paste0(diver_elev[i,2], "/", diver_elev[i,3]), 1]/100 + diver_elev[i,1]
      # add diver elevation and validity interval to comments
      comment(GWL) <- append(comment(GWL), paste(names(diver_elev), collapse = "   "))
      comment(GWL) <- append(comment(GWL), paste(diver_elev[i,1],
                                                 format(diver_elev[i,2:3], "%Y-%m-%d %H:%M:%S"),
                                                 collapse = " "))
      #END if
    } else {
      # consecutive intervals are take without start time stamp
      GWL <- rbind(GWL, WC[paste0(diver_elev[i,2], "/", diver_elev[i,3]), 1][-1] / 100 + diver_elev[i,1])
      # add diver elevation and validity interval to comments
      comment(GWL) <- append(comment(GWL), paste(diver_elev[i,1],
                                                 format(diver_elev[i,2:3], "%Y-%m-%d %H:%M:%S"),
                                                 collapse = " "))
    }#END else
  }#END for

  # Change time series name
  names(GWL) <- "GWL(masl)"

  return(GWL)
}#END function absGWL

# ============================= depthGWL ======================================#

#' Calculate depth to groundwater table
#'
#' @description Function converts absolute groundwater levels to depth to ground
#'   water from the ground surface.
#'
#' @param absGWL univariate xts object of absolute groundwater levels.
#' @param ground absolute elevation of ground surface at the station.
#'
#' @return an xts object.
#' @export
#'
DepthGWL <- function(absGWL, ground){

  # input check
  if (length(ground) != 1) stop("'ground' is not a number")
  if (ncol(absGWL) != 1) stop("'absGWL' is not an univariate xts object")

  # depth to ground water
  depth <- absGWL - ground

  # comments
  names(depth) <- "b(m)"
  comment(depth) <- append(comment(depth), paste("  Ground surface elevation=", ground, "  m a.m.s.l."))

  return(depth)
}#End function depthGWL

#============================= spike.rm ======================================#

#' Spike removal
#'
#' @description removes spikes in the xts objects by comparing original time series and series smoothed by rolling mean
#'
#' @param x xts object containing spikes. Only first collumn will be corrected.
#' @param d numeric value threshold difference to rolling mean.
#' @param k window length 2*k+1 in indices, default is 30.
#' @param t0 threshold, default is 3 (Pearson's rule). A high threshold makes
#'   the filter more forgiving, a low one will declare more points to be outliers.
#'
#' @note Works best on Baro compensated data - otherwise a lot of smoothing occurs.
#'
#' @return a corrected xts object.
#' @export
#'
#' @importFrom pracma hampel
spike.rm <- function(x, d, k = 30,t0 = 7){

  # indices of NA values
  na.ind <- which(is.na(x[,1]))
  # 'x' without NAs - works better with 'hampel'
  x.na.rm <- x[!x[,1] %in% x[na.ind,1],1]

  if (length(x[,1]) < (2*k+1)) stop("'x' is shorter than 2*'k'+1!")

  ## spikes in the middle of time series
  # Median absolute deviation (MAD) outlier in Time Series
  ha <- pracma::hampel(x.na.rm[,1],k, t0 = t0)
  x[zoo::index(ha$y),1] <- ha$y
  ind <- which(zoo::index(x) %in% zoo::index(ha$y[ha$ind]))

  ## spikes on the edges of time series
  # series length
  len <- length(index(x))
  # start of the series - next observation carried backward
  for (i in k:1){
    if (!is.na(x[i,1]) & abs(x[i,1]-median(x[i:(i+k),1], na.rm = T)) > d){
      x[i,1] <- NA
      x[i,1] <- zoo::na.locf(x[i:(i+k),1], fromLast = T)[zoo::index(x[i])]
      ind <- append(ind, i)
    }#END if
  }#END for


  # end of the series - last observation carried foreward
  len <- length(zoo::index(x))
  for (i in (len - k):len){
    if (!is.na(x[i,1]) & abs(x[i,1] - mean(x[(len - k):len,1], na.rm = T)) > d){
      x[i,1] <- NA
      x[i,1] <- zoo::na.locf(x[(len - k):len,1])[zoo::index(x[i])]
      ind <- append(ind, i)
    }#END if
  }#END for

  ## for diagnostics
  #dygraph(merge(mon3[,1],x[,1])) %>%
  #  dyEvent(index(x[ind]), rep("", length(ind)), labelLoc = "bottom")%>%
  #  dyOptions(useDataTimezone = TRUE)


  # comments and messages
  # ind.br <- c(0,which(diff(ind) != 1, length(ind)))
  # ind.list <-  sapply(seq(length(ind.br)-1),
  #                    function(i) ind[(ind.br[i]+1):(ind.br[i+1])])
  # for (i in 1:length(ind.list)){
  #   ind.int[i] <- ifelse(length(ind.list[[i]]) == 1,
  #                     as.character(index(x[ind.list[[i]]])),
  #                     as.character(interval(start = index(x[ind.list[[i]][1]]),
  #                              end = index(x[ind.list[[i]][length(ind.list[[i]])]]))))
  # }
  #
  # comment(x) <- append(comment(x), "  Interpolated points     =")
  # if (length(ind.int) < 20) comment(x) <- append(comment(x), paste(ind.int))
  # else {
  #   comment(x) <- append(comment(x), paste(length(ind), "points"))
  #   warning(paste(length(ind), "outlier points were found. Verify the results!"))
  # }
  return(x)

}#End function spike.rm

