#==============================================================================#
## Diver data ####
#==============================================================================#

#' Read .MON file
#'
#' Read .MON files produced by VanEssen Diver as a xts object
#'
#' @param filename Path to the .MON file to be read.
#' @param tzone Time zone of the data. Defaults to "Etc/GMT-1" (UTC+1).
#'
#' @return An \code{xts} object with two columns: pDiver - pressure [cmH20];
#' TDiver - temperature [*C]
#'
#' @export
#'
MONread <- function(filename, tzone = "Etc/GMT-1"){
  # reading MON file to data frame
  df <- data.table::fread(file = filename, header = F, data.table = F,
              skip = 54, nrows = (length(readLines(filename))-54-1))

  # column names
  colnames(df) <- c("date", "time", "pDiver", "TDiver")
  # change value type
  df$pDiver <- as.numeric(df$pDiver)
  df$TDiver <- as.numeric(df$TDiver)
  # add a column with combined date and time (as POSIXct)
  df$datetime <- as.POSIXct(paste(df[,1],df[,2]), tz = tzone,
                            format = "%Y/%m/%d %H:%M:%S")

  #convert data frame to xts object
  df <- xts::xts(df[,3:4], order.by = df[,5], tzone = tzone)

  # adding metadata
  attr(df, "comment") <- base::readLines(filename, n = 51)

  # returns the xts object
  return(df)
}

#==============================================================================#

#' Check xts from .MON file
#'
#' @description
#' Function that check if the \code{xts} object produced from MON file has
#' constant time step and if it falls to full minutes.
#'
#' @param mon An \code{xts} object created by \code{\link{MONread}}.
#'
#' @return An \code{xts} object with same number of columns as \code{mon}.
#'
#' @export
MONcheck <- function(mon){

  # check periodicity of the xts object. Produces a list object
  per <- xts::periodicity(mon)

  ## Round time stamps (e.g. 15:10:00, 15:15:00)

  # rounded first time stamp of the object to the nearest periodicity time stamp
  # e.g. "2016-12-20 12:13:09 +01" + periodicity "5 mins" -> "2016-12-20 12:15:00 +01"
  st <- lubridate::round_date(index(mon[1]), paste(per$frequency, per$units))
  #rounded last time stamp
  ed <- lubridate::round_date(index(xts::last(mon)), paste(per$frequency, per$units))

  # check if time stamps are not already rounded
  if (zoo::index(mon[1]) != st){
    # rounded time stamps grid
    g <- seq(st, ed, by = paste(per$frequency, per$units))
    # linear interpolation to rounded grid
    mon <- zoo::na.approx(mon, xout = g, maxgap = 5)
    # round values
    mon[, 1] <- round(mon[, 1], digits = 1) #pDiver
    mon[, 2] <- round(mon[, 2], digits = 2) #TDiver
  }

  return(mon)
}

#==============================================================================#

#' Combining two xts objects
#'
#' Function that combines two xts objects created from MON files
#'
#' @param mon1,mon2 An \code{xts} object.
#' @param FUN A function for data aggregation (e.g. \code{mean} or \code{sum})
#'
#' @return \code{xts} objects.
#'
#' @details Intention of this function is to add combine two \code{xts} object
#' created from consecutive MON files of a Diver.
#'
#' The function converts both data sets to a common time step determined from
#' \code{mon1}. If the time step of any datasets is 5 min this time step is
#' carried forward. The function also tries to fill the gap between two datasets.
#'
#' Combined dataset is regularized using \code{\link{xtsreg}}.
#'
#' @examples
#' \dontrun{
#'   H01.2017 <- MONread(H01_201601010000_201701010000.MON)
#'   H01.2018 <- MONread(H01_201701010000_201801010000.MON)
#'   H01 <- MONagg(mon1, mon2, FUN = mean)
#' }
#'
#' @export

MONagg <- function(mon1, mon2, FUN = mean){
  per1 <- as.double(median(diff(index(mon1)),na.rm = T), units = "secs")
  per2 <- as.double(median(diff(index(mon2)),na.rm = T), units = "secs")

  start1 <- zoo::index(xts::first(mon1))
  start2 <- zoo::index(xts::first(mon2))

  end1 <- zoo::index(xts::last(mon1))
  end2 <- zoo::index(xts::last(mon2))


  # Different dimensions
  if (ncol(mon1) != ncol(mon2)) stop (paste("Dimension error:/n
                                            'mon1' has", ncol(mon1), "columns,
                                            'mon2' has", ncol(mon2), "columns."))

  # Overlaping time-series
  if( start1 == start2) stop(paste("'mon1' and 'mon2' have same start time:",start1))
  if (end1 > start2) stop(paste("Error: 'mon2' starts (",start2,")
                                before 'mon1' ends (",end1,")!"))

  # different periodicity
  if (per1 != per2){
    # interpolate to prefered periodicity of 5 min
    if (any(per1 == 300 | per2 == 300)) {
      x <- HOAL::to_timestep(rbind(mon1, mon2), by = 300, FUN = FUN)
    }else {
      x <- HOAL::to_timestep(rbind(mon1,mon2),
                             by = max(c(per1, per2)),
                             FUN = FUN)
    }
  }else x <- rbind(mon1, mon2)

  # Interpolation of missing values
  x <- zoo::na.approx(x, maxgap = round(3600/max(c(per1, per2))))

  # Regularize
  if(!zoo::is.regular(x, strict = T)) x <- HOAL::xtsreg(x)

  return(x)
}

#------------------------------------------------------------------------------#

#' Search for MON files
#'
#' Returns a list of .MON files in the \code{path} based on the search criteria
#'
#' @param path A character of full path name.
#' @param keywords An optional character vector with one or more search keywords
#'   to be used in file search (e.g. diver code).
#    File will be used if all keywords are found. Case sensitive.
#' @param start,end optional character or POSIXct in "ymd" format. Filters found
#' files to only include those containing measurements after \code{start} and
#' before \code{end}.
#'
#' @return Data frame with file:
#' \item{}{Full file path.}
#' \item{}{Start datetime of measurements in the file}
#' \item{}{End datetime of measurements in the file}
#'
#' @examples
#' \dontrun{
#'  MONfiles(path = "GWdata", keywords = "H01",
#'            start = "2017-01-01", end = "2018-01-01")
#' }
#'
#' @export
#'
MONfiles <- function(path,keywords,start,end){
  if(missing(path))stop("'path' is missing and has no default.")
  if(missing(keywords)) keywords <- c()
  if(missing(start)) start <- as.POSIXct("1900-01-01", tz="Etc/GMT-1")
  if(missing(end)) end <- as.POSIXct("2900-01-01", tz="Etc/GMT-1")
  start <- as.POSIXct(start, tz="Etc/GMT-1")
  end <- as.POSIXct(end, tz = "Etc/GMT-1")
  #if(substr(path, nchar(path), nchar(path))  != "/") path <- paste0(path, "/", collapse = "")
  if(substr(path, nchar(path), nchar(path))  == "/") path <- substr(path, 1, nchar(path)-1)

  # file list
  files <- data.frame(filenames = list.files(path, pattern = paste0("(", keywords, ").*([[:digit:]]{12}.MON)$"),
                                             recursive = T, full.names = T),
                      stringsAsFactors = F)
  if(nrow(files) == 0){
    warning(paste("no files found matching the search parameters: path =", path,
                  "; keywords =", keywords, "; start =", start, "; end =", end))
    #stopQuietly()
    return()
  }
  files$start <- as.POSIXct(regmatches(x = files$filenames,
                                       regexpr("[[:digit:]]{12}(?=_till_)", files$filenames, perl = T)),
                            format = "%Y%m%d%H%M", tz = "Etc/GMT-1")
  files$end <- as.POSIXct(regmatches(x = files$filenames,
                                     regexpr("[[:digit:]]{12}(?!_till_)", files$filenames, perl = T)),
                          format = "%Y%m%d%H%M", tz = "Etc/GMT-1")

  # sorting
  files <- files[order(files$start), ]

  # extract only files within desired period
  files <- files[files$end > start & files$start < end,]

  return (files)

}# End function MONfiles

#------------------------------------------------------------------------------#

#' Importing multiple MON files
#'
#' MONimport reads and aggregates data of all MON files find by search citeria.
#'
#' @inheritParams MONfiles
#'
#' @return An \code{xts} object of all the measurements found in the \code{path}
#'
#' @details This function uses \code{\link{MONfiles}} to find MON files,
#' \code{\link{MONread}} and \code{\link{MONcheck}} to read them and finally
#' \code{\link{MONagg}} to combine them to one file.
#'
#' @examples
#' \dontrun{
#'   H01.2017 <- MONimport(path = "GWdata", keywords = "H01",
#'                         start = "2017-01-01", end = "2018-01-01")
#'   summary(H01.2017)
#'}
#'
#'@export

MONimport <- function(path, keywords, start, end){
  # input data
  files <- MONfiles(path, keywords, start, end)

  # remove duplicated files
  files <- files[!duplicated(files$end),]

  # read and aggregate data

  if(nrow(files) > 0){
    # 'x' is an 'xts' object
    x <- MONcheck(MONread(files[1, 1]))
    j = 2

    # Aggregate data
    while(j <= nrow(files)){
      x <- MONagg(x, MONcheck(MONread(files[j,1])))
      j = j + 1
    } # END while

    # Keep only data within the time interval 'date_int'
    x <- x[index(x) >= start & index(x) <= end]

    # regulize
    x <- xtsreg(x)

    return(x)

  }# END if
}# END function MONimport

#------------------------------------------------------------------------------#

#' Shifts the time-steps in MON data file
#'
#' @param filename character or character vector of paths to MON files to change.
#' @param shift time shift in seconds.
#' @param unit unit of time as character. Defaults to "hours".
#'
#' @return A MON file with changed data and a backup of the original data.
#' @export
#'
MONshift <- function(filename, shift, unit = "hours"){
  for(i in 1:length(filename)){
    # read the MON file
    ln <-readLines(filename[i])
    len <- length(ln)
    # line where data starts
    data_ln <- grep(pattern = "(\\[Data\\])", x = ln) + 2
    # number of lines of data
    data_len <- as.numeric(ln[data_ln - 1])

    # read the index and shift it
    ind <- as.POSIXct(substr(ln[data_ln:(data_ln + data_len -1)],1,19), tz= "Etc/GMT-1")
    ind <- ind + as.difftime(shift, units = unit)

    # replace the index with the shifted one
    substr(ln[data_ln:(data_ln + data_len -1)],1,19) <- as.character(ind, format ="%Y/%m/%d %H:%M:%S")

    # New filename
    filename_new <- paste(grep(x = strsplit(basename(filename[i]), split = "_")[[1]],
                               pattern = "([[:alpha:]+])([[:digit:]+])", value = T),
                          collapse = "_")
    start <- format(ind[1], format = "%Y%m%d%H%M")
    end <- format(ind[data_len], format = "%Y%m%d%H%M")

    ## Replace timesteps in file header
    # filename
    ln_fn <- grep(pattern = "(FILENAME)", x = ln)
    substr(ln[ln_fn],
           start = stringr::str_locate(ln[ln_fn], "\\s*:\\s*")[2] + 1,
           stop = nchar(ln[ln_fn])) <-
      paste0(filename_new, "_", start, "_till_", end, ".MON")
    # Start time
    ln_st <- grep(pattern = "(Start date / time)", x = ln)
    substr(ln[ln_st],
           start = stringr::str_locate(ln[ln_st], "\\s*=")[2] + 1,
           stop = nchar(ln[ln_st])) <-
      format(ind[1], format = "%S:%M:%H %d/%m/%y")
    # End time
    ln_ed <- grep(pattern = "(End date / time)", x = ln)
    substr(ln[ln_ed],
           start = stringr::str_locate(ln[ln_ed], "\\s*=")[2] + 1,
           stop = nchar(ln[ln_ed])) <-
      format(ind[data_len], format = "%S:%M:%H %d/%m/%y")

    ## write files ##
    # backup original data
    file.copy(filename[i], paste0(filename[i], "_backup"))
    # new file
    writeLines(ln, paste0(dirname(filename[i]), "/", filename_new, "_", start, "_till_", end, ".MON"))

    # remove original file
    fs::file_delete(filename[i])
  }
}
