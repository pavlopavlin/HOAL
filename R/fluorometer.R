#==============================================================================#
##        Fluorometer functions
#==============================================================================#

#' Reads procesed .ppb file from fluorometer software
#'
#' @param file character. Filepath of .ppb file
#'
#' @return xts object with columns:
#'  \describe{
#'    \item{Tracer1}{Tracer from class 1 (normaly Uranine) [ppb].}
#'    \item{Tracer2}{Tracer from class 2 (normaly Suphurhodamine B) [ppb].}
#'    \item{Tracer3}{Tracer from class 3 (normaly Tinopal) [ppb].}
#'    \item{Turbidity}{Turbidity [ppb]}
#'    \item{Temperature}{[°C]}
#'  }
#' @export
ReadPpb <- function(file){
  df <- read.csv(file = file, header = F, skip = 3, sep ="",
                 stringsAsFactors = F, skipNul =T,
                 col.names = c("ID", "Date", "Tracer1","Tracer2","Tracer3",
                               "Turbidity", "Temperature"),
                 colClasses = c("numeric", "character",rep("numeric",5)))
  df.xts <- xts(df[,3:7], order.by = as.POSIXct(df$Date, tz = "Etc/GMT-1",
                                                format = "%d/%m/%y-%H:%M:%S"))

  return(df.xts)
}

##============================================================================##

#' Reads raw fluorometer file .mv
#'
#' @param file character. Filepath of .mv file
#'
#' @return xts object with columns:
#'  \describe{
#'    \item{Tracer1}{Tracer from class 1 (normaly Uranine) [ppb].}
#'    \item{Tracer2}{Tracer from class 2 (normaly Suphurhodamine B) [ppb].}
#'    \item{Tracer3}{Tracer from class 3 (normaly Tinopal) [ppb].}
#'    \item{Turbidity}{Turbidity [ppb]}
#'    \item{Baseline}{}
#'    \item{BatteryV}{Battery voltage [V]}
#'    \item{Temperature}{[°C]}
#'    \item{conductivity}{Electrical conductivity[]}
#'  }
#' @export
ReadMv <- function(file){
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
