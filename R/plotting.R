#==============================================================================#
##  Plotting  ##
#==============================================================================#

#==============================================================================#
## dygraphs ####

#' dygraphs plugin - highlights one time series
#'
#' @param dygraph dygraph object.
#'
#' @return dygraph object.
#' @export
#'
dyHighlighted <- function(dygraph){
  dyHighlight(dygraph, highlightSeriesOpts = list(strokeWidth = 3)) %>%
    dyCSS(textConnection("
                         .dygraph-legend > span { display: inline; }
                         .dygraph-legend > span.highlight { border-radius: 5px 5px 5px;
                         border: 1px solid red;}
                         "))%>%
    dyCrosshair(direction = "vertical")%>%
    dyOptions(useDataTimezone = T)
}

##----------------------------------------------------------------------------##

#' Shading effect for dygraph
#'
#' @description Extends \link[dygraphs]{dyShading} to allow multiple areas to be ploted at once.
#'
#' @param dygraph Dygraph to add shading to.
#' @param from Date/time or numeric vector to shade from (for date/time this must be a as.POSIXct object or another object convertible via as.POSIXct).
#' @param to Date/time or numeric vector to shade to (for date/time this must be a as.POSIXct object or another object convertible via as.POSIXct).
#' @param color Color of shading. This can be of the form "#AABBCC" or "rgb(255,100,200)" or "yellow". Defaults to a very light gray.
#' @param axis Axis to apply shading. Choices are "x" or "y".
#'
#' @import dygraphs
#' @note See the
#'   \href{https://rstudio.github.io/dygraphs/gallery-annotations.html}{online
#'   documentation} for additional details and examples.
#' @return A dygraph with the specified shading
#' @export
dyShadings <- function(dygraph, from, to, color = "#EFEFEF", axis = "x"){
  if(length(from) != length(to))
    stop(paste0("'from' and 'to' have to be of the same length!
                \n length(from) = ",length(from), ", length(to) = ", length(to)))
  if(length(from) < 1) stop("'from' and 'to' must be at least of length 1!")

  for(ii in 1:length(from)){
    dygraph <- dygraph %>%
      dyShading(from = from[ii], to = to[ii],
                color = color, axis = axis)
  }
  dygraph
}

##----------------------------------------------------------------------------##

#' Plot of lag times from ccf as ellipses
#'
#' Plot of lag times between all stations in the data in the form of collored ellipses
#'
#' @param lag.df melted data frame of lag times.
#' @param title optional plot main title.
#' @param subtitle An optional plot subtitle
#' @param file An optional file name to which the plot should be saved.
#' @param type The type of the plot. Defaults to \code{"ellipse"}.
#' @param invert should the values in the lower triangle be inverted?
#'
#' @details
#' \code{lag.df} should only have data for lower triangel (below diagonal) and
#' have the following columns:
#' \describe{
#'   \item{Var1}{Name of the station on the x axis.}
#'   \item{Var2}{Name of the station on the y axis.}
#'   \item{value}{Numeric value of the lag time in hours}
#' }
#' There are two types of this plot: \code{"ellipse"} plots the whole matrix as
#' ellipses, where \code{"mixed"} plots numerical values above the diagonal
#' instead.
#'
#' The plotted lag is positive when column station laggs behind the row station.
#'
#' @return ggplot object
#'
#' @examples \dontrun{
#' LagCorrPlot(ccf.lag.sig.melt,
#'             title = "Mean lag time of the column stations to the row stations",
#'             subtitle = "CCF on differentiated 1h time-series, CCF > 0.4",
#'             file = paste0(dir.plot, "/ccf_lag_sig0.4_diff_1h_ellipse.png"))
#' }
#'
#' @export
#' @import ggplot2
LagCorrPlot <- function(lag.df, title = "Mean lag time of the column stations to the row stations",
                        subtitle = "CCF on differentiated 1h time-series, CCF > 0.4", file = NULL,
                        type = c("ellipse", "mixed"), invert = T) {

  ell.dat <- function(rho, length = 99, max.rho = 25) {
    rho[rho >= max.rho] <- max.rho
    rho[rho <= -max.rho] <- -max.rho
    rho <- rho / (max.rho + .01*max.rho)
    k <- seq(0, 2 * pi, length = length)
    if(is.na(rho)){
      x <- rep(0,length(k))
      y <- rep(0,length(k))
    }else{
      x <- cos(k + acos(rho)/2)/2
      y <- cos(k - acos(rho)/2)/2
    }
    cbind(rbind(x, y), c(NA, NA))
  }

  type = match.arg(type)

  dat <- cbind(lag.df,
               ID =1:nrow(lag.df),
               x = as.numeric(lag.df$Var1),
               y = as.numeric(lag.df$Var2))

  ELL.dat <- lapply(dat$value, ell.dat, max.rho=24)
  ELL.dat2 <- 0.85 * matrix(unlist(ELL.dat), ncol = 2,
                            byrow = TRUE)
  ELL.dat2 <- ELL.dat2 + dat[rep(1:nrow(dat), each = 100), 5:6]
  ELL.dat2 <- cbind(ID = dat[rep(1:nrow(dat), each = 100), 4], ELL.dat2)
  dat.plot <- merge(dat[1:4], ELL.dat2, by = "ID")
  dat.plot$Discrete = cut(dat.plot$value, c(-24,-12,-6,-3,-1.5,-0.1,0.1, 1.5, 3,6,12,24),include.lowest=T)


  ## ellipse
  if(type == "ellipse"){
    g1 <- ggplot()+
      geom_abline(slope = 1, intercept = 0)+
      geom_polygon(data= dat.plot, aes(x=x, y=y, group = ID, fill = value))
    if (invert == T){
      g1 <- g1 + geom_polygon(data= dat.plot,
                              aes(x = 2*as.numeric(Var2)-y, y=x,
                                  group = ID, fill = -value))
    }else{
      g1 <- g1 + geom_polygon(data= dat.plot,
                              aes(x = 2*as.numeric(Var2)-y, y=x,
                                  group = ID, fill = value))
    }
    g1 <- g1 +
      theme(aspect.ratio = 1,
            axis.text.x.bottom =  element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.x.top =  element_text(angle = 90, hjust = 0, vjust = 0.5))+
      scale_fill_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral"))(11),
                           na.value = "transparent", name = "Lag [h]",
                           breaks = c(-24,-12,-6,-3,-1.5,0, 1.5, 3,6,12,24),
                           limits = c(-24,24),
                           labels = c("<24",-12,-6,-3, -1.5, 0, 1.5, 3,6,12,">24"),
                           values = c(0,0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6,1),
                           oob = scales::squish)+
      scale_x_continuous(expand = expand_scale(add=0.5),
                         breaks = 1:nrow(ccf.all.mean),
                         minor_breaks = NULL,
                         labels = rownames(ccf.all.mean),
                         name = NULL,
                         sec.axis = dup_axis())  +
      scale_y_continuous(expand = expand_scale(add=0.5),
                         breaks = 1:nrow(ccf.all.mean),
                         minor_breaks = NULL,
                         labels = rownames(ccf.all.mean),
                         name = NULL,
                         sec.axis = dup_axis())+
      guides(fill= guide_colorbar(barheight=20)) +
      ggtitle(label = title, subtitle = subtitle)


    ## ellipse + text
  }else if(type == "mixed"){
    g1 <- ggplot()+
      geom_abline(slope = 1, intercept = 0)+
      geom_polygon(data= dat.plot, aes(x=x, y=y, group = ID, fill = value)) +
      {if(invert ==T){
        geom_text(data = dat.plot,
                  aes(x=as.numeric(Var2), y = as.numeric(Var1),
                      label = round(-value,1)),
                  show.legend = F, size = 2, color = "black")
      }else{
        geom_text(data = dat.plot,
                  aes(x=as.numeric(Var2), y = as.numeric(Var1),
                      label = round(value,1)),
                  show.legend = F, size = 2, color = "black")
      }
      }+
      theme(aspect.ratio = 1,
            axis.text.x.bottom =  element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.text.x.top =  element_text(angle = 90, hjust = 0, vjust = 0.5))+
      scale_fill_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral"))(11),
                           na.value = "transparent", name = "Lag [h]",
                           breaks = c(-24,-12,-6,-3,-1.5,0, 1.5, 3,6,12,24),
                           limits = c(-24,24),
                           labels = c("<-24",-12,-6,-3, -1.5, 0, 1.5, 3,6,12,">24"),
                           values = c(0,0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6,1),
                           oob = scales::squish)+
      scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral"))(11),
                            na.value = "transparent", name = "Lag [h]",
                            breaks = c(-24,-12,-6,-3,-1.5,0, 1.5, 3,6,12,24),
                            limits = c(-24,24),
                            labels = c("< -24",-12,-6,-3, -1.5, 0, 1.5, 3,6,12,">24"),
                            values = c(0,0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6,1),
                            oob = scales::squish)+
      scale_x_continuous(expand = expand_scale(add=0.5),
                         breaks = 1:nrow(ccf.all.mean),
                         minor_breaks = NULL,
                         labels = rownames(ccf.all.mean),
                         name = NULL,
                         sec.axis = dup_axis())  +
      scale_y_continuous(expand = expand_scale(add=0.5),
                         breaks = 1:nrow(ccf.all.mean),
                         minor_breaks = NULL,
                         labels = rownames(ccf.all.mean),
                         name = NULL,
                         sec.axis = dup_axis())+
      guides(fill= guide_colorbar(barheight=20)) +
      ggtitle(label = title, subtitle = subtitle)
  }

  if(is.character(file)) {
    ggsave(file ,plot = g1, device = "png",width = 7, height = 7, dpi = 1000)
  }
  g1
} #end LagCorrPlot

##----------------------------------------------------------------------------##
#' Save dygraph to file
#'
#' @param dygraph dygraph object
#' @param filename File name to create on disk. Without extension.
#' @param device Device to use. Supported: "png", "jpg","jpeg", "pdf"
#' @param width,height image width/height.
#' @param units character. Units of height and width.
#' @param dpi integer. Image resolution in DPI.
#' @param zoom numeric. Functions as zoom in a web browser.
#'    Higher number increases the labels' size.
#' @import fs
#' @export
#'
#' @import dygraphs
#'
#' @examples
#' \dontrun{
#'   ## Save a plot of yearly sunspot data
#'   x <- xts::as.xts(sunspot.year)
#'   d1 <- dygraphs::dygraph(x)
#'   dySave(d1, "sunspots_yearly", device = "jpg", width = 8, height = 6)
#'   }
#'

dySave <- function(dygraph, filename, device = "png", width = NA, height = NA,
                   units = c("in", "cm", "mm"), dpi = 300, zoom = 3.75){
  # get extension and set it as 'device'
  ext <- fs::path_ext(filename)
  if(ext != "") device <- ext

  if(!device %in% c("png", "jpg","jpeg", "pdf", "html")){
    stop(paste("unsupported driver:", device))
  }

  # get absolute path without extension
  ifelse(fs::is_absolute_path(filename),
         file <- fs::path_ext_remove(filename),
         file <- fs::path_wd(fs::path_ext_remove(filename)))

  on.exit({if(dir.exists(paste0(file, "_lib/")))  try(fs::dir_delete(paste0(file, "_lib/")))})


  htmlwidgets::saveWidget(widget = dygraph,
                          file = paste0(file, ".html"),
                          libdir = paste0(fs::path_file(file), "_lib/"),
                          selfcontained = T)

  if(device != "html"){
  dim <- plot_dim(dim = c(width,height), scale = 1, units = units)
  dim_pix <- floor(dim * dpi)

  webshot::webshot(url = paste0(file, ".html"),
                   file = paste0(file, ".", device),
                   vwidth = dim_pix[1],
                   vheight = dim_pix[2],
                   zoom = zoom)
  }
}


#' Gets plot dimensions for dySave
#'
#' @param dim character vector of width and height.
#' @param scale numeric. Scaling factor.
#' @param units character. Units of input dimensions.
#' @param limitsize logical. if TRUE (default) the size of output image is limited to 50'.
#'
#' @return numeric vector of length 2: width, height in inches.
#'
plot_dim <- function (dim = c(NA, NA), scale = 1, units = c("in", "cm", "mm"), limitsize = T) {
  units <- match.arg(units)
  to_inches <- function(x) x/c(`in` = 1, cm = 2.54, mm = 2.54 *
                                 10)[units]
  from_inches <- function(x) x * c(`in` = 1, cm = 2.54, mm = 2.54 *
                                     10)[units]
  dim <- to_inches(dim) * scale
  if (any(is.na(dim))) {
    default_dim <- c(21.4,10.7)
    dim[is.na(dim)] <- default_dim[is.na(dim)]
    dim_f <- prettyNum(from_inches(dim), digits = 3)
    message("Saving ", dim_f[1], " x ", dim_f[2], " ", units,
            " image")
  }

  if (limitsize && any(dim >= 50)) {
    stop("Dimensions exceed 50 inches (height and width are specified in '",
         units, "' not pixels). If you're sure you want a plot that big, use ",
         "`limitsize = FALSE`.", call. = FALSE)
  }

  dim
}

#==============================================================================#

#' Extract limits from ggplot
#'
#' @param plot ggplot object
#'
#' @return matrix with four rows (xmin, xmax, ymin, ymax) and number of columns
#'   equal to the number of panels (facets)
#' @export
#'
#' @examples
#'   \dontrun{
#'   library(ggplot)
#'   g1 <-
#'   ggplot(mpg, aes(hwy,cty)) +
#'     geom_point() +
#'     facet_wrap(drv~., scales = "free_y")
#'   gglimits(g1)
#'   }
gglimits <- function(plot) {
  gb = ggplot_build(plot)
  sapply(gb$layout$panel_params, function(ls){
    c(xmin = ls$x.range[1],
      xmax = ls$x.range[2],
      ymin = ls$y.range[1],
      ymax = ls$y.range[2])
  })
}
