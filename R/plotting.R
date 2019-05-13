#==============================================================================#
#Plotting ####

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
                           oob = squish)+
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
                           oob = squish)+
      scale_color_gradientn(colors = colorRampPalette(brewer.pal(11,"Spectral"))(11),
                            na.value = "transparent", name = "Lag [h]",
                            breaks = c(-24,-12,-6,-3,-1.5,0, 1.5, 3,6,12,24),
                            limits = c(-24,24),
                            labels = c("< -24",-12,-6,-3, -1.5, 0, 1.5, 3,6,12,">24"),
                            values = c(0,0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6,1),
                            oob = squish)+
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
