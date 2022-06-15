#' plot PRS distribution
#'
#' Generates simple plots of PRS distribution. For quantitative trait, generates a scatter plot of TRAITxPRS.
#' For binary trait, generates a desnity plot of PRS stratified by status.
#' @importFrom ggplot2 ggplot aes geom_point theme_bw ylab xlab theme geom_density scale_fill_manual geom_smooth
#' @importFrom MASS kde2d
#' @importFrom stats density
#' @importFrom rlang .data
#' @param trait trait vector. Can be binary or quantitative
#' @param PRS PRS vector
#' @param add.density If quantitative trait, color points by densiy? Defaults to TRUE
#' @export

plotPRS <- function(trait, PRS, add.density=TRUE){
  dat <- data.frame(TRAIT=trait, PRS=PRS, stringsAsFactors = FALSE)

  if(length(unique(trait)) > 2){
    if(add.density == TRUE){
      #code from https://slowkow.com/notes/ggplot2-color-by-density/
      get_density <- function(x, y, ...) {
        dens <- kde2d(x, y, ...)
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        return(dens$z[ii])
      }
      dat$density <- get_density(dat$TRAIT, dat$PRS, n=100)

      g <- ggplot(data=dat, aes(.data$TRAIT, .data$PRS, color=density)) +
        geom_point(size=1)+
        theme_bw()
      #g <-  g + theme(legend.position='none')
      g <- g+ xlab("TRAIT") + ylab("PRS") +
        geom_smooth(method='lm', color="red")
    }else{
      g <- ggplot(data=dat, aes(.data$TRAIT, .data$PRS)) +
        geom_point(size=1)+
        theme_bw()
      g <- g+ xlab("TRAIT") + ylab("PRS") +
        geom_smooth(method='lm', color="red")
    }
  }else{
    dat$TRAIT <- as.factor(ifelse(dat$TRAIT == 1, "CASE", ifelse(dat$TRAIT == 0, "CONTROL", dat$TRAIT)))

    g <- ggplot(data=dat, aes(.data$PRS, fill=.data$TRAIT)) +
      geom_density(alpha = 0.5)+
      theme_bw()
    g <- g + scale_fill_manual(values = c("#E69F00", "#0072B2"))
    g <- g+ xlab("TRAIT") + ylab("PRS")
  }

  print(g)
}
