#' @title
#' Plot standard maps
#'

#' @export
plot_legend <-
function( Zlim,
          legend_x = c(0,0.05),
          legend_y = c(0.05,0.45),
          cex.legend = 1,
          col = viridisLite::viridis,
          legend_digits = 1,
          f = identity ){

      # Boundaries
      xl = (1-legend_x[1])*par('usr')[1] + (legend_x[1])*par('usr')[2]
      xr = (1-legend_x[2])*par('usr')[1] + (legend_x[2])*par('usr')[2]
      yb = (1-legend_y[1])*par('usr')[3] + (legend_y[1])*par('usr')[4]
      yt = (1-legend_y[2])*par('usr')[3] + (legend_y[2])*par('usr')[4]

      # Logic
      if( diff(legend_y) > diff(legend_x) ){
        align = c("lt","rb")[2]
        gradient = c("x","y")[2]
      }else{
        align = c("lt","rb")[1]
        gradient = c("x","y")[1]
      }

      # Make plot
      plotrix::color.legend( xl=xl,
                             yb=yb,
                             xr=xr,
                             yt=yt,
                             legend = round(f(seq(Zlim[1],Zlim[2],length=4)),legend_digits),
                             rect.col=col(1000),
                             cex=cex.legend,
                             align=align,
                             gradient=gradient)
}
