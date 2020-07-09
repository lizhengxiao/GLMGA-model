## difine the qq-plot function for data set

qqplot.gg <- function(x, y, titlegg = titlegg, ...){
  newtheme <-   theme_bw() + theme(axis.line = element_line(colour = "black"),
                                   axis.text.x = element_text(margin = margin(t = 0.25, unit = "cm")),
                                   axis.text.y = element_text(margin = margin(r = 0.25, unit = "cm"),
                                                              size = 10, 
                                                              vjust = 0.5, 
                                                              hjust = 0.5),
                                   axis.ticks.length = unit(-0.1, "cm"),
                                   plot.title = element_text(hjust = 0.5),
                                   legend.direction = 'vertical',
                                   legend.position = c(.95, .99),
                                   legend.justification = c("right", "top"),
                                   legend.box.just = "right",
                                   legend.text = element_text(size = 10)) 
  x <- sort(na.omit(x))
  y <- sort(na.omit(y))
  qy <- quantilefun(y)
  m <- length(x)
  n <- length(y)
  N <- m + n
  M <- m * (n/N)
  K <- 1.36
  p <- (1:m - 1)/(m - 1)
  yq <- qy(p)
  yl <- qy(p - K/sqrt(M))
  yu <- qy(p + K/sqrt(M))
  df <- data.frame(
    observed = x,
    expected = y,
    clower   = yl,
    cupper   = yu
  )
  
  p0 <-  ggplot(df, aes(x = observed, y = expected)) + newtheme + 
    geom_point(aes(observed, expected)) +
    geom_abline(intercept = 0, slope = 1, alpha = 1, col = 'red', lwd = 1) +
    labs(title = titlegg, x = "Sample Quantiles", y = "Theoretical Quantiles")
  p0
}