# Functions for biomass callibration project

plot_obs_vs_pred <- function(obs, mdl, data=plant_correlation, img_name="", xlab="", ylab="")
{
  # Description: Plot objserved data in function of prediction for input model
  
  # Initialise plot and save image (if set)
  par(par_init)
  if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
  do.call(par, c(list(mar=c(4, 4, 1, 1)), par_args))
  # Check available species levels in current dataframe
  avail_species <- levels(data$species)
  
  # Plot prediction
  plot(obs ~ predict(mdl),
       col=species_color[data$species], pch=16,
       xlab=xlab, ylab=ylab)
  legend("bottomright", legend=avail_species, col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)
  
  # Add 1:1 line
  range_min <- min(obs)
  range_max <- max(obs)
  lines(c(range_min, range_max), c(range_min, range_max), lty=3)
  # Add R squared
  r_squared <- round(as.numeric(summary(mdl)["r.squared"]), digit=3)
  formula   <- Reduce(paste, deparse(mdl$call$formula))
  mtext(paste0(" ", formula), side=3, adj=0, line=-1.25, cex=0.8)
  mtext(paste0(" R²=", r_squared), side=3, adj=0, line=-2, cex=0.8)
  if(save_plot){dev.off()}
}

model_stats <- function(mdl, plot_res=TRUE)
{
  # Description display stats for input model
  if(plot_res)
  {
    # Reset graphic before simulate residuals
    par(par_init)
    simulateResiduals(mdl, plot=TRUE)
  }
  # Save stats
  mdl_sum <- summary(mdl)
  mdl_aic <- AIC(mdl)
  mdl_drop <- drop1(mdl, test="F")
  # Return stats as list
  return(list(summary=mdl_sum, AIC=mdl_aic, drop1=mdl_drop))
}

transformation_mapping <- function(transf)
{
  # Based on input transformation, change some of the display and predict value
  if(transf == "")
  {
    log_axis    <- ""
    back_transf <- function(x) x
  }
  else if(transf == "log")
  {
    log_axis    <- "xy"
    back_transf <- exp
  }
  else if(transf == "sqrt")
  {
    log_axis    <- ""
    back_transf <- function(x) x^2
  }
  else
  {
    warning("Transformation unknown, use no transformation")
    log_axis    <- ""
    back_transf <- function(x) x
  }
  # Return axis and back transformation function
  return(list(log_axis=log_axis, back_transf=back_transf))
}