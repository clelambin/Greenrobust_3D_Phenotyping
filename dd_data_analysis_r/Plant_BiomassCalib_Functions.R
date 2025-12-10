# Functions for biomass callibration project

# Functions
plot_obs_vs_pred <- function(mdl, img_name="", data=NA, obs=NA, xlab=NA, ylab=NA)
{
  # Description: Plot observed data in function of prediction for input model
  
  # If observation not specified, use the response variable from the model
  if(all(is.na(obs))){obs <- mdl$model[,1]}
  # If labels not specified, use the names of the response variable in the model
  if(is.na(xlab)){xlab <- paste("Predicted", names(mdl$model)[1])}
  if(is.na(ylab)){ylab <- paste("Measured", names(mdl$model)[1])}
  # If data not specified, take from model
  if(all(is.na(data))){data <- eval(mdl$call$data)}
  
  # Initialise plot and save image (if set)
  par(par_init)
  if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
  do.call(par, c(list(mar=c(4, 4, 1, 1)), par_args))
  # Check available species levels in current dataframe
  avail_species <- levels(data$species)
  
  # Plot prediction
  plot(obs ~ predict(mdl, type="response"),
       col=species_color[data$species], pch=16,
       xlab=xlab, ylab=ylab)
  legend("bottomright", legend=avail_species, col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)
  
  # Add 1:1 line
  range_min <- min(obs)
  range_max <- max(obs)
  lines(c(range_min, range_max), c(range_min, range_max), lty=3)
  # Add R squared (extracted from summary for lm, computed for glm)
  if(is.null(mdl$deviance)){
    r_squared <- round(as.numeric(summary(mdl)["r.squared"]), digit=3)
  } else {
    r_squared <- round(1 - mdl$deviance/mdl$null.deviance, digit=3)
  }
  
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
  # Set test type based on model family
  if(is.null(mdl$family)){
    test="F"
  } else {
    test="Chisq"
  }
  mdl_drop <- drop1(mdl, test=test)
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