# ---- Header ----
# Name: Plant_BiomassCalib
# Description: Calibrate plant biomass based on data extracted from 3D model

# ---- Preparation ----

#setting working directory & clean up environment
rm(list=ls())
# Reset display and save initial par
if(length(dev.list()!=0)){dev.off()}
par_init <- par(no.readonly = TRUE)

# Script variables
input_from3d  <- "data//Plant_data_20250717_All_Pot13cm.csv"
input_biomass <- "data//20250623_biomass_data.csv"
output_label  <- "20250717"
output_folder <- "output"
save_plot     <- TRUE
# Set default arguments for plot display (par) and image export (jpeg)
jpeg_args     <- list(height=4, width=6, units="in", res=300)
# mgp=c(title, labels, line) set the distance for the axis (default is (3, 2, 0))
par_args      <- list(cex = 1, cex.axis=0.7, mgp=c(2, 0.8, 0))

# Libraries
library(stringr)      # Regular expression
library(DHARMa)       # Check model assumption
library(Hmisc)        # Cross-correlation with significance test
library(corrplot)     # Correlation plot

# ---- User function ----
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
  # Return axis and back transformation funtion
  return(list(log_axis=log_axis, back_transf=back_transf))
}

plot_predict <- function(dataframe, img_name="", main="", print_model=FALSE, predictor="temperature", transf="log")
{
  # Description: Plot biomass in function of volume and fit linear model based on input predictor
  
  # Map axis and back transformation based on transformation value
  transf_map  <- transformation_mapping(transf)
  log_axis    <- transf_map$log_axis
  back_transf <- transf_map$back_transf
  
  # Check number of level of input predictor
  predictor_level = levels(dataframe[,predictor])
  number_of_level = length(predictor_level)
  # Plot biomass in function of volume
  if(img_name != ""){do.call(jpeg, c(filename=img_name, jpeg_args))}
  do.call(par, c(list(mar=c(3.5, 3.5, 3, 1)), par_args))
  color_palette = colorRampPalette(c("blue", "yellow", "red"))
  color_predictor = color_palette(number_of_level)
  plot(biomass~Volume,
       data=dataframe,
       col=color_predictor[dataframe[,predictor]],
       xlab="Volume (m3)", ylab="Biomass(g)",
       main=main, cex.lab=1, pch=16, log=log_axis)
  
  # Fit linear model
  lm_biomass <- lm(paste0("log(biomass)~log(Volume)+", predictor), data=dataframe)
  # simulateResiduals(lm_biomass, plot=TRUE)
  if(print_model)
  {
    print(summary(lm_biomass))
    print(drop1(lm_biomass, test="F"))
  }
  
  # Compute prediction for each predictor level
  for(i in 1:length(predictor_level))
  {
    predictor_value <- predictor_level[i]
    volume_for_predict <- dataframe[dataframe[,predictor]==predictor_value, "Volume"]
    pred_biomass <- data.frame(Volume=seq(from=min(volume_for_predict),
                                         to=max(volume_for_predict),
                                         length.out=20))
    pred_biomass[,predictor] <- rep(predictor_value, 20)
    pred_biomass$biomass <- back_transf(predict(lm_biomass, newdata = pred_biomass))
    lines(pred_biomass$Volume, pred_biomass$biomass, col=color_predictor[i], lty=2)
  }
  # Add legend and save image
  legend(x="bottomright", legend=levels(dataframe[,predictor]), col=color_predictor, cex=0.8, pch=16, title=predictor)
  if(img_name != ""){dev.off()}
}

plot_predict_multivar <- function(dataframe, display_var, linear_model,
                                  static_vars= c("Volume", "Dim_Z", "Top_Area"),
                                  predictor="Species", plot_legend=FALSE, plot_conf=FALSE,
                                  legend_x=0.2, legend_y=35, transf="log")
{
  # Description: plot biomass in function of display var and fit linear model (keeping constant non input variable)
  
  # Map axis and back transformation based on transformation value
  transf_map  <- transformation_mapping(transf)
  log_axis    <- transf_map$log_axis
  back_transf <- transf_map$back_transf
  
  # Color point based on species
  color_palette <- colorRampPalette(c("blue", "yellow", "red"))
  predictor_level <- levels(dataframe[,predictor])
  number_of_level <-  length(predictor_level)
  color_predictor <- color_palette(number_of_level)
  
  # Plot base points
  plot(biomass~ dataframe[,display_var], data=dataframe, log=log_axis, pch=16,
       col=color_predictor[dataframe[,predictor]], xlab="", ylab="")
  
  # Compute prediction for input linear model for each predictor levels
  for(i in 1:number_of_level)
  {
    # Create subset dataframe with only values of given predictor
    predictor_value <- predictor_level[i]
    var_for_predict <- dataframe[dataframe[,predictor]==predictor_value, display_var]
    # Initialise predict dataframe, require at least one variable with fixed name
    pred_biomass <- data.frame(dummy_var = rep(NA, 20))
    pred_biomass[,predictor] <- rep(predictor_value, 20)
    for(static_var in static_vars)
    {
      pred_biomass[,static_var] = rep(mean(dataframe[,static_var]), 20)
    }
    pred_biomass[,display_var] <- seq(from=min(var_for_predict),
                                      to=max(var_for_predict),
                                      length.out=20)
    # Compute predictor and overlay plot with prediction line
    lm_predict <- predict(linear_model, newdata = pred_biomass, se.fit=TRUE, interval="confidence", level=0.95)
    lm_predict <- back_transf(lm_predict$fit)
    # If plot_conf, plot both sup and inf confidence interval line, otherwise, only the prediction line
    nb_line <- ifelse(plot_conf, ncol(lm_predict), 1)
    for(j in 1:nb_line)
    {
      # Set different line type and width for prediction line and confidence interval
      line_type  <- ifelse(j==1, 1, 2)
      line_width <- ifelse(j==1, 2, 1)
      lines(pred_biomass[,display_var], lm_predict[,j], col=color_predictor[i], lty=line_type, lwd=line_width)
    }
  }
  if(plot_legend)
  {
    legend(legend_x, legend_y, legend=species_level, col=color_predictor, pch=16, bty="n", horiz=TRUE, xjust=0.5, cex=0.8)
  }
}

# ---- Data preparation ----
## ==== Data extraction ====

# Load data
plant_biomass <- read.csv(input_biomass)
plant_volume  <- read.csv(input_from3d)

# Extract label from model name
plant_volume$label <- str_extract(plant_volume$Plant.name, "[A-Z]{2}[0-9]{3}_*[A-Z][0-9]{2}")
plant_volume$label <- str_remove_all(plant_volume$label, "_")

# Add species as variable
species_label <- str_extract(plant_volume$label, "^[A-Z]{2}")
species_corresp <- c(BR= "Brassica rapa",
                     HS="Hordeum spontaneum",
                     HV="Hordeum vulgare",
                     NB="Nicotiana benthamiana",
                     SD="Solanum dulcamara",
                     SL="Solanum lycopersicum",
                     TA="Thlaspi arvense")
#plant_correlation$Species <- species_corresp[species_label]
plant_volume$Species <- species_label
plant_volume$Species <- as.factor(plant_volume$Species)

# Look for matching label
plant_correlation <- plant_volume
#plant_correlation$biomass <- plant_biomass[plant_biomass$pot_id == plant_correlation$label,"Biomass_g"]
plant_correlation$biomass <- NA
for(i in 1:nrow(plant_correlation))
{
  single_label <- plant_correlation[i, "label"]
  single_biomass <- plant_biomass[plant_biomass$pot_id == single_label, "Biomass_g"]
  print(single_label)
  print(single_biomass)
  if(length(single_biomass > 0))
  {
    plant_correlation[i, "biomass"] <- single_biomass
  }
}

# Save temperature as separate treatment
plant_correlation$temperature <- substr(str_extract(plant_correlation$label, "[0-9]S[0-9]{2}"), 2,4)
plant_correlation$temperature <- factor(plant_correlation$temperature)

## ==== Data correction ====
# In 3D extraction, scale all pot to 13cm side, need to correct for plant with actual pot size of 11cm
species_11cm <- c("BR", "HS", "HV", "NB", "TA")
corrector    <- 0.11 / 0.13
# Apply correction to length (with factor 1), area (factor 2) and volume (factor 3)
to_correct   <- c("Dim_X", "Dim_Y", "Dim_Z", "Cumul_Area", "Volume")
corr_factor  <- c(1,       1,       1,       2,            3)
# Apply corrector to measures for given species
corr_species <- plant_correlation$Species %in% species_11cm
for(i in 1:length(to_correct))
{
  plant_correlation[corr_species, to_correct[i]] <- plant_correlation[corr_species, to_correct[i]] * corrector^corr_factor[i]
}

## ==== Projected area ====
# Compute projected area from X, Y and Z Dim and area ratio
# BBox area is total area of the bounding box
plant_correlation$Front_BBox <- plant_correlation$Dim_Y * plant_correlation$Dim_Z
plant_correlation$Left_BBox  <- plant_correlation$Dim_X * plant_correlation$Dim_Z
plant_correlation$Top_BBox   <- plant_correlation$Dim_X * plant_correlation$Dim_Y
# Projected area is BBox * Area ratio
plant_correlation$Front_Area <- plant_correlation$Front_BBox * plant_correlation$Font_Area_Ratio
plant_correlation$Left_Area  <- plant_correlation$Left_BBox  * plant_correlation$Left_Area_Ratio
plant_correlation$Top_Area   <- plant_correlation$Top_BBox   * plant_correlation$Top_Area_Ratio

## ==== Data cleanup ====

# Remove outliers
plant_correlation <- plant_correlation[plant_correlation$Volume<0.01,]
plant_correlation <- plant_correlation[plant_correlation$Volume!=0,]

# Remove all plants with a code error
plant_correlation <- plant_correlation[plant_correlation$Code_Error == 0,]

# Remove data without biomass info
plant_correlation <- plant_correlation[!is.na(plant_correlation$biomass),]

# Cleanup unused factors
plant_correlation$Species <- droplevels(plant_correlation$Species)

## ==== Confidence assessment ====
plant_confidence <- data.frame(Model.issue = plant_correlation$Model.issue)
# Define confidence rate, from 1 to 4 based on model issue
# 1 : No issue
# 2 : "Slight"
# 3 : "Bad"
# 4 : "Very bad"
# Set default confidence to 1
plant_confidence$confidence = 1
# If above pattern present, update confidence
confidence_pattern = c("Slight", "Bad", "Very bad")
for(conf in 1:length(confidence_pattern))
{
  pattern = confidence_pattern[conf]
  has_patern = grepl(pattern, plant_confidence$Model.issue, fixed = TRUE)
  plant_confidence$confidence[has_patern] = conf + 1
}
# Use confidence to filter plant correlation dataframe
plant_correlation_no_issue = plant_correlation[plant_confidence$confidence==1,]
plant_correlation_no_bad   = plant_correlation[plant_confidence$confidence <3,]

## ==== Overview ====
# Processed plant overview (per species)
table_init <- table(plant_volume$Species)
table_corr <- table(plant_correlation$Species)
table_nobad   <- table(plant_correlation_no_bad$Species)
table_noissue <- table(plant_correlation_no_issue$Species)

table_combine <- rbind(c(table_init, sum(table_init)),
                       c(table_corr, sum(table_corr)),
                       c(table_nobad, sum(table_nobad)),
                       c(table_noissue, sum(table_noissue)))
colnames(table_combine)[ncol(table_combine)] <- "sum"
rownames(table_combine) <- c("Processed data", "No reported error", "No extra bad error", "No extra error")
table_combine

# Check final prepared data
summary(plant_correlation)
str(plant_correlation)

# ---- Data analysis ----
## ==== Cross-correlation ====
# From: https://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# Note: Pearson vs Spearman (https://www.reddit.com/r/statistics/comments/76iw0w/how_do_you_decide_which_correlation_coefficient/)
#       Pearson measures the linearity of data
#       Spearman is computed on ranks and measures the monotonicity of data (doesn't necessarily have to be linear)
correlation_data <- as.matrix(subset(plant_correlation, select=-c(Path, Plant.name, Code_Error, Model.issue, Cause, label, Species, temperature)))
summary(correlation_data)
correlation_mat <- rcorr(correlation_data, type="pearson")
# (Set diagonal of significance to 0, not computed during correlation)
diag(correlation_mat$P) <- 0
# Plot correlation table
img_name <- paste0(output_folder, "//Correlation_matrix_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
# Reset graphic
par(par_init)
corrplot(correlation_mat$r, type = "upper", order = "hclust", 
         p.mat = correlation_mat$P, sig.level = 0.01, insig = "blank",
         tl.col = "Black")
if(save_plot){dev.off()}

## ==== Variable distribution ====
# Relevant variable to study distribution from
variable_distr <- c("biomass", "Volume", "Dim_Z", "Top_Area")
variable_label <- c("Biomass (g)", "Volume (m3)", "Plant height (m)", "Top area (m2)")

# Check variable distribution to check if data transformation should be used
img_name <- paste0(output_folder, "//VariableDistribution_noTransf_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(2, 2), mar=c(3.5, 1.5, 0.5, 0.5), oma=c(0, 2.5, 0, 0), xpd=NA), par_args))
for(i in 1:length(variable_distr))
{
  hist(plant_correlation[,variable_distr[i]],
       xlab=variable_label[i], ylab="", main="")
}
mtext("Frequency", side=2, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Check variable distribution with log transformation
img_name <- paste0(output_folder, "//VariableDistribution_withLog_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(2, 2), mar=c(3.5, 1.5, 0.5, 0.5), oma=c(0, 2.5, 0, 0), xpd=NA), par_args))
for(i in 1:length(variable_distr))
{
  hist(log(plant_correlation[,variable_distr[i]]),
       xlab=paste0("Log of ", variable_label[i]), ylab="", main="")
}
mtext("Frequency", side=2, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Check variable distribution with sqrt transformation
img_name <- paste0(output_folder, "//VariableDistribution_withSqrt_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(2, 2), mar=c(3.5, 1.5, 0.5, 0.5), oma=c(0, 2.5, 0, 0), xpd=NA), par_args))
for(i in 1:length(variable_distr))
{
  hist(sqrt(plant_correlation[,variable_distr[i]]),
       xlab=paste0("Sqrt of ", variable_label[i]), ylab="", main="")
}
mtext("Frequency", side=2, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

## ==== (log) Biomass ~ Volume + Temperature ====
# Reset graphic
par(par_init)
# Create prediction for the different filter level to see impact of outliers
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Temperature_LogAxis_", output_label, "_AllPoints.jpg"), "")
plot_predict(plant_correlation, img_name=img_name, print_model=TRUE,
             main="Temperature treatment - All datapoints")
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Temperature_LogAxis_", output_label, "_NoBadIssue.jpg"), "")
plot_predict(plant_correlation_no_bad, img_name=img_name, print_model=TRUE,
             main="Temperature treatment - No bad datapoints")
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Temperature_LogAxis_", output_label, "_NoModelIssue.jpg"), "")
plot_predict(plant_correlation_no_issue, img_name=img_name, print_model=TRUE,
             main="Temperature treatment - No issue datapoint")

## ==== (log) Biomass ~ Volume + Species ====
# Create prediction for the different filter level to see impact of outliers
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Species_LogAxis_", output_label, "_AllPoints.jpg"), "")
plot_predict(plant_correlation, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - All datapoints")
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Species_LogAxis_", output_label, "_NoBadIssue.jpg"), "")
plot_predict(plant_correlation_no_bad, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - No bad datapoints")
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Species_LogAxis_", output_label, "_NoModelIssue.jpg"), "")
plot_predict(plant_correlation_no_issue, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - No issue datapoint")

## ==== (log) Biomass ~ Volume + Dim_Z + Top_Area + Species ====
lm_biomass_volume_top <- lm(log(biomass) ~ log(Volume) + log(Dim_Z) + log(Top_Area) + Species, data=plant_correlation)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_volume_top, plot=TRUE)
summary(lm_biomass_volume_top)
drop1(lm_biomass_volume_top, test="F")

# Set Species color
color_palette <- colorRampPalette(c("blue", "yellow", "red"))
species_level <- levels(plant_correlation$Species)
species_color <- color_palette(length(species_level))

# Plot Volume, Dim_Z and Top_Area in function of each other
img_name <- paste0(output_folder, "//DependantVar_Plot_LogAxis_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(2, 2), mar=c(1.5, 1.5, 1, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot(Volume~Dim_Z, data=plant_correlation, log="xy", pch=16, col=species_color, xlab="", ylab="")
plot(Volume~Top_Area, data=plant_correlation, log="xy", pch=16, col=species_color, xlab="", ylab="")
plot.new()
legend("right", legend=species_level, col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)
plot(Dim_Z~Top_Area, data=plant_correlation, log="xy", pch=16, col=species_color, xlab="", ylab="")
mtext("Volume (m3)", side=2, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=2, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Plot Biomass in function of Volume, Dim_Z and Top_Area
img_name <- paste0(output_folder, "//Biomass_Volume_DimZ_TopArea_LogAxis_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(1, 3), mar=c(1.5, 1.5, 1.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot_predict_multivar(plant_correlation, "Volume", lm_biomass_volume_top)
plot_predict_multivar(plant_correlation, "Dim_Z", lm_biomass_volume_top, plot_legend=TRUE)
plot_predict_multivar(plant_correlation, "Top_Area", lm_biomass_volume_top)
mtext("Biomas (g)", side=2, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Volume (m3)", side=1, adj=0.1, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.9, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

## ==== (log) Biomass ~ Dim_Z + Top_Area + Species ====
lm_biomass_top <- lm(log(biomass) ~ log(Dim_Z) + log(Top_Area) + Species, data=plant_correlation)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_top, plot=TRUE)
summary(lm_biomass_top)
drop1(lm_biomass_top, test="F")

# Observed vs predicted plot
# Reset graphic before simulate residuals
par(par_init)
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_LinearModel_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mar=c(4, 4, 1, 1)), par_args))
plot(log(plant_correlation$biomass) ~ predict(lm_biomass_top),
     col=species_color[plant_correlation$Species], pch=16,
     xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
legend("bottomright", legend=species_level, col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)
lines(c(-2, 3), c(-2, 3), lty=3)
text(-1.8, 2.5, "R²=0.89", cex=0.8)
if(save_plot){dev.off()}

# Observed vs predicted plot (all species)
# Reset graphic before simulate residuals
par(par_init)
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_LinearModel_AllSpecies_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mar=c(4, 4, 1, 1)), par_args))
plot(log(plant_correlation$biomass) ~ predict(lm_biomass_top), pch=20, col="darkgrey",
     xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
lines(c(-2, 3), c(-2, 3), lty=3)
text(-1.8, 2.5, "R²=0.89", cex=0.8)
if(save_plot){dev.off()}

# Plot Biomass in function of Dim_Z and Top_Area only
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_LogAxis_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(1, 2), mar=c(1.5, 1.5, 1.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot_predict_multivar(plant_correlation, "Dim_Z", lm_biomass_top, plot_legend=TRUE, legend_x = 1)
plot_predict_multivar(plant_correlation, "Top_Area", lm_biomass_top)
mtext("Biomas (g)", side=2, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

## ==== (log with interact) Biomass ~ Dim_Z + Top_Area + Species ====
# Test with interaction
lm_biomass_top_with_interact <- lm(log(biomass) ~ (log(Dim_Z) + log(Top_Area)) * Species, data=plant_correlation)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_top_with_interact, plot=TRUE)
summary(lm_biomass_top_with_interact)
drop1(lm_biomass_top_with_interact, test="F")
plot(log(plant_correlation$biomass) ~ predict(lm_biomass_top_with_interact))

# To check: visreg
# Partial residual plot

## ==== (sqrt) Biomass ~ Volume + Species ====
# Reset graphic before simulate residuals
par(par_init)
# Create prediction for the different filter level to see impact of outliers
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Species_Sqrt_", output_label, "_AllPoints.jpg"), "")
plot_predict(plant_correlation, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - All datapoints", transf="sqrt")
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Species_Sqrt_", output_label, "_NoBadIssue.jpg"), "")
plot_predict(plant_correlation_no_bad, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - No bad datapoints", transf="sqrt")
img_name = ifelse(save_plot, paste0(output_folder, "//Biomass_Volume_Species_Sqrt_", output_label, "_NoModelIssue.jpg"), "")
plot_predict(plant_correlation_no_issue, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - No issue datapoint", transf="sqrt")

## ==== (sqrt) Biomass ~ Dim_Z + Top_Area + Species ====
lm_biomass_top_sqrt <- lm(sqrt(biomass) ~ sqrt(Dim_Z) + sqrt(Top_Area) + Species, data=plant_correlation)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_top_sqrt, plot=TRUE)
summary(lm_biomass_top_sqrt)
drop1(lm_biomass_top_sqrt, test="F")

# Plot Biomass in function of Dim_Z and Top_Area only
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_Sqrt_", output_label, "_AllPoints.jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(1, 2), mar=c(1.5, 1.5, 1.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot_predict_multivar(plant_correlation, "Dim_Z", lm_biomass_top_sqrt, plot_legend=TRUE,
                      legend_x=1.5, legend_y=18, transf="sqrt")
plot_predict_multivar(plant_correlation, "Top_Area", lm_biomass_top_sqrt, transf="sqrt")
mtext("Biomas (g)", side=2, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Prediction excluding bad undetected errors
lm_biomass_top_sqrt <- lm(sqrt(biomass) ~ sqrt(Dim_Z) + sqrt(Top_Area) + Species, data=plant_correlation_no_bad)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_top_sqrt, plot=TRUE)
summary(lm_biomass_top_sqrt)
drop1(lm_biomass_top_sqrt, test="F")

# Plot Biomass in function of Dim_Z and Top_Area only
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_Sqrt_", output_label, "_NoBadIssue.jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(1, 2), mar=c(1.5, 1.5, 1.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot_predict_multivar(plant_correlation_no_bad, "Dim_Z", lm_biomass_top_sqrt, plot_legend=TRUE,
                      legend_x=1.5, legend_y=18, transf="sqrt")
plot_predict_multivar(plant_correlation_no_bad, "Top_Area", lm_biomass_top_sqrt, transf="sqrt")
mtext("Biomas (g)", side=2, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Prediction excluding all undetected errors
lm_biomass_top_sqrt <- lm(sqrt(biomass) ~ sqrt(Dim_Z) + sqrt(Top_Area) + Species, data=plant_correlation_no_issue)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_top_sqrt, plot=TRUE)
summary(lm_biomass_top_sqrt)
drop1(lm_biomass_top_sqrt, test="F")

# Plot Biomass in function of Dim_Z and Top_Area only
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_Sqrt_", output_label, "_NoModelIssue.jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(1, 2), mar=c(1.5, 1.5, 1.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot_predict_multivar(plant_correlation_no_issue, "Dim_Z", lm_biomass_top_sqrt, plot_legend=TRUE,
                      legend_x=1.5, legend_y=18, transf="sqrt")
plot_predict_multivar(plant_correlation_no_issue, "Top_Area", lm_biomass_top_sqrt, transf="sqrt")
mtext("Biomas (g)", side=2, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Repeat by removing outlier
plant_correlation_filt_out <- plant_correlation_no_issue[plant_correlation_no_issue$biomass < 15,]
lm_biomass_top_sqrt <- lm(sqrt(biomass) ~ sqrt(Dim_Z) + sqrt(Top_Area) + Species, data=plant_correlation_filt_out)
# Reset graphic before simulate residuals
par(par_init)
simulateResiduals(lm_biomass_top_sqrt, plot=TRUE)
summary(lm_biomass_top_sqrt)
drop1(lm_biomass_top_sqrt, test="F")

# Plot Biomass in function of Dim_Z and Top_Area only
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_Sqrt_", output_label, "_RmOutlier.jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
do.call(par, c(list(mfrow=c(1, 2), mar=c(1.5, 1.5, 1.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA), par_args))
plot_predict_multivar(plant_correlation_filt_out, "Dim_Z", lm_biomass_top_sqrt, plot_legend=TRUE,
                      legend_x=0.9, legend_y=11.5, transf="sqrt")
plot_predict_multivar(plant_correlation_filt_out, "Top_Area", lm_biomass_top_sqrt, transf="sqrt")
mtext("Biomas (g)", side=2, adj=0.5, line=1, cex=1, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1, col="black", outer=TRUE)
if(save_plot){dev.off()}

# ---- Reset ----
# Reset graphic display
par(par_init)
