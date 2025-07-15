# ---- Header ----
# Name: Plant_BiomassCalib
# Description: Calibrate plant biomass based on data extracted from 3D model

# ---- Preparation ----

#setting working directory & clean up environment
rm(list=ls())
setwd("C:/Users/cleme/Documents/Hohenheim/00_Courses/320_Landscape_and_Plant_Ecology/MSc_3D_Plant_Characterisation/Analysis")
# Save initial par
par_init <- par()

# Script variables
input_from3d  <- "Plant_data_20250717_All_Pot13cm.csv"
input_biomass <- "20250623_biomass_data.csv"
output_label  <- "20250717"
save_plot     <- FALSE
plot_height   <- 480
plot_width    <- 640

# Libraries
library(stringr)      # Regular expression
library(DHARMa)       # Check model assumption
library(Hmisc)        # Cross-correlation with significance test
library(corrplot)     # Correlation plot

# ---- User function ----
plot_predict <- function(dataframe, img_name="", main="", print_model=FALSE, predictor="temperature")
{
  # Description: Plot biomass in function of volume and fit linear model based on input predictor
  
  # Check number of level of input predictor
  predictor_level = levels(dataframe[,predictor])
  number_of_level = length(predictor_level)
  # Plot biomass in function of volume
  if(img_name != ""){jpeg(img_name, height=plot_height, width=plot_width)}
  color_palette = colorRampPalette(c("blue", "yellow", "red"))
  color_predictor = color_palette(number_of_level)
  plot(biomass~Volume,
       data=dataframe,
       col=color_predictor[dataframe[,predictor]],
       xlab="Volume (m3)", ylab="Biomass(g)",
       main=main,
       cex.lab=1.2,
       pch=16,
       log="xy")
  
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
    pred_biomass$biomass <- exp(predict(lm_biomass, newdata = pred_biomass))
    lines(pred_biomass$Volume, pred_biomass$biomass, col=color_predictor[i], lty=2)
  }
  # Add legend and save image
  legend(x="bottomright", legend=levels(dataframe[,predictor]), col=color_predictor, cex=0.8, pch=16, title=predictor)
  if(img_name != ""){dev.off()}
}

plot_predict_multivar <- function(dataframe, display_var, linear_model,
                                  static_vars= c("Volume", "Dim_Z", "Top_Area"),
                                  predictor="Species", plot_legend=FALSE,
                                  legend_x=0.2, legend_y=20)
{
  # Description: plot biomass in function of display var and fit linear model (keeping constant non input variable)
  
  # Color point based on species
  color_palette <- colorRampPalette(c("blue", "yellow", "red"))
  predictor_level = levels(dataframe[,predictor])
  number_of_level = length(predictor_level)
  color_predictor = color_palette(number_of_level)
  
  # Plot base points
  plot(biomass~ dataframe[,display_var],
       data=dataframe, log="xy", pch=16,
       col=species_color, xlab="", ylab="")
  
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
    pred_biomass$biomass = exp(predict(linear_model, newdata = pred_biomass))
    lines(pred_biomass[,display_var], pred_biomass$biomass, col=color_predictor[i], lty=2)
  }
  if(plot_legend)
  {
    legend(legend_x, legend_y, legend=species_level, col=species_color, pch=16, bty="n", horiz=TRUE, xjust=0.5)
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

# Add species as variable
species_label <- str_extract(plant_correlation$label, "^[A-Z]{2}")
species_corresp <- c(BR= "Brassica rapa",
                     HS="Hordeum spontaneum",
                     HV="Hordeum vulgare",
                     NB="Nicotiana benthamiana",
                     SD="Solanum dulcamara",
                     SL="Solanum lycopersicum",
                     TA="Thlaspi arvense")
#plant_correlation$Species <- species_corresp[species_label]
plant_correlation$Species <- species_label
plant_correlation$Species <- as.factor(plant_correlation$Species)

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
# Initial data
print(paste0("Initial prepared model  = ", nrow(plant_volume)))
print(paste0("No error prepared model = ", nrow(plant_correlation)))
print(paste0("No bad prepared model   = ", nrow(plant_correlation_no_bad)))
print(paste0("No issue prepared model = ", nrow(plant_correlation_no_issue)))
# Check final prepared data
summary(plant_correlation)
str(plant_correlation)

# ---- Data analysis ----
## ==== Biomass ~ Volume + Temperature ====
img_name = ifelse(save_plot, paste0("Biomass_Volume_Temperature_", output_label, "_AllPoints.jpg"), "")
plot_predict(plant_correlation, img_name=img_name, print_model=TRUE,
             main="Temperature treatment - All datapoints")
img_name = ifelse(save_plot, paste0("Biomass_Volume_Temperature_", output_label, "_NoBadIssue.jpg"), "")
plot_predict(plant_correlation_no_bad, img_name=img_name, print_model=TRUE,
             main="Temperature treatment - No bad datapoints")
img_name = ifelse(save_plot, paste0("Biomass_Volume_Temperature_", output_label, "_NoBadIssue.jpg"), "")
plot_predict(plant_correlation_no_issue, img_name=img_name, print_model=TRUE,
             main="Temperature treatment - No issue datapoint")

## ==== Biomass ~ Volume + Species ====
img_name = ifelse(save_plot, paste0("Biomass_Volume_Species_", output_label, "_AllPoints.jpg"), "")
plot_predict(plant_correlation, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - All datapoints")
img_name = ifelse(save_plot, paste0("Biomass_Volume_Species_", output_label, "_NoBadIssue.jpg"), "")
plot_predict(plant_correlation_no_bad, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - No bad datapoints")
img_name = ifelse(save_plot, paste0("Biomass_Volume_Species_", output_label, "_NoModelIssue.jpg"), "")
plot_predict(plant_correlation_no_issue, predictor="Species", img_name=img_name, print_model=TRUE,
             main="Species treatment - No issue datapoint")

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
if(save_plot){jpeg(paste0("Correlation_matrix_", output_label, ".jpg"), height=plot_height, width=plot_width)}
corrplot(correlation_mat$r, type = "upper", order = "hclust", 
         p.mat = correlation_mat$P, sig.level = 0.01, insig = "blank",
         tl.col = "Black")
if(save_plot){dev.off()}

## ==== Biomass ~ Volume + Dim_Z + Top_Area + Species ====
lm_biomass_volume_top <- lm(log(biomass) ~ log(Volume) + log(Dim_Z) + log(Top_Area) + Species, data=plant_correlation)
simulateResiduals(lm_biomass_volume_top, plot=TRUE)
summary(lm_biomass_volume_top)
drop1(lm_biomass_volume_top, test="F")

# Set Species color
color_palette <- colorRampPalette(c("blue", "yellow", "red"))
species_level <- levels(plant_correlation$Species)
species_color <- color_palette(length(species_level))

# Plot Volume, Dim_Z and Top_Area in function of each other
if(save_plot){jpeg(paste0("DependantVar_Plot_", output_label, ".jpg"), height=plot_height, width=plot_width)}
par(mfrow=c(2, 2), mar=c(1.5, 1.5, 1, 1), oma=c(2.5, 2.5, 0, 0))
plot(Volume~Dim_Z, data=plant_correlation, log="xy", pch=16, col=species_color)
plot(Volume~Top_Area, data=plant_correlation, log="xy", pch=16, col=species_color)
plot.new()
legend("bottom", legend=species_level, col=species_color, pch=16, bty="n", horiz=TRUE)
plot(Dim_Z~Top_Area, data=plant_correlation, log="xy", pch=16, col=species_color)
mtext("Volume (m3)", side=2, adj=0.8, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Plant height (m)", side=2, adj=0.2, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1.2, col="black", outer=TRUE)
if(save_plot){dev.off()}

## ==== Biomass ~ Dim_Z + Top_Area + Species ====
lm_biomass_top <- lm(log(biomass) ~ log(Dim_Z) + log(Top_Area) + Species, data=plant_correlation)
simulateResiduals(lm_biomass_top, plot=TRUE)
summary(lm_biomass_top)
drop1(lm_biomass_top, test="F")

# Plot Biomass in function of Volume, Dim_Z and Top_Area
if(save_plot){jpeg(paste0("Biomass_Volume_DimZ_TopArea_", output_label, ".jpg"), height=plot_height, width=plot_width)}
par(mfrow=c(1, 3), mar=c(1.5, 1.5, 2.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA)
plot_predict_multivar(plant_correlation, "Volume", lm_biomass_volume_top)
plot_predict_multivar(plant_correlation, "Dim_Z", lm_biomass_volume_top, plot_legend=TRUE, legend_y=25)
plot_predict_multivar(plant_correlation, "Top_Area", lm_biomass_volume_top)
mtext("Biomas (g)", side=2, agj=0.5, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Volume (m3)", side=1, adj=0.1, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.5, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.9, line=1, cex=1.2, col="black", outer=TRUE)
if(save_plot){dev.off()}

# Plot Biomass in function of Dim_Z and Top_Area only
if(save_plot){jpeg(paste0("Biomass_DimZ_TopArea_", output_label, ".jpg"), height=plot_height, width=plot_width)}
par(mfrow=c(1, 2), mar=c(1.5, 1.5, 2.5, 1), oma=c(2.5, 2.5, 0, 0), xpd=NA)
plot_predict_multivar(plant_correlation, "Dim_Z", lm_biomass_top, plot_legend=TRUE, legend_x = 1, legend_y=28)
plot_predict_multivar(plant_correlation, "Top_Area", lm_biomass_top)
mtext("Biomas (g)", side=2, agj=0.5, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Plant height (m)", side=1, adj=0.2, line=1, cex=1.2, col="black", outer=TRUE)
mtext("Top area (m2)", side=1, adj=0.8, line=1, cex=1.2, col="black", outer=TRUE)
if(save_plot){dev.off()}

# ---- Reset ----
# Reset graphic display
par(par_init)
