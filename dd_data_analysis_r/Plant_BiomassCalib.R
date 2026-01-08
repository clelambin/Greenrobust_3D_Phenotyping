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
output_folder <- "output//Pilot01"
save_plot     <- FALSE
# Set default arguments for plot display (par) and image export (jpeg)
jpeg_args     <- list(height=4, width=6, units="in", res=300)
# mgp=c(title, labels, line) set the distance for the axis (default is (3, 2, 0))
par_args      <- list(cex = 1, cex.axis=0.7, mgp=c(2, 0.8, 0))

# Libraries
library(stringr)      # Regular expression
library(DHARMa)       # Check model assumption
library(Hmisc)        # Cross-correlation with significance test
library(corrplot)     # Correlation plot
library(stats)        # For AIC

# ---- User function ----
source("Plant_BiomassCalib_Functions.R")

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
#plant_correlation$species <- species_corresp[species_label]
plant_volume$species <- species_label
plant_volume$species <- as.factor(plant_volume$species)

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
corr_species <- plant_correlation$species %in% species_11cm
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
plant_correlation$species <- droplevels(plant_correlation$species)

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
table_init <- table(plant_volume$species)
table_corr <- table(plant_correlation$species)
table_nobad   <- table(plant_correlation_no_bad$species)
table_noissue <- table(plant_correlation_no_issue$species)

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
## ==== Data overwrite ====
# (Optional, used to compare model with filtered data)
# plant_correlation <- plant_correlation_no_issue

## ==== Cross-correlation ====
# From: https://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# Note: Pearson vs Spearman (https://www.reddit.com/r/statistics/comments/76iw0w/how_do_you_decide_which_correlation_coefficient/)
#       Pearson measures the linearity of data
#       Spearman is computed on ranks and measures the monotonicity of data (doesn't necessarily have to be linear)
correlation_data <- as.matrix(subset(plant_correlation, select=-c(Path, Plant.name, Code_Error, Model.issue, Cause, label, species, temperature)))
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

## ==== Graph settings ====
# Set species color
color_palette <- colorRampPalette(c("blue", "yellow", "red"))
species_level <- levels(plant_correlation$species)
species_color <- color_palette(length(species_level))

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

## ==== Plant trait distribution ====
# Plot distribution of the different plant traits per species
img_name <- paste0(output_folder, "//Biomass_distribution_per_species_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(biomass) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Biomass (g))")
if(save_plot){dev.off()}

# Plot distribution of the different plant traits per treatment
img_name <- paste0(output_folder, "//Biomass_distribution_per_treatment_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(biomass) ~ temperature, data=plant_correlation, col=c("darkblue", "lightblue", "yellow", "orangered", "darkred"),
        xlab="", ylab="log(Biomass (g))", names=c("22°C", "26°C", "30°C", "34°C", "38°C"))
if(save_plot){dev.off()}

# Biomass stats in function of treatment and species
by(plant_correlation, plant_correlation$temperature, function(subframe)
  {
    stat<-c(mean(subframe$biomass), var(subframe$biomass), sd(subframe$biomass))
    names(stat)<-c("mean", "var", "std dev")
    stat
  }, simplify = TRUE)
by(plant_correlation, plant_correlation$species, function(subframe)
  {
    stat<-c(mean(subframe$biomass), var(subframe$biomass), sd(subframe$biomass))
    names(stat)<-c("mean", "var", "std dev")
    stat
  }, simplify = TRUE)

## ==== Volume, Dim_Z, Top_Area ====
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

## ==== (log) Biomass ~ species ====
# Fit linear model only dependent on species to use as baseline for comparison
lm_biomass_species_log <- lm(log(biomass) ~ species, data=plant_correlation)
model_stats(lm_biomass_species_log)
img_name = paste0(output_folder, "//Biomass_Species_LogAxis_", output_label, "_AllPoints.jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass), mdl=lm_biomass_species_log,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ Volume + species ====
# Create prediction for the different filter level to see impact of outliers
lm_biomass_volume_log <- lm(log(biomass) ~ log(Volume) + species, data=plant_correlation)
model_stats(lm_biomass_volume_log)
img_name = paste0(output_folder, "//Biomass_Volume_Species_LogAxis_", output_label, "_AllPoints.jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass), mdl=lm_biomass_volume_log,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ Dim_Z + Top_Area + species ====
lm_biomass_top <- lm(log(biomass) ~ log(Dim_Z) + log(Top_Area) + species, data=plant_correlation)
model_stats(lm_biomass_top)
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_LinearModel_", output_label, ".jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass), mdl=lm_biomass_top,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log with interact) Biomass ~ Dim_Z + Top_Area + species ====
# Test with interaction
lm_biomass_top_with_interact <- lm(log(biomass) ~ (log(Dim_Z) + log(Top_Area)) * species, data=plant_correlation)
model_stats(lm_biomass_top_with_interact)
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_Interaction_LinearModel_", output_label, ".jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass), mdl=lm_biomass_top_with_interact,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
# To check: visreg
# Partial residual plot

## ==== (sqrt) Biomass ~ Dim_Z + Top_Area + species ====
lm_biomass_top_sqrt <- lm(sqrt(biomass) ~ sqrt(Dim_Z) + sqrt(Top_Area) + species, data=plant_correlation)
model_stats(lm_biomass_top_sqrt)
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_LinearModel_Sqrt_", output_label, ".jpg")
plot_obs_vs_pred(obs=sqrt(plant_correlation$biomass), mdl=lm_biomass_top_sqrt,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (sqrt) Biomass ~ Dim_Z + Top_Area ====
lm_bmass_z_top_sqrt_nospecies <- lm(sqrt(biomass) ~ sqrt(Dim_Z) + sqrt(Top_Area), data=plant_correlation)
model_stats(lm_bmass_z_top_sqrt_nospecies)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_DimZ_TopArea_NoSpecies_LinearModel_Sqrt_", output_label, ".jpg")
plot_obs_vs_pred(obs=sqrt(plant_correlation$biomass), mdl=lm_bmass_z_top_sqrt_nospecies,
                 img_name=img_name, xlab="Predicted sqrt(biomass(g))", ylab="Measured sqrt(biomass(g))")

## ==== Model comparison ====
model_to_compare <- list(lm_biomass_volume_log, lm_biomass_top, lm_biomass_top_sqrt, lm_biomass_top_with_interact)
model_names <- c("Biomass~volume (log)", "Biomass~height+area (log)", "Biomass~height+area (sqrt)", "Biomass~height+area (with interact)")
model_AIC <- sapply(model_to_compare, FUN=AIC)
names(model_AIC) <- model_names
model_AIC

# ---- Reset ----
# Reset graphic display
par(par_init)
