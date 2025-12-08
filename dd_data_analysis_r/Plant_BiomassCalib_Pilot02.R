# ---- Header ----
# Name: plant_harvestCalib_Pilot02
# Description: Calibrate plant biomass based on data extracted from 3D model from 2nd pilot

# ---- Preparation ----

#setting working directory & clean up environment
rm(list=ls())
# Reset display and save initial par
if(length(dev.list()!=0)){dev.off()}
par_init <- par(no.readonly = TRUE)

# Script variables
input_from3d  <- "data//Plant_data_20251208_pilot02.csv"
input_harvest <- "data//20251208_ppa_drought.csv"
input_organes <- "data//Plant_organes_count_20251202.csv"
output_label  <- "20251208"
output_folder <- "output//pilot02"
save_plot     <- FALSE
# Set default arguments for plot display (par) and image export (jpeg)
jpeg_args     <- list(height=4, width=6, units="in", res=300)
# mgp=c(title, labels, line) set the distance for the axis (default is (3, 2, 0))
par_args      <- list(cex = 1, cex.axis=0.7, mgp=c(2, 0.8, 0))
# Regular expression pattern
indiv_label   <- "[A-Z]{2}[0-9]{2}I[CD]"
# Relevant columns to use
filter_collect <- c("pot_id", "species", "treatment", "aboveground_dry_biomass_g", "Height..cm.", "Nb.leaves", "Nb.Flowers", "Leaf.area..cm2.")
filter_from3d  <- c("Volume", "Dim_X", "Dim_Y", "Dim_Z", "Height", "Cumul_Area", "Font_Area_Ratio", "Left_Area_Ratio", "Top_Area_Ratio")


# Libraries
library(stringr)      # Regular expression
library(DHARMa)       # Check model assumption
library(Hmisc)        # Cross-correlation with significance test
library(corrplot)     # Correlation plot
library(stats)        # For AIC

# ---- User function ----
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

# ---- Data preparation ----
## ==== Data extraction ====
# Load data
plant_harvest <- read.csv(input_harvest)
plant_organes <- read.csv(input_organes)
plant_from3d  <- read.csv(input_from3d)

# Extract label from model name
plant_from3d$label <- str_extract(plant_from3d$Plant.name, indiv_label)
# Exclude population pot from the analysis (script not handling them properly yet)
plant_from3d <- plant_from3d[!is.na(plant_from3d$label),]
plant_harvest <- plant_harvest[plant_harvest$orga_level!="Population",]

## ==== Handle duplicate ====
# If a plant has been scanned multiple time, use the last time only
# - plant_from3d
# (first extract datecode and reorder by descending date, because duplicated takes the first element)
plant_from3d$datecode <- str_extract(plant_from3d$Plant.name, "2025[0-9]{4}")
plant_from3d <- plant_from3d[order(plant_from3d$datecode, decreasing = TRUE),]
# (then use duplicated to remove all but the 1st entry of each duplicated species)
plant_from3d <- plant_from3d[!duplicated(plant_from3d$label),]
# - plant_organes
harvest_date_list <- strsplit(plant_organes$Harvest.date, "/")
plant_organes$datecode <- unlist(lapply(harvest_date_list, function(l) paste(l[c(3,2,1)], collapse="")))
plant_organes <- plant_organes[order(plant_organes$datecode, decreasing = TRUE),]
plant_organes <- plant_organes[!duplicated(plant_organes$Plant),]

## ==== Merge dataframes ====
# Merge all 3 dataframe
plant_datacollection <- merge(plant_harvest, plant_organes, by.x="pot_id", by.y="Plant", all=FALSE)
# plant_datacollection <- subset(plant_datacollection, select = filter_collect)
plant_correlation <- merge(plant_datacollection, plant_from3d, by.x="pot_id", by.y="label", all=FALSE)
plant_correlation <- subset(plant_correlation, select = c(filter_collect, filter_from3d))

# Rename variables
names(plant_correlation)[names(plant_correlation)=="aboveground_dry_biomass_g"] <- "biomass_g"

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

# Save a backup before data cleanup to compare
plant_init        <- plant_correlation
# Remove outliers
# (To check why outliers occurs in script)
plant_correlation <- plant_correlation[plant_correlation$Dim_X > 0 & plant_correlation$Dim_Y > 0 & plant_correlation$Dim_Z > 0,]
# Remove missing data
plant_correlation <- plant_correlation[!is.na(plant_correlation$biomass_g),]

## ==== Convert to factors ====
plant_harvest$species <- as.factor(plant_harvest$species)
# plant_organes$species <- as.factor(plant_organes$Plant)
# plant_from3d$label <- as.factor(plant_from3d$label)
plant_correlation$species <- as.factor(plant_correlation$species)
plant_correlation$treatment <- as.factor(plant_correlation$treatment)

## ==== Overview ====
# Data summary
summary(plant_correlation)

# Processed plant overview (per species)
all_species   <- levels(plant_harvest$species)
table_harvest <- table(factor(plant_harvest$species, levels=all_species))
# table_organes <- table(factor(plant_organes$Plant, levels=all_species))
# table_from3d  <- table(factor(plant_from3d$label, levels=all_species))
table_init    <- table(factor(plant_init$species, levels=all_species))
table_corr    <- table(factor(plant_correlation$species, levels=all_species))
table_combine <- rbind(c(table_harvest, sum(table_harvest)),
                       # c(table_organes, sum(table_organes)),
                       # c(table_from3d, sum(table_from3d)),
                       c(table_init, sum(table_init)),
                       c(table_corr, sum(table_corr)))
colnames(table_combine)[ncol(table_combine)] <- "sum"
# rownames(table_combine) <- c("Havest", "Organes counting", "3D data", "Merged data", "Merged, no outliers")
rownames(table_combine) <- c("Havest", "Merged data", "Merged, no missing")
table_combine

# ---- Data analysis global ----
## ==== Cross-correlation ====
# From: https://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# Note: Pearson vs Spearman (https://www.reddit.com/r/statistics/comments/76iw0w/how_do_you_decide_which_correlation_coefficient/)
#       Pearson measures the linearity of data
#       Spearman is computed on ranks and measures the monotonicity of data (doesn't necessarily have to be linear)
correlation_data <- as.matrix(subset(plant_correlation, select=-c(pot_id, species, treatment, Leaf.area..cm2.,
                                                                  Font_Area_Ratio, Top_Area_Ratio, Left_Area_Ratio,
                                                                  Front_BBox, Top_BBox, Left_BBox)))
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
# Set Species color
color_palette <- colorRampPalette(c("blue", "yellow", "red"))
species_level <- levels(plant_correlation$species)
species_color <- color_palette(length(species_level))

## ==== Height ~ Height ====
lm_height <- lm(Height..cm. ~ Height, data=plant_correlation)
model_stats(lm_height)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Height_Predict_", output_label, ".jpg")
plot_obs_vs_pred(obs=plant_correlation$Height..cm., mdl=lm_height,
                 img_name=img_name, xlab="Predicted height(cm)", ylab="Measured height(cm)")

## ==== Leaf_Area ~ Cumul_Area + Species ====
plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.),])
# plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.) & plant_correlation$species=="Nicotiana benthamiana",])
# plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.) & plant_correlation$species=="Hordeum vulgare",])
lm_area_log <- lm(Leaf.area..cm2. ~ Cumul_Area + species, data=plant_lf_area)
# lm_area_log <- lm(Leaf.area..cm2. ~ Cumul_Area, data=plant_lf_area)
model_stats(lm_area_log)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Leaf_Area_Species_Predict_", output_label, ".jpg")
plot_obs_vs_pred(obs=plant_lf_area$Height..cm., mdl=lm_area_log, data=plant_lf_area,
                 img_name=img_name, xlab="Predicted leaf area(cm2)", ylab="Measured leaf area(cm2)")

# ---- Data analysis biomass ----
## ==== (log) Biomass ~ Volume + Species ====
lm_bmass_vol_log <- lm(log(biomass_g) ~ log(Volume) + species, data=plant_correlation)
model_stats(lm_bmass_vol_log)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Volume_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass_g), mdl=lm_bmass_vol_log,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ Dim_Z + Top_Area + Species ====
lm_bmass_z_top_log <- lm(log(biomass_g) ~ log(Dim_Z) + log(Top_Area) + species, data=plant_correlation)
model_stats(lm_bmass_z_top_log)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass_g), mdl=lm_bmass_z_top_log,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (sqrt) Biomass ~ Dim_Z + Top_Area + Species ====
lm_bmass_z_top_sqrt <- lm(sqrt(biomass_g) ~ sqrt(Dim_Z) + sqrt(Top_Area) + species, data=plant_correlation)
model_stats(lm_bmass_z_top_sqrt)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Species_Predict_Sqrt_", output_label, ".jpg")
plot_obs_vs_pred(obs=sqrt(plant_correlation$biomass_g), mdl=lm_bmass_z_top_sqrt,
                 img_name=img_name, xlab="Predicted sqrt(biomass(g))", ylab="Measured sqrt(biomass(g))")

## ==== (log) Biomass ~ Dim_Z + Top_Area ====
lm_bmass_z_top_log_nospecies <- lm(log(biomass_g) ~ log(Dim_Z) + log(Top_Area), data=plant_correlation)
model_stats(lm_bmass_z_top_log_nospecies)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(obs=log(plant_correlation$biomass_g), mdl=lm_bmass_z_top_log_nospecies,
                 img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (sqrt) Biomass ~ Dim_Z + Top_Area ====
lm_bmass_z_top_sqrt_nospecies <- lm(sqrt(biomass_g) ~ sqrt(Dim_Z) + sqrt(Top_Area), data=plant_correlation)
model_stats(lm_bmass_z_top_sqrt_nospecies)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Predict_Sqrt_", output_label, ".jpg")
plot_obs_vs_pred(obs=sqrt(plant_correlation$biomass_g), mdl=lm_bmass_z_top_sqrt_nospecies,
                 img_name=img_name, xlab="Predicted sqrt(biomass(g))", ylab="Measured sqrt(biomass(g))")

## ==== (sqrt) Biomass ~ Dim_Z * Top_Area ====
lm_bmass_z_top_sqrt_nospecies_interact <- lm(sqrt(biomass_g) ~ sqrt(Dim_Z) * sqrt(Top_Area), data=plant_correlation)
model_stats(lm_bmass_z_top_sqrt_nospecies_interact)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_inter_TopArea_Predict_Sqrt_", output_label, ".jpg")
plot_obs_vs_pred(obs=sqrt(plant_correlation$biomass_g), mdl=lm_bmass_z_top_sqrt_nospecies_interact,
                 img_name=img_name, xlab="Predicted sqrt(biomass(g))", ylab="Measured sqrt(biomass(g))")

# Analysis per species (skipped because bad correlation)
if(FALSE)
{
  for(species in species_level)
  {
    # ---- Data analysis per species ----
    # Create filtered dataframe
    plant_perspecies <- droplevels(plant_correlation[plant_correlation$species==species,])
    # Replace space by underscore for file output
    species_string <- replace(species, " ", "_")
    
    ## ==== (log) Biomass ~ Dim_Z + Top_Area ====
    # Model
    lm_bmass_z_top_log_perspecies <- lm(log(biomass_g) ~ log(Dim_Z) + log(Top_Area), data=plant_perspecies)
    model_stats(lm_bmass_z_top_log_perspecies, plot_res=FALSE)
    # Observed vs predicted plot
    img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Predict_Sqrt_", output_label, ".jpg")
    plot_obs_vs_pred(obs=log(plant_perspecies$biomass_g), mdl=lm_bmass_z_top_log_perspecies, data=plant_perspecies,
                     img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
    
    ## ==== (sqrt) Biomass ~ Dim_Z + Top_Area ====
    # Model
    lm_bmass_z_top_sqrt_perspecies <- lm(sqrt(biomass_g) ~ sqrt(Dim_Z) + sqrt(Top_Area), data=plant_perspecies)
    model_stats(lm_bmass_z_top_sqrt_perspecies, plot_res=FALSE)
    # Observed vs predicted plot
    img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Predict_Sqrt_", output_label, ".jpg")
    plot_obs_vs_pred(obs=sqrt(plant_perspecies$biomass_g), mdl=lm_bmass_z_top_sqrt_perspecies, data=plant_perspecies,
                     img_name=img_name, xlab="Predicted sqrt(biomass(g))", ylab="Measured sqrt(biomass(g))")
    
  }
}
