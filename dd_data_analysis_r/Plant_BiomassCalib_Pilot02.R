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
input_from3d  <- "data//Plant_data_20251210_pilot02_stickrm.csv"
input_harvest <- "data//20251208_ppa_drought.csv"
input_organes <- "data//Plant_organes_count_20251202.csv"
output_label  <- "20251210"
output_folder <- "output//Pilot02"
save_plot     <- FALSE
plot_dharma   <- FALSE
distrib_plot  <- TRUE
species_plot  <- FALSE
# Set default arguments for plot display (par) and image export (jpeg)
jpeg_args     <- list(height=4, width=6, units="in", res=300)
# mgp=c(title, labels, line) set the distance for the axis (default is (3, 2, 0))
par_args      <- list(cex = 1, cex.axis=0.7, mgp=c(2, 0.8, 0))
# Regular expression pattern
indiv_label   <- "[A-Z]{2}[0-9]{2}I[CD]"
# Relevant columns to use
filter_collect <- c("pot_id", "species", "treatment", "aboveground_dry_biomass_g", "Height..cm.", "Nb.leaves", "Nb.Flowers", "Leaf.area..cm2.")
filter_from3d  <- c("Volume", "Dim_X", "Dim_Y", "Dim_Z", "Height", "Cumul_Area", "Font_Area_Ratio", "Left_Area_Ratio", "Top_Area_Ratio", "Convex_hull", "Convex_hull_40", "Convex_hull_60")
# Apply filter based on height prediction deviation (in %)
height_filter_thresh <- 0.2

# Script functions
label_var  <- c(
  "Cumul_Area"="Cumulative area (m2)",
  "Volume"="Volume (m3)",
  "Convex_hull"="Volume of convex hull (m3)",
  "Convex_hull_40"="Volume of 40% convex hull (m3)",
  "Convex_hull_60"="Volume of 60% convex hull (m3)",
  "Top_Area"="Projected top area (m2)",
  "Height"="Plant height (m)"
)
fct_label <- c(ident="", log="log of ", sqrt="sqrt of ")


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
plant_correlation <- plant_correlation[plant_correlation$Dim_X > 0 & plant_correlation$Dim_Y > 0 & plant_correlation$Dim_Z > 0 & plant_correlation$Height > 0,]
# Remove missing data
plant_correlation <- plant_correlation[!is.na(plant_correlation$biomass_g),]

## ==== Convert to factors ====
species_corresp <- c("Arabidopsis thaliana"="AT",
                     "Brachypodium distachyon"="BD",
                     "Brassica rapa"="BR",
                     "Hordeum spontaneum"="HS",
                     "Hordeum vulgare"="HV",
                     "Nicotiana benthamiana"="NB",
                     "Solanum dulcamara"="SD",
                     "Solanum lycopersicum"="SL",
                     "Thlaspi arvense"="TA")
treatment_corresp <- c("drought"="D",
                       "control"="C")
plant_harvest$species <- as.factor(species_corresp[plant_harvest$species])
plant_init$species    <- as.factor(species_corresp[plant_init$species])
plant_correlation$species <- as.factor(species_corresp[plant_correlation$species])
plant_correlation$treatment <- as.factor(treatment_corresp[plant_correlation$treatment])
# plant_harvest$species <- as.factor(plant_harvest$species)
# plant_correlation$species <- as.factor(plant_correlation$species)
# plant_correlation$treatment <- as.factor(plant_correlation$treatment)

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

## ==== Height ~ Species ====
lm_height_species_only <- lm(Height..cm. ~ species , data=plant_correlation)
model_stats(lm_height_species_only, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Height_Species_only_Predict_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_height_species_only, img_name=img_name, xlab="Predicted height(cm)", ylab="Measured height(cm)")

## ==== Height ~ Height ====
lm_height <- lm(Height..cm. ~ Height, data=plant_correlation)
model_stats(lm_height, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Height_Predict_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_height, img_name=img_name, xlab="Predicted height(cm)", ylab="Measured height(cm)")

## ==== Height ~ Height + species ====
lm_height_species <- lm(Height..cm. ~ Height+species, data=plant_correlation)
model_stats(lm_height_species, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Height_species_Predict_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_height_species, img_name=img_name, xlab="Predicted height(cm)", ylab="Measured height(cm)")

## ==== Height ~ Height * species ====
lm_height_species_inter <- lm(Height..cm. ~ Height*species, data=plant_correlation)
model_stats(lm_height_species_inter, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Height_inter_species_Predict_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_height_species_inter, img_name=img_name, xlab="Predicted height(cm)", ylab="Measured height(cm)")

## ==== Model comparison ====
mdl_compare(list(lm_height, lm_height_species_only, lm_height_species, lm_height_species_inter))

## ==== Plant trait distribution ====
# Plot distribution of the different plant traits per species
img_name <- paste0(output_folder, "//Trait_distribution_per_species_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(2, 2), mar=c(2, 3, 1, 1)), par_args))
boxplot(Height..cm. ~ species, data=plant_correlation, col=species_color, xlab="", ylab="Height (cm)")
boxplot(log(biomass_g) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Biomass (g))")
boxplot(log(Nb.leaves) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Leaf count)")
boxplot(log(Nb.Flowers) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Flower count)")
if(save_plot){dev.off()}

img_name <- paste0(output_folder, "//Trait_distribution_per_treatment_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(2, 2), mar=c(2, 3, 1, 1)), par_args))
boxplot(Height..cm. ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="Height (cm)", names=c("Control", "Drought"))
boxplot(log(biomass_g) ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="log(Biomass (g))", names=c("Control", "Drought"))
boxplot(log(Nb.leaves) ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="log(Leaf count)", names=c("Control", "Drought"))
boxplot(log(Nb.Flowers) ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="log(Flower count)", names=c("Control", "Drought"))
if(save_plot){dev.off()}

# Individual plots
# Plot distribution of the different plant traits per species
img_name <- paste0(output_folder, "//Height_distribution_per_species_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(Height..cm. ~ species, data=plant_correlation, col=species_color, xlab="", ylab="Height (cm)")
if(save_plot){dev.off()}
img_name <- paste0(output_folder, "//Biomass_distribution_per_species_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(biomass_g) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Biomass (g))")
if(save_plot){dev.off()}
img_name <- paste0(output_folder, "//NbLeaves_distribution_per_species_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(Nb.leaves) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Leaf count)")
if(save_plot){dev.off()}
img_name <- paste0(output_folder, "//NbFlowers_distribution_per_species_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(Nb.Flowers) ~ species, data=plant_correlation, col=species_color, xlab="", ylab="log(Flower count)")
if(save_plot){dev.off()}

# Plot distribution of the different plant traits per treatment
img_name <- paste0(output_folder, "//Height_distribution_per_treatment_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(Height..cm. ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="Height (cm)", names=c("Control", "Drought"))
if(save_plot){dev.off()}
img_name <- paste0(output_folder, "//Biomass_distribution_per_treatment_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(biomass_g) ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="log(Biomass (g))", names=c("Control", "Drought"))
if(save_plot){dev.off()}
img_name <- paste0(output_folder, "//NbLeaves_distribution_per_treatment_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(Nb.leaves) ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="log(Leaf count)", names=c("Control", "Drought"))
if(save_plot){dev.off()}
img_name <- paste0(output_folder, "//NbFlowers_distribution_per_treatment_", output_label, ".jpg")
if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
par(c(list(mfrow=c(1, 1), mar=c(2, 3, 1, 1)), par_args))
boxplot(log(Nb.Flowers) ~ treatment, data=plant_correlation, col=c("lightblue", "orangered"),
        xlab="", ylab="log(Flower count)", names=c("Control", "Drought"))
if(save_plot){dev.off()}

# Biomass stats in function of treatment and species
by(plant_correlation, plant_correlation$treatment, function(subframe)
  {
    stat<-c(mean(subframe$biomass_g), var(subframe$biomass_g), sd(subframe$biomass_g))
    names(stat)<-c("mean", "var", "std dev")
    stat
  }, simplify = TRUE)
by(plant_correlation, plant_correlation$species, function(subframe)
  {
    stat<-c(mean(subframe$biomass_g), var(subframe$biomass_g), sd(subframe$biomass_g))
    names(stat)<-c("mean", "var", "std dev")
    stat
  }, simplify = TRUE)

# ## ==== Leaf_Area ~ Cumul_Area + Species ====
# plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.),])
# # plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.) & plant_correlation$species=="Nicotiana benthamiana",])
# # plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.) & plant_correlation$species=="Hordeum vulgare",])
# lm_area_log <- lm(Leaf.area..cm2. ~ Cumul_Area + species, data=plant_lf_area)
# # lm_area_log <- lm(Leaf.area..cm2. ~ Cumul_Area, data=plant_lf_area)
# model_stats(lm_area_log, plot_res=plot_dharma)
# 
# # Observed vs predicted plot
# img_name <- paste0(output_folder, "//Leaf_Area_Species_Predict_", output_label, ".jpg")
# plot_obs_vs_pred(obs=plant_lf_area$Leaf.area..cm2., mdl=lm_area_log, data=plant_lf_area,
#                  img_name=img_name, xlab="Predicted leaf area(cm2)", ylab="Measured leaf area(cm2)")
# 
# 
# ## ==== Leaf_Area ~ Cumul_Area + Species ====
# plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.),])
# # plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.) & plant_correlation$species=="Nicotiana benthamiana",])
# # plant_lf_area <- droplevels(plant_correlation[!is.na(plant_correlation$Leaf.area..cm2.) & plant_correlation$species=="Hordeum vulgare",])
# lm_area_log <- lm(Leaf.area..cm2. ~ Cumul_Area * species, data=plant_lf_area)
# # lm_area_log <- lm(Leaf.area..cm2. ~ Cumul_Area, data=plant_lf_area)
# model_stats(lm_area_log, plot_res=plot_dharma)
# 
# # Observed vs predicted plot
# img_name <- paste0(output_folder, "//Leaf_Area_Species_Predict_", output_label, ".jpg")
# plot_obs_vs_pred(obs=plant_lf_area$Leaf.area..cm2., mdl=lm_area_log, data=plant_lf_area,
#                  img_name=img_name, xlab="Predicted leaf area(cm2)", ylab="Measured leaf area(cm2)")

# ---- Data analysis biomass ----
## ==== (log) Biomass ~ Species ====
lm_bmass_spec_log <- lm(log(biomass_g) ~ species, data=plant_correlation)
model_stats(lm_bmass_spec_log, plot_res=plot_dharma)
img_name <- paste0(output_folder, "//Biomass_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_spec_log, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ Volume + Species ====
lm_bmass_vol_log <- lm(log(biomass_g) ~ log(Volume) + species, data=plant_correlation)
model_stats(lm_bmass_vol_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Volume_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_vol_log, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ Height + Top_Area + Species ====
lm_bmass_z_top_log <- lm(log(biomass_g) ~ log(Height) + log(Top_Area) + species, data=plant_correlation)
model_stats(lm_bmass_z_top_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_z_top_log, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ (Height + Top_Area) * Species ====
lm_bmass_z_top_inter_log <- lm(log(biomass_g) ~ (log(Height) + log(Top_Area)) * species, data=plant_correlation)
model_stats(lm_bmass_z_top_inter_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_inter_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_z_top_inter_log, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ Height + Top_Area ====
lm_bmass_z_top_log_nospecies <- lm(log(biomass_g) ~ log(Height) + log(Top_Area), data=plant_correlation)
model_stats(lm_bmass_z_top_log_nospecies, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_z_top_log_nospecies, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== (log) Biomass ~ 60% Convex Hull * Species ====
lm_bmass_chull_inter_log <- lm(log(biomass_g) ~ log(Convex_hull_60) * species, data=plant_correlation)
model_stats(lm_bmass_chull_inter_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Chull60_inter_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_chull_inter_log, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

## ==== Model comparison ====
mdl_compare(list(lm_bmass_vol_log, lm_bmass_z_top_log, lm_bmass_z_top_inter_log, lm_bmass_z_top_log_nospecies, lm_bmass_chull_inter_log, lm_bmass_spec_log))

# ---- Data analysis height filter ----
# Filter data where the height prediction is too far from the measured height
height_deviation <- abs(predict(lm_height) - plant_correlation$Height..cm.) / plant_correlation$Height..cm.
paste0("Ratio of points with deviation below ", height_filter_thresh, ": ", round(mean(height_deviation<height_filter_thresh), 3))
# Create dataframe with deviation above threshold
plant_filtered <- plant_correlation
plant_filtered$deviation <- height_deviation
plant_filtered <- droplevels(plant_filtered[height_deviation<height_filter_thresh,])
# Sort dataframe by deviation
plant_filtered <- plant_filtered[sort(plant_filtered$deviation, decreasing=TRUE, index.return=TRUE)$ix,]

## ==== Height ~ Height ====
lm_height_filt <- lm(Height..cm. ~ Height, data=plant_filtered)
model_stats(lm_height_filt, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Height_Predict_Filtered_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_height_filt, img_name=img_name, xlab="Predicted height(cm)", ylab="Measured height(cm)")

## ==== (log) Biomass ~ Height + Top_Area + Species ====
lm_bmass_z_top_log_filt <- lm(log(biomass_g) ~ log(Height) + log(Top_Area) + species, data=plant_filtered)
model_stats(lm_bmass_z_top_log_filt, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Species_Predict_Filtered_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_z_top_log_filt, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")


## ==== (log) Biomass ~ (Height + Top_Area) * Species ====
lm_bmass_z_top_inter_log_filt <- lm(log(biomass_g) ~ (log(Height) + log(Top_Area)) * species, data=plant_filtered)
model_stats(lm_bmass_z_top_inter_log_filt, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_inter_Species_Predict_Filtered_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=lm_bmass_z_top_inter_log_filt, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")

# ---- Data analysis leaf count ----
## ==== Trait distribution ====
# Initialise variables used for plotting the different trait distribution
dataframe <- plant_correlation
ident <- function(x) x
transf_fct <- c(log=log,          log=log,  log=log,    log=log)
indep_var  <- c("Convex_hull_60", "Height", "Top_Area", "Cumul_Area")

# Loop through all variable and plot distribution
if (distrib_plot)
{
  for(index in 1:length(indep_var))
  {
    fct   <- transf_fct[[index]]
    var_x <- indep_var[index]
    fct_name <- names(transf_fct)[index]
    img_name <- paste0(output_folder, "//NbLeaves_", var_x, "_", fct_name, "_", output_label, ".jpg")
    if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
    par(c(list(mar=c(4, 4, 1, 1)), par_args))
    plot(fct(dataframe$Nb.leaves)~fct(dataframe[,var_x]),
         pch=16, col=species_color[dataframe$species],
         xlab=paste0(fct_label[fct_name], label_var[var_x]),
         ylab=paste0(fct_label[fct_name], "Number of leaves"))
    legend("bottomright", legend=levels(dataframe$species), col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)
    if(save_plot){dev.off()}
  }
}

# # Individual plot
# par(c(list(mar=c(4, 4, 1, 1)), par_args))
# plot(Nb.leaves~Cumul_Area, data=plant_correlation,
#      pch=16, col=species_color[plant_correlation$species],
#      xlab="Cumulative area (m2)", ylab="Number of leaves")
# legend("bottomright", legend=levels(plant_correlation$species), col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)

## ==== Nb Leaves ~ Species ====
glm_nbleaves_spec <- glm(Nb.leaves ~ species,
                         data=plant_correlation, family = poisson())
model_stats(glm_nbleaves_spec, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbLeaves_Species_Predict_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbleaves_spec, img_name=img_name, xlab="Predicted leaf count", ylab="Measured leaf count")

## ==== Nb Leaves ~ 60% Convex Hull + Species ====
lm_nbleaves_chull_60 <- lm(Nb.leaves ~ Convex_hull_60 + species, data=plant_correlation)
model_stats(lm_nbleaves_chull_60, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbLeaves_Chull60_Species_Predict_Sqrt_", output_label, ".jpg")
plot_obs_vs_pred(obs=plant_correlation$Nb.leaves, mdl=lm_nbleaves_chull_60,
                 img_name=img_name, xlab="Predicted leaf count", ylab="Measured leaf count")

## ==== (log) Nb Leaves ~ 60% Convex Hull * Species ====
glm_nbleaves_chull_60_inter_log <- glm(Nb.leaves ~ log(Convex_hull_60) * species,
                                       data=plant_correlation, family = poisson())
model_stats(glm_nbleaves_chull_60_inter_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbLeaves_Chull60_inter_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbleaves_chull_60_inter_log, img_name=img_name, xlab="Predicted leaf count", ylab="Measured leaf count")

## ==== (log) Nb Leaves ~ Height + Top_Area + Species ====
glm_nbleaves_z_top_log <- glm(Nb.leaves ~ log(Height) + log(Top_Area) + species,
                              data=plant_correlation, family = poisson())
model_stats(glm_nbleaves_z_top_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbLeaves_Height_TopArea_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbleaves_z_top_log, img_name=img_name, xlab="Predicted leaf count", ylab="Measured leaf count")

## ==== (log) Nb Leaves ~ (Height + Top_Area) * Species ====
glm_nbleaves_z_top_inter_log <- glm(Nb.leaves ~ (log(Height) + log(Top_Area)) * species,
                                    data=plant_correlation, family = poisson())
model_stats(glm_nbleaves_z_top_inter_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbLeaves_Height_TopArea_inter_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbleaves_z_top_inter_log, img_name=img_name, xlab="Predicted leaf count", ylab="Measured leaf count")

## ==== Model comparison ====
mdl_compare(list(glm_nbleaves_chull_60_inter_log, glm_nbleaves_z_top_log, glm_nbleaves_z_top_inter_log, glm_nbleaves_spec))

# ---- Data analysis flower count ----
## ==== Trait distribution ====
# Initialize variables used for plotting the different trait distribution
dataframe <- plant_correlation
ident <- function(x) x
transf_fct <- c(log=log,          log=log,  log=log,    log=log)
indep_var  <- c("Convex_hull_60", "Height", "Top_Area", "Cumul_Area")

# Loop through all variable and plot distribution
if(distrib_plot)
{
  for(index in 1:length(indep_var))
  {
    fct   <- transf_fct[[index]]
    var_x <- indep_var[index]
    fct_name <- names(transf_fct)[index]
    img_name <- paste0(output_folder, "//NbFlowers_", var_x, "_", fct_name, "_", output_label, ".jpg")
    if(save_plot){do.call(jpeg, c(filename=img_name, jpeg_args))}
    par(c(list(mar=c(4, 4, 1, 1)), par_args))
    plot(fct(dataframe$Nb.Flowers)~fct(dataframe[,var_x]),
         pch=16, col=species_color[dataframe$species],
         xlab=paste0(fct_label[fct_name], label_var[var_x]),
         ylab=paste0(fct_label[fct_name], "Number of flowers"))
    legend("bottomright", legend=levels(dataframe$species), col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)
    if(save_plot){dev.off()}
  }
}

# # Individual plot
# par(c(list(mar=c(4, 4, 1, 1)), par_args))
# plot(Nb.Flowers~Cumul_Area, data=plant_correlation,
#      pch=16, col=species_color[plant_correlation$species],
#      xlab="Cumulative area (m2)", ylab="Number of flowers")
# legend("bottomright", legend=levels(plant_correlation$species), col=species_color, pch=16, bty="y", horiz=FALSE, cex=0.8)

## ==== Nb Flower ~ Species ====
glm_nbflowers_species <- glm(Nb.Flowers ~ species, data=plant_correlation, family=poisson())
model_stats(glm_nbflowers_species, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbFlowers_Species_Predict_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbflowers_species, img_name=img_name, xlab="Predicted flower count", ylab="Measured flower count")

## ==== (log) Nb Flowers ~ 60% Convex Hull * Species ====
glm_nbflowers_chull_60_inter_log <- glm(Nb.Flowers ~ log(Convex_hull_60) * species, data=plant_correlation, family=poisson())
model_stats(glm_nbflowers_chull_60_inter_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbFlowers_Chull60_inter_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbflowers_chull_60_inter_log, img_name=img_name, xlab="Predicted flower count", ylab="Measured flower count")

## ==== (log) Nb Flowers ~ Height + Top_Area + Species ====
glm_nbflowers_z_top_log <- glm(Nb.Flowers ~ log(Height) + log(Top_Area) + species,
                               data=plant_correlation, family = poisson())
model_stats(glm_nbflowers_z_top_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbFlowers_Height_TopArea_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbflowers_z_top_log, img_name=img_name, xlab="Predicted flower count", ylab="Measured flower count")

## ==== (log) Nb Flowers ~ (Height + Top_Area) * Species ====
glm_nbflowers_z_top_inter_log <- glm(Nb.Flowers ~ (log(Height) + log(Top_Area)) * species,
                                     data=plant_correlation, family = poisson())
model_stats(glm_nbflowers_z_top_inter_log, plot_res=plot_dharma)

# Observed vs predicted plot
img_name <- paste0(output_folder, "//NbFlowers_Height_TopArea_inter_Species_Predict_Log_", output_label, ".jpg")
plot_obs_vs_pred(mdl=glm_nbflowers_z_top_inter_log, img_name=img_name, xlab="Predicted flower count", ylab="Measured flower count")

## ==== Model comparison ====
mdl_compare(list(glm_nbflowers_z_top_inter_log, glm_nbflowers_z_top_log, glm_nbflowers_species))

# = = = = = = =
# Analysis per species (skipped because bad correlation)
# dataframe = plant_correlation
dataframe = plant_filtered
if(species_plot)
{
  for(species in species_level)
  {
    # ---- Data analysis per species ----
    # Create filtered dataframe
    plant_perspecies <- droplevels(dataframe[dataframe$species==species,])
    # Replace space by underscore for file output
    species_string <- replace(species, " ", "_")
    
    ## ==== (log) Biomass ~ Dim_Z + Top_Area ====
    # Model
    lm_bmass_z_top_log_perspecies <- lm(log(biomass_g) ~ log(Dim_Z) + log(Top_Area), data=plant_perspecies)
    model_stats(lm_bmass_z_top_log_perspecies, plot_res=FALSE)
    # Observed vs predicted plot
    img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_Predict_Log_", species, "_", output_label, ".jpg")
    plot_obs_vs_pred(mdl=lm_bmass_z_top_log_perspecies, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
    
    ## ==== (log) Biomass ~ 60% Convex hull + Species ====
    # Model
    lm_bmass_chull_perspecies <- lm(log(biomass_g) ~ log(Convex_hull_60), data=plant_perspecies)
    model_stats(lm_bmass_chull_perspecies, plot_res=FALSE)
    # Observed vs predicted plot
    img_name <- paste0(output_folder, "//Biomass_ConvexHull60_Predict_Log_", species, "_", output_label, ".jpg")
    plot_obs_vs_pred(mdl=lm_bmass_chull_perspecies, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
    
    ## ==== (log) Biomass ~ Dim_Z + Top_Area + 60% Convex hull + Species ====
    # Model
    lm_bmass_z_top_chull_perspecies <- lm(log(biomass_g) ~ log(Dim_Z) + log(Top_Area) + log(Convex_hull_60), data=plant_perspecies)
    model_stats(lm_bmass_z_top_chull_perspecies, plot_res=FALSE)
    # Observed vs predicted plot
    img_name <- paste0(output_folder, "//Biomass_Heigth_TopArea_ConvexHull60_Predict_Log_", species, "_", output_label, ".jpg")
    plot_obs_vs_pred(mdl=lm_bmass_z_top_chull_perspecies, img_name=img_name, xlab="Predicted log(biomass(g))", ylab="Measured log(biomass(g))")
    
  }
}

# ---- Reset ----
# Reset graphic display
par(par_init)
