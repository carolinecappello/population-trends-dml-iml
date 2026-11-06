# FINAL MODEL SCNENARIOS

# create lists of features to go in the propensity, conditional mean, and treatment models


#### LIST FEATURES #############################################################

# observation process features
effort_features <- c("day_of_year",
                     "solar_noon_diff",                                      
                     "effort_hrs", "effort_distance_km",
                     "effort_speed_kmph", "is_stationary",
                     "num_observers")

# elevation and topography features
elevation_features <- c("eastness",     
                        "eastness_sd",  
                        "northness",     
                        "northness_sd",  
                        "elevation",     
                        "elevation_sd")  

# habitat features
static_features <- c("astwbd_c2_pland",   # lake class
                     "astwbd_c3_pland",   # river class
                     "road_density_c1",   # highways  
                     "road_density_c2",   # primary 
                     "road_density_c3",   # secondary
                     "road_density_c4",   # tertiary
                     "road_density_c5")   # local
       

static_features_combined <- c("astwbd_pland",   # lakes and rivers COMBINED  
                              "road_density")   # road length per unit of area COMBINED (km of road per km^2 of area)

eco_regions <- c("eco_1_6", 
                 "eco_1_9", 
                 "eco_1_10", 
                 "eco_1_11", 
                 "eco_1_13") #"eco_1_12"

yearly_features_modis <- c("mcd12q1_lccs1_c1_pland",      # barren
                           "mcd12q1_lccs1_c2_pland",      # perm. snow and ice
                           "mcd12q1_lccs2_c25_pland",     # forest/cropland mosaic -- often an ecotone between herbaceous cropland and forest
                           "mcd12q1_lccs2_c36_pland",     # herbaceous cropland -- wheat,corn,barley,soy,alfalfa,hay,potatoes,root veg // planted & harvested seasonally or annually // no woody stems // land is managed
                           "mcd12q1_lccs3_c27_pland",     # woody wetlands 
                           "mcd12q1_lccs3_c50_pland",     # herbaceous wetlands
                           "mcd12q1_lccs3_c51_pland")     # tundra 

yearly_modis_reduced <- c("mcd12q1_lccs1_c1_pland",      # barren
                          "mcd12q1_lccs1_c2_pland",      # perm. snow and ice
                          "mcd12q1_lccs2_c25_pland",     # forest/cropland mosaic -- often an ecotone between herbaceous cropland and forest
                          "mcd12q1_lccs2_c36_pland",     # herbaceous cropland -- wheat,corn,barley,soy,alfalfa,hay,potatoes,root veg // planted & harvested seasonally or annually // no woody stems // land is managed
                          "mcd12q1_wetlands_pland")      # wetlands (woody & herbaceous combined)

yearly_features <- c("ntl",            # nighttime lights
                     "ppt",            # precipitation
                     "tmax",           # temperature
                     "rap_annual",     # annual grasses and forbs      
                     "rap_perennial",  # perennial grasses and forbs
                     "rap_shrub",      # shrubs
                     "rap_pj",         # pinyon-juniper pland
                     "rap_nonpj",      # non-piyon-juniper tree
                     "fire")           # fire in the last 10yr

land_config <- c("pj_contagion", "pj_ed") # pj woodland contagion, pinyon-juniper edge

disturbance <- c("disturb_clearcut_all",
                 "disturb_insect_all",         
                 "disturb_mechanical_all", # all other mechanical removal        
                 "disturb_thinning_all")
                 #"disturb_removal" # combines clearcut, thinning, other mech


#### SORT FEATURES #############################################################

# > model scenario 1 -----------------------------------------------------------
# "reduced var model, with updated cci and contagion"

propensity_features_1 <- c(effort_features, 
                        "cci",
                        elevation_features,
                        static_features,
                        eco_regions,
                        paste0(yearly_features_modis, "_baseline"),
                        paste0(yearly_features, "_baseline"),
                        paste0(disturbance, "_baseline"),
                        paste0(land_config[1], "_baseline")
                        )

marginal_features_1 <- c(effort_features, 
                       "cci",
                       elevation_features,
                       static_features, 
                       eco_regions,
                       paste0(yearly_features_modis, "_baseline"),
                       paste0(yearly_features, "_baseline"),
                       paste0(disturbance, "_baseline"),
                       paste0(land_config[1], "_baseline")
                       )

treatment_features_1 <- c(elevation_features, 
                          static_features_combined, 
                          eco_regions, 
                          paste0(yearly_modis_reduced, "_baseline"),
                          paste0(yearly_modis_reduced, "_change"),
                          paste0(yearly_features, "_baseline"), 
                          paste0(yearly_features, "_change"),
                          paste0(disturbance, "_baseline"), 
                          paste0(disturbance, "_change"),
                          paste0(land_config[1], "_baseline"),
                          paste0(land_config[1], "_change"),
                          "pdsi")

