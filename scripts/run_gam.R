# RUN GAM 

#### SET-UP ####################################################################

# > packages -------------------------------------------------------------------
packages <- c("arrow","beepr","dplyr","ebirdst","fs","glue","here","mgcv","sf",
              "terra","tictoc","tidyverse")
invisible(lapply(packages, library, character.only = TRUE))
rm(packages)


# > parameters and run options -------------------------------------------------
dml_label <- "hpc2" 
j <- 1
min_ss <- 30


# > directories ----------------------------------------------------------------
dir_data <- here("data","tidy")
dir_res <- here("results", dml_label)


# > load tidy data ------------------------------------------------------------------
vi <- arrow::read_parquet(here(dir_res, "vis.parquet")) %>% 
  filter(., model == "treatment", scenario == j) %>% 
  mutate(response = ifelse(response == "occurrence", "occ", response)); tail(vi)

srd_import <- readRDS(here(dir_data, "srd_tidy_2026-03-10.rds")) %>%  
  select(., c("srd_id", "year", "latitude", "longitude", unique(vi$predictor))) 


# load eBird relative abundance to join to final partial-effect dataframe
ebirdst_data_dir()
abd_rast <- load_raster("pinjay",
                        product = "abundance",
                        period = "seasonal",
                        metric = "mean",
                        resolution = "3km")

abd <- as.data.frame(abd_rast, xy = TRUE, cells = TRUE) %>% 
  rename(srd_id = cell, abd = resident)
head(abd)


#### GET TOP ENV PREDICTORS  ###################################################

# list of predictors to not consider for GAM analysis
nuisance <- c("eastness", "northness", "elevation", "elevation_sd",
              "northness_sd", "eastness_sd")


for (r in c("occ", "count")) {
  
  cat(r, "\n")
  
  top <- dplyr::filter(vi, response == r, !(predictor %in% nuisance)) %>%
    group_by(predictor) %>%
    summarize(med = median(importance_sc)) %>%
    filter(med > 0) %>%
    arrange(desc(med)) %>%
    head(15)
  
  pred_names <- top$predictor
  

#### PREPARE DATA & DIRECTORIES FOR GAM ANALYSIS ###############################

  experiment_label <- glue("{r}/sc{j}") #(experiment_label <- glue("full_{r}_2")) # comment this back on
  #experiment_label <- glue("{r}/test_public_repo")
  rfull <- if (r == "occ") { "occurrence" } else { "count"}
  
  srd <- select(srd_import, c("srd_id", "latitude", "longitude", pred_names)) 
  

  # > directories ----------------------------------------------------------------
  dir_int <- here("intermediates",dml_label, "gam", experiment_label)

  fs::dir_create(dir_int)


    # > correlations ---------------------------------------------------------------
    cor_mat <- srd %>% 
      select(., -c(srd_id, latitude, longitude)) %>% 
      cor(., method = "pearson") 
    
    corpiv <- pivot_longer(as_tibble(cor_mat) %>% mutate(var=row.names(cor_mat)), cols = 1:ncol(cor_mat)) %>% 
      mutate(value = as.numeric(value)) %>% 
      rename(var1 = var, var2 = name) %>%
      filter(var1 != var2, var1 != "fold", var2 != "fold") %>% 
      mutate(pair_id = paste(pmin(var1, var2), pmax(var1, var2), sep = "_")) %>%
      distinct(pair_id, .keep_all = TRUE) %>%
      select(-pair_id) %>% 
      arrange(desc(abs(value)))
    
    rm(cor_mat)
    head(corpiv)
    
    maxcor <- nrow(filter(corpiv, value >= 0.6))
    cat(r, "preds with corr > 0.6 =", maxcor, "\n")

    
    # > check for sufficient variation ---------------------------------------------
    pred_unique <- as.data.frame(NA)
    
    for (iii in 1:ncol(srd)){
      pred_unique[iii] <- length(unique(srd[ , names(srd)[iii]]))
    }
    colnames(pred_unique) <- colnames(srd)
    pred_unique
    
    cat(r, "preds with insuff var =",
    pred_unique %>%
      select(where(~ any(.x < min_ss))) %>%
      ncol(),
    "\n"
    )


#### RUN GAM ###################################################################

  model_name <- NULL
  for (iii in 1:length(pred_names)){                      # for each predictor
    model_name[iii] <- 
      paste0( "s( ",pred_names[iii],",bs=\"ds\", k=4, m=c(1,0)) + ")         
  }      
  model_name <- c("cate_pred ~", model_name,
                  "s(longitude, latitude, bs=\"ds\",k=200, m=c(1,0.5)) +",  
                  "1")
  (model_name <- paste(model_name, collapse = " "))
  
  rm(corpiv, top, iii, maxcor) 
  
  # set up objects to receive outputs
  gam.summaries <- as.data.frame(matrix(nrow = 0, ncol = length(pred_names)+3))
  colnames(gam.summaries) <- c("fold", "metric", pred_names, "latlon")
  r.sq <- NULL
  
  gc()
  
  # > begin for loops ------------------------------------------------------------
  
  folds <- c(1:100)
  
  cates <- arrow::read_parquet(here(dir_res, glue("{r}_cates.parquet"))) %>% 
    filter(scenario == j, fold %in% folds)
  
  for (i in unique(cates$fold)) {
    
    cat("fold ", i, "\n")
    
    tic()
    cates_i <- filter(cates, fold == i)
    
    # > join cates and srd -------------------------------------------------------
    D <- inner_join(select(cates_i, c(fold, srd_id, cate_pred)), srd, by = "srd_id")
    head(D)
    
    colSums(is.na(D))[colSums(is.na(D)) > 0]
    
    # run the model. note selection and extra gamma penalization
    d.gam <- mgcv::bam( 
      as.formula(model_name), 
      gamma = 5,                       # increase above 1 to force smoother fits
      select = T,                      # should selection penalties be added to the smooth effects, so that they can in principle be penalized out of the model? 
      discrete = T, 
      data = D)
    
    #save(d.gam, file = here(dir_int, "model.RData"))
    
    
    # > summarize ------------------------------------------------------------------
    
    # F and p value from each fold
    d.sum <- summary(d.gam)
    d.add <- as.data.frame(d.sum$s.table) %>% select(-c(Ref.df)) %>%
      rownames_to_column("predictor") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      mutate(fold = rep(i), .before = metric) %>% 
      slice(-1)
    names(d.add) <- names(gam.summaries)                   # assign column names
    gam.summaries <- rbind(gam.summaries, d.add)           # add gam results for fold i to full gam.summaries table
    
    # R2 and deviance explained from each fold
    r.sq.i <- cbind(fold = i, r.sq = d.sum$r.sq, dev.expl = d.sum$dev.expl * 100)
    r.sq <- rbind(r.sq, r.sq.i)
    
    gc()
    
    
    #### PREDICT INDIVIDUAL FEATURE EFFECTS ###########################################################
    
    pred  <- predict.gam(d.gam,type = "iterms", se.fit = TRUE) # iterms gives the parts of the prediction ("How much is s(x1) contributing here?"); type = "response" gives total prediction (sum of all contributions).
    
    d.effect_plotting <- cbind(data.frame(srd_id = D$srd_id), data.frame(pred)) # marginal effects and SE for this fold
    d.effect_plotting <- d.effect_plotting %>% mutate(fold = i, .before = srd_id)
    
    if(i == min(folds)) {
      effect_plotting <- d.effect_plotting
      #rm(d.effect_plotting)
    } else{
      effect_plotting <- rbind(effect_plotting, d.effect_plotting)
    }
    
    rm(d.effect_plotting, pred, d.gam, d.add)
    gc()
    toc()
    
  } # {end i loop} ----

beep(1)



#### SAVE RESULTS ##############################################################

saveRDS(effect_plotting, here(dir_int, glue("effect_plotting_folds{min(folds)}-{max(folds)}.rds")))
saveRDS(gam.summaries, here(dir_int, glue("gam-summaries_folds{min(folds)}-{max(folds)}.rds")))
saveRDS(r.sq, here(dir_int, glue("r-sq_folds{min(folds)}-{max(folds)}.rds")))


rm(r.sq.i, D, d.sum, cates)
gc()


#} # {end b loop} ----

rm(gam.summaries, r.sq)


#### PREPARE MEAN DATAFRAME ####################################################

cat("prepare mean dataframe\n")

effect_1D_mean <- effect_plotting %>% 
  group_by(srd_id) %>%
  summarize(across(starts_with("fit"), list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}"))
#head(effect_1D_mean)


# > add rel abund from eBird Status 2023 ---------------------------------------

effect_1D_mean <- left_join(effect_1D_mean, abd, by = "srd_id")


# > save mean dataframe --------------------------------------------------------

saveRDS(effect_1D_mean, here(dir_int, "effect_1D_mean.rds"))

rm(effect_1D_mean, effect_plotting)
gc()


} # end r loop




#### END OF CODE FOR MANUSCRIPT ################################################

#### FOR CREATING THE MEAN EFFECT DF ONLY ######################################

  # > load fold-wise gam results -------------------------------------------------
  
  (gam_files <- list.files(dir_int, pattern = "effect_plotting"))
  
  all_gam <- readRDS(here(dir_int, "effect_plotting_folds1-100.rds"))
  
  beep(1)
  head(all_gam)
  rm(gam_files, g)
  
  
  # > prepare mean dataframe -----------------------------------------------------
  
  effect_1D_mean <- all_gam %>% 
    group_by(srd_id) %>%
    summarize(across(starts_with("fit"), list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"))
  head(effect_1D_mean)
  
  
  # > add rel abund from eBird Status 2023 ---------------------------------------
  ebirdst_data_dir()
  abd_rast <- load_raster("pinjay",
                          product = "abundance",
                          period = "seasonal",
                          metric = "mean",
                          resolution = "3km")
  
  abd <- as.data.frame(abd_rast, xy = TRUE, cells = TRUE) %>% 
    rename(srd_id = cell, abd = resident)
  head(abd)
  
  effect_1D_mean <- left_join(effect_1D_mean, abd, by = "srd_id")
  
  rm(all_gam, abd, abd_rast)
  gc()


saveRDS(effect_1D_mean, here(dir_int, "effect_1D_mean.rds"))

