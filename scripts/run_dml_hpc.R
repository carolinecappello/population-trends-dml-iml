
# slurm variables
i <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
n_threads=64
i=1

experiment_label <- "pinjay_hpc_2"
file_label <- ""
data_notes <- "hpc run"


#### LOAD ######################################################################

# > packages -------------------------------------------------------------------
packages <- c("ebirdst", "fields", "fs", "glue", "grf", "lubridate", "jsonlite",
              "ranger", "rnaturalearth", "sf", "terra", "tictoc", "tidyverse")
invisible(lapply(packages, library, character.only = TRUE))
rm(packages)


# > directories and scripts ----------------------------------------------------
source("scripts/model_scenarios.R")
source("R/partial-effects_cdc.R")

message("Processing task ", i, " with ", n_threads, " cores.")
options(Ncpus = n_threads, ranger.num.threads = n_threads)


# > data -----------------------------------------------------------------------
erd_tidy <- readRDS("data/tidy/erd_tidy_2026-03-10.rds")
srd_tidy <- readRDS("data/tidy/srd_tidy_2026-03-10.rds")

# srd template
srd_template <- rast("data/srd_3km_raster_template.tif")

# state boundaries for mapping -- are these used for the trend map below? If not, remove
states <- rnaturalearth::ne_states(iso_a2 = c("CA", "US", "MX")) |>
  st_geometry() |>
  st_transform(st_crs(srd_template))
states_sf <- st_sf(geometry = states)



#### SET UP ####################################################################

# > global details -------------------------------------------------------------

season_date_range <- c("2023-02-01", "2023-7-01") |> ymd()
start_year <- 2012L # baseline year
end_year <- 2023L

# model parameters
n_trees_rf <- 2000
n_trees_cf <- 2000
use_simulated_data <- FALSE
save_intermediates <- TRUE

propensity_model <- "regression"
propensity_package <- "ranger"
model_notes <- "importance = impurity"

n_folds <- 100
n_scenarios <- 2
set.seed(1) # does not change
(seeds <- sample(1:2000, 100, replace = F)) # 1017 1860  679  129  930...


# > select features ------------------------------------------------------------

for (j in 1:n_scenarios) {   # {begin scenario loop, j} ----

  tic()
  cat(">>> scenario ", j, "<<<\n")

  sub_label <- paste0("sc",j)
  fig_dir <- path("figures", experiment_label, sub_label)
  dir_create(fig_dir)
  int_dir <- path("intermediates", experiment_label, sub_label)
  dir_create(int_dir)

  propensity_features <- get(glue("propensity_features_{j}"))
  marginal_features <- get(glue("marginal_features_{j}"))
  treatment_features <- get(glue("treatment_features_{j}"))

  all_variables <- unique(c("checklist_id", "obs_count", "observer_id",
                            "latitude", "longitude",
                            "day_of_year", "year",
                            propensity_features,
                            marginal_features,
                            treatment_features))

  all_srd_variables <- unique(c("srd_id", "latitude", "longitude", "year",
                                propensity_features,
                                marginal_features,
                                treatment_features))

  erd <- erd_tidy[, all_variables]
  srd <- srd_tidy


  # > train and test datasets ----------------------------------------------------

  occ_summ <- vector(mode = "list", length = 9)
  count_summ <- vector(mode = "list", length = 9)
  vi <- NULL

  #for (i in 1:n_folds){  # {begin replicate loop, i} ---- # if running locally, not on the hpc

  #cat(">> fold ", i, "<<\n")


  # > spatiotemporal sampling ----------------------------------------------------

  # decide which dataset to use for training
  if (use_simulated_data) {
    train_data <- sim
  } else {
    train_data <- erd  
  }

  # grid-sampled training data for occurrence model
  set.seed(seeds[i])
  train_occ_all <- grid_sample_stratified(train_data,                               # sample points from xy, using a grid (raster) as stratification
                                          coords = c("longitude", "latitude", "day_of_year"),
                                          #obs_column = "obs_count",                # toggle on if case_control = TRUE
                                          case_control = FALSE,                     # whether presence and absence should be sampled independently
                                          sample_size_per_cell = 1,                 # default
                                          cell_sample_prop = 1,                     # proportion ⁠(0-1]⁠; if less than 1, only this proportion of cells will be randomly selected for sampling
                                          by_year = TRUE)                           # default


  # grid-sampled training data for count model
  set.seed(seeds[i])
  train_cnt_all <- train_data[train_data$obs_count > 0, ] %>%                   # filter to rows that have observations
    grid_sample_stratified(.,                                                   # sample points from xy, using a grid (raster) as stratification
                           coords = c("longitude", "latitude", "day_of_year"),
                           case_control = FALSE,                                # whether presence and absence should be sampled independently
                           sample_size_per_cell = 1,                            # up to n points are sampled from each stratum (cell). 1=default
                           cell_sample_prop = 1,                                # proportion ⁠(0-1]⁠; if less than 1, only this proportion of cells will be randomly selected for sampling
                           by_year = TRUE)                                      # default

  train_occ <- train_occ_all
  train_cnt <- train_cnt_all

  message("Training data size: ", nrow(train_data))
  message("Occurrence sample size: ", nrow(train_occ))
  message("Count sample size: ", nrow(train_cnt))


  # > srd prediction unit --------------------------------------------------------
  # set prediction unit for srd to median of sampled training data (for cci and effort features)

  pred_unit <- summarize(train_occ_all,
                         across(all_of(c(effort_features, "cci")), median))

  srd_pred_unit <- srd |>                                                       # combine pred_unit with srd df
    filter(year == end_year) |>
    bind_cols(pred_unit)
  srd_pred_unit <- srd_pred_unit[,all_srd_variables]

  rm(pred_unit)
  gc()


  #### OCCURRENCE MODEL ################################################################

  # 1) propensity score ----------------------------------------------------------

  cat("OCCURRENCE MODEL\n")
  cat("propensity: ")

  W <- train_occ$year
  X <- train_occ[, propensity_features]

  # probability forest for observation year
  if (propensity_package == "ranger") {
    cat("ranger, ")
    if (propensity_model == "classification") {
      cat("classification\n")
      W_factor <- factor(W, levels = seq(end_year, start_year))                     # convert years to factor format
      occ_propensity <- ranger(x = X,                                               # predictor data (propensity features)
                               y = W_factor,                                        # response vector (year)
                               num.trees = n_trees_rf,
                               probability = TRUE,                                  # grow a probability forest (Malley et al 2012)
                               importance = "impurity",                             # variable importance mode; impurity = Gini index for classification, the variance of the responses for regression
                               num.threads = n_threads)                             # for parallelization. NULL uses number of cores available on computer. 1 disables parallelization

      # probability each checklist was assigned to each treatment (i.e. year)
      W_hat_occ_probs <- occ_propensity$predictions                                 # each row sums to 1. each row is a sampled ebird checklist
      # most likely treatment (i.e. year) that a checklist was assigned to
      W_hat_occ <- colnames(W_hat_occ_probs)[max.col(W_hat_occ_probs)] |>
        as.integer()

      # variable importance
      vi <- tibble(response = "occurrence", model = "propensity",
                   predictor = names(occ_propensity$variable.importance),
                   importance = unname(occ_propensity$variable.importance)) |>
        bind_rows(vi)

    } else if (propensity_model == "regression") {
      cat("regression\n")
      occ_propensity <- ranger(x = X,                                               # regression tree
                               y = W,
                               num.trees = n_trees_rf,
                               importance = "impurity",
                               num.threads = n_threads)
      #  predict year of each checklist (continuous)
      W_hat_occ <- occ_propensity$predictions

      # variable importance
      vi <- tibble(fold = i, response = "occurrence", model = "propensity",
                   predictor = names(occ_propensity$variable.importance),
                   importance = unname(occ_propensity$variable.importance),
                   importance_sc = occ_propensity$variable.importance / sum(occ_propensity$variable.importance)) |>
        bind_rows(vi)
    }

  } else if (propensity_package == "grf") {
    cat("grf\n")
    occ_propensity <- grf::regression_forest(X = X,
                                             Y = W,
                                             num.trees = n_trees_rf,
                                             tune.parameters = "all",
                                             num.threads = n_threads)
    #  predict year of each checklist (continuous)
    W_hat_occ <- stats::predict(occ_propensity)$predictions
    # variable importance
    vi <- tibble(response = "occurrence", model = "propensity",
                 predictor = names(occ_propensity$X.orig),
                 importance = variable_importance(occ_propensity)[, 1]) |>
      bind_rows(vi)

  }

  gc()


  # 2) marginal occurrence -------------------------------------------------------
  #probability of detection marginalized across years
  cat("marginal\n")

  Y_occ <- as.numeric(train_occ[["obs_count"]] > 0)                          # binary vector, was the focal species observed in this checklist?
  Y_occ_fct <- factor(Y_occ, levels = c(0, 1))
  X <- train_occ[, marginal_features]                                        # predictor variables (observer details, weather details, unchanging landscape features)

  # observed proportion and number of detections across all checklists
  #prop_detections <- mean(Y_occ)
  n_detections <- sum(Y_occ)

  # check if factors would have been coded incorrectly
  if (Y_occ[1] == 1) {
    warning("Error in factor coding!! j = ", j)
  }
  # this is essentially the status base occurrence model without year
  occ_marginal <- ranger(x = X,                                                 # predictor variables (observer details, weather details, unchanging landscape features)
                         y = Y_occ_fct,                                         # response; observed yes or no
                         num.trees = n_trees_rf,
                         probability = TRUE,
                         #sample.fraction = c(prop_detections, prop_detections), # fraction of observations to sample. For classification, this can be a vector of class-specific values. So, we're sampling equal numbers of checklists with and without detections.
                         min.node.size = ceiling(n_detections * 0.01),           # minimal node size to split at (i.e., how many checklists with focal-species detection needed).
                         importance = "impurity",
                         num.threads = n_threads)
  Y_hat_occ <- occ_marginal$predictions[, "1", drop = TRUE]                     # predicted probability of each checklist having an occurrence
  Y_hat_occ <- pmax(0, pmin(1, Y_hat_occ))                                      # bounds the predictions between 0 and 1. Values below 0 become 0. Values above 1 become 1


  # > predict marginal to srd ----------------------------------------------------
  # given, the features of interest (i.e., marginal features) in each srd, what is probability of species occurence in each srd?
  occ_pred <- predict(occ_marginal,
                      data = srd_pred_unit[, marginal_features],
                      num.threads = n_threads)                                  # number of CPUs
  occ_pred <- occ_pred$predictions[, "1", drop = TRUE]                          # save predictions. drop = TRUE creates a simplified output (i.e., vector instead of single-column matrix)
  occ_pred <- pmax(0, pmin(1, occ_pred))                                        # bounds the predictions between 0 and 1
  # variable importance
  vi <- tibble(fold = i, response = "occurrence", model = "marginal",
               predictor = names(occ_marginal$variable.importance),
               importance = unname(occ_marginal$variable.importance),
               importance_sc = occ_marginal$variable.importance / sum(occ_marginal$variable.importance)) |>
    bind_rows(vi)


  # 3) occurrence treatment effects, i.e. occurrence trend -------------------------------
  # treatment is continuous here, so this model will estimate the conditional
  # average treatment effect per unit treatment assignment (W)
  # in this case treatment is year, so the model is estimating the difference in occ per year
  cat("treatment\n")

  X <- train_occ[, treatment_features]
  occ_treatment <- grf::causal_forest(X,                                        # pland baseline & change and elev/topo/water/road features
                                      Y_occ,                                    # binary vector, was the focal species observed in this checklist?
                                      W,                                        # vector of true checklist years (treatment)
                                      W.hat = W_hat_occ,                        # estimates of the treatment propensities E[W | Xi].from our propensity model: the most likely treatment (i.e. year) that a checklist was assigned to
                                      Y.hat = Y_hat_occ,                        # estimates of the expected responses E[Y | Xi], marginalizing over treatment.from our marginal model: predicted probability of each checklist having an occurrence
                                      num.trees = n_trees_cf,
                                      tune.parameters = "all",                  # vector of parameter names to tune. "all" = all tunable params are tuned by cross-validation: ("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves", "alpha", "imbalance.penalty")
                                      num.threads = n_threads)
  # predict conditional average treatment effects for the srd
  occ_cate <- predict(occ_treatment,
                      newdata = srd_pred_unit[, treatment_features],
                      num.threads = n_threads)
  occ_cate <- occ_cate$predictions
  # variable importance: how often a variable 𝑋𝑗was split on.
  vi <- tibble(fold = i, response = "occurrence", model = "treatment",
               predictor = names(occ_treatment$X.orig),
               importance = variable_importance(occ_treatment)[, 1],
               importance_sc = variable_importance(occ_treatment)[, 1] / sum(variable_importance(occ_treatment)[, 1])) |>
    bind_rows(vi)
  gc()


  #### COUNT MODEL #########################################################

  cat("COUNT MODEL\n")

  # 1) propensity score ----
  cat("propensity: ")

  # code is the same as occurrence
  W <- train_cnt$year
  X <- train_cnt[, propensity_features]

  # probability forest for observation year
  if (propensity_package == "ranger") {
    cat("ranger, ")

    if (propensity_model == "classification") {
      cat("classification\n")
      W_factor <- factor(W, levels = seq(end_year, start_year))
      count_propensity <- ranger(x = X,
                                 y = W_factor,
                                 num.trees = n_trees_rf,
                                 probability = TRUE,
                                 importance = "impurity",
                                 num.threads = n_threads)
      # the probability each checklists was assigned to each treatment (i.e. year)
      W_hat_count_probs <- count_propensity$predictions
      # the most likely treatment (i.e. year) that a checklists was assigned to
      W_hat_count <- colnames(W_hat_count_probs)[max.col(W_hat_count_probs)] |>
        as.integer()

      # variable importance
      vi <- tibble(response = "count", model = "propensity",
                   predictor = names(count_propensity$variable.importance),
                   importance = unname(count_propensity$variable.importance)) |>
        bind_rows(vi)


    } else if (propensity_model == "regression") {
      cat("regression\n")
      count_propensity <- ranger(x = X,
                                 y = W,
                                 num.trees = n_trees_rf,
                                 probability = FALSE,
                                 importance = "impurity",
                                 num.threads = n_threads)
      #  predict year of each checklist (continuous)
      W_hat_count <- count_propensity$predictions

      # variable importance
      vi <- tibble(fold = i, response = "count", model = "propensity",
                   predictor = names(count_propensity$variable.importance),
                   importance = unname(count_propensity$variable.importance),
                   importance_sc = count_propensity$variable.importance / sum(count_propensity$variable.importance)) |>
        bind_rows(vi)
    }

  } else if (propensity_package == "grf") {
    cat("grf\n")
    count_propenisty <- grf::regression_forest(X, W,
                                               num.trees = n_trees_rf,
                                               tune.parameters = "all",
                                               num.threads = n_threads)
    #  predict year of each checklist (continuous)
    W_hat_count <- stats::predict(count_propenisty)$predictions
    # variable importance
    vi <- tibble(response = "count", model = "propensity",
                 predictor = names(count_propenisty$X.orig),
                 importance = variable_importance(count_propenisty)[, 1]) |>
      bind_rows(vi)

  }


  # 2) marginal count, estimated count marginalized across years ----
  cat("marginal\n")

  Y_log <- log(train_cnt$obs_count + 1)                                         # add 1 and take the log. add 1 because log(1) = 0 and we can't have a 0 count because count is given occurrence = 1
  X <- train_cnt[, marginal_features]
  count_marginal <- ranger(x = X,                                               # regression tree
                           y = Y_log,
                           num.trees = n_trees_rf,
                           importance = "impurity",
                           #keep.inbag = TRUE,
                           num.threads = n_threads)

  Y_hat_count <- count_marginal$predictions                                     # save predictions
  # predict to srd
  count_pred <- predict(count_marginal,
                        data = srd_pred_unit[, marginal_features],
                        num.threads = n_threads)
  count_pred <- count_pred$predictions


  # variable importance
  vi <- tibble(fold = i, response = "count", model = "marginal",
               predictor = names(count_marginal$variable.importance),
               importance = unname(count_marginal$variable.importance),
               importance_sc = count_marginal$variable.importance / sum(count_marginal$variable.importance)) |>
    bind_rows(vi)

  # 3) count treatment effects, i.e. count trend ----
  # treatment is continuous here, so this model will estimate the conditional
  # average treatment effect per unit treatment assignment (W)
  # in this case treatment is year, so the model is estimating the difference in
  # log(count + 1) per year
  cat("treatment\n")

  X <- train_cnt[, treatment_features]
  count_treatment <- causal_forest(X,
                                   Y_log,
                                   W,
                                   W.hat = W_hat_count,
                                   Y.hat = Y_hat_count,
                                   tune.parameters = "all",
                                   num.trees = n_trees_cf,
                                   num.threads = n_threads)
  # predict conditional average treatment effects for the srd
  count_cate <- predict(count_treatment,                                        # dataframe n_srd x 1
                        newdata = srd_pred_unit[, treatment_features],
                        num.threads = n_threads)
  count_cate <- count_cate$predictions                                          # convert to vector
  # variable importance
  vi <- tibble(fold = i, response = "count", model = "treatment",
               predictor = names(count_treatment$X.orig),
               importance = variable_importance(count_treatment)[, 1],
               importance_sc = variable_importance(count_treatment)[, 1] / sum(variable_importance(count_treatment)[, 1])) |> # same as importance. Looks like grf scales these for you.
    bind_rows(vi)


  #### OUTPUTS ###################################################################

  # > cates etc. -------------------------------------------------------------
  get_cates <- function(treatment_mod, train, test, W_hat, summ_list, response) {

    cat("cates, ")
    predict_cates <- predict(treatment_mod,
                             newdata = get(glue("srd_pred_unit"))[, treatment_features],
                             estimate.variance = TRUE,
                             num.threads = n_threads,
                             drop = TRUE)                                       # default: output simplified to vector, not kept as a matrix

    cate_df <-  as.data.frame(predict_cates) %>%
      dplyr::mutate(fold = i, srd_id = get(glue("srd_pred_unit"))$srd_id,
                    response = "cate",
                    std_err = sqrt(variance.estimates),
                    cate_upper = predictions + 1.2816 * std_err,                # ~80% two-sided CI, z-score cut-off (critical value) = 1.2816
                    cate_lower = predictions - 1.2816 * std_err,                # ~80% two-sided CI
                    signif = sign(cate_upper) == sign(cate_lower)) %>%
      rename(cate_pred = predictions, variance = variance.estimates) %>%
      relocate(fold, srd_id, response)
    #head(cate_df)


    cat("ates, ")
    ate_df <- grf::average_treatment_effect(treatment_mod) %>%
      {data.frame(fold = i, ate = .[1], ate_std_err = .[2])}
    rownames(ate_df) <- NULL


    cat("predict propensity, ")
    prop_pred_df <- as.data.frame(cbind(checklist_id = train$checklist_id,
                                        prop_pred = W_hat)) %>%
      mutate(fold = i, response = "prop_pred") %>%
      relocate(fold, response)
    #head(prop_pred_df)


    cat("predict marginal, ")
    marg_pred_df <- as.data.frame(cbind(checklist_id = train$checklist_id,
                                        pred = get(paste0("Y_hat_", response)))) %>%
      mutate(fold = i, response = "marg_pred") %>%
      relocate(fold, response)
    #head(marg_pred_df)


    marg_srd_pred_df <- as.data.frame(cbind(srd_id = get(glue("srd_pred_unit"))$srd_id,
                                            marg_srd_pred = get(paste0(response,"_pred")))) %>%
      mutate(fold = i, response = "marg_srd_pred") %>%
      relocate(fold, response)
    #head(marg_srd_pred_df)


    cat("test calibration, ")
    treatment_calibration <- grf::test_calibration(treatment_mod)
    treatment_calibration <- as.data.frame(rbind(NA,treatment_calibration)[-1,]) %>%
      rownames_to_column(var = "term") %>%
      mutate(fold = i) %>%
      relocate(fold)
    colnames(treatment_calibration) <- c("fold", "term", "estimate", "std_err", "t", "p")
    # head(treatment_calibration)


    cat(" prop performance\n")
    prop_mod <- get(glue("{response}_propensity"))

    # propensity training data
    train_mse_baseline <- mean((train$year - mean(train$year))^2)


    perform_prop <- as.data.frame(cbind(fold = i,
                                        pred_error = prop_mod$prediction.error,
                                        mse_baseline = train_mse_baseline,
                                        r_sq = prop_mod$r.squared,
                                        n_checklists = nrow(train))) 
    #head(perform_prop)


    cat(" prop marg performance\n")
    marg_mod <- get(glue("{response}_marginal"))
    if (response == "count") {
      baseline_mse <- mean((Y_log - mean(Y_log)) ^ 2)
    } else {
      baseline_mse <- NULL
    }

    perform_marg <- as.data.frame(cbind(fold = i,
                                        pred_error = marg_mod$prediction.error, # brier score for prob and mse for regression trees
                                        mse_baseline = baseline_mse,            # null for probability trees
                                        r_sq = marg_mod$r.squared))             # null for probability trees
    #head(perform_marg)



    summ_list[[1]] <- rbind(summ_list[[1]], cate_df)
    summ_list[[2]] <- rbind(summ_list[[2]], ate_df)
    summ_list[[3]] <- rbind(summ_list[[3]], prop_pred_df)
    summ_list[[4]] <- rbind(summ_list[[4]], marg_pred_df)
    summ_list[[5]] <- rbind(summ_list[[5]], marg_srd_pred_df)
    summ_list[[6]] <- rbind(summ_list[[6]], perform_prop)
    summ_list[[7]] <- rbind(summ_list[[7]], perform_marg)
    summ_list[[8]] <- rbind(summ_list[[8]], treatment_calibration)

    cat("dev explained\n")
    for (m in c("propensity", "marginal", "treatment")) {
      if (response == "occ") {
        t <- train_occ
      } else {
        t <- train_cnt
      }

      partials_gam <- partial_effects(data = t, model = get(glue("{response}_{m}")))
      options(scipen = 999)
      r_squared <- round(summary(partials_gam)$r.sq,4)
      dev_expl <- round(summary(partials_gam)$dev.expl,4)
      dev_explained_mr <- as.data.frame(cbind(fold = i, model = glue("{response}_{m}"),
                                              dev_explained = dev_expl, adj_r2 = r_squared))
      summ_list[[9]] <- rbind(summ_list[[9]], dev_explained_mr)
    }

    names(summ_list) <- c("cates", "ates", "prop_pred", "marg_pred", "marg_srd",
                          "perf_prop", "perf_marg", "grf_checks", "dev_explained")
    return(summ_list)
  }


  occ_summ <- get_cates(treatment_mod = occ_treatment, summ_list = occ_summ,
                        train = train_occ, test = test_occ, W_hat = W_hat_occ,
                        response = "occ")
  count_summ <- get_cates(treatment_mod = count_treatment, summ_list = count_summ,
                          train = train_cnt, test = test_cnt, W_hat = W_hat_count,
                          response = "count")


  # > trend maps ------------------------------------------------------------------------

  cat("maps\n")

  occ_ppy <- occ_cate / occ_pred                                                # relative change per year                
  count_ppy <- count_cate


  # occ_pred: given X (i.e., marginal features), what is the prob of species occurrence in each srd?
  srd_pred <- data.frame(srd_id = srd_pred_unit$srd_id,
                         srd_yr = srd_pred_unit$year,
                         occ_pred = occ_pred,                                   
                         # back-transform count from log scale
                         count = exp(count_pred) - 1,                           
                         abd = occ_pred * count_pred,
                         occ_ppy = occ_ppy,
                         # chain rule
                         abd_ppy = occ_ppy + count_cate,
                         # count is already a rate because it was log-transformed
                         count_ppy = count_ppy
  ) |>
    mutate(across(everything(), ~ ifelse(is.finite(.x), .x, NA_real_)))

  ppy <- rast(srd_template, nlyr = 2, names = c("occ_ppy", "count_ppy"))

  for (l in seq_along(names(ppy))) {
    set.values(ppy, cells = srd_pred$srd_id,
               values = srd_pred[[names(ppy)[l]]],
               layer = l)
  }

  ppy3 <- trim(ppy)

  # trend breaks and palette
  breaks <- seq(-0.1, 0.1, by = 0.02)  
  breaks_legend <- breaks
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  pal <- ebirdst_palettes(n = length(breaks) - 1, type = "trends")


  km <- 3
  # maps
  for (n in names(ppy)) {
    n_name <- case_match(n,
                         "occ_ppy" ~ "Occurrence",
                         "count_ppy" ~ "Count",
                         "occ_ppy_sig" ~ "Significant Occurrence",
                         "count_ppy_sig" ~ "Significant Count")
    png(path(fig_dir, glue("map_trend_{n}{file_label}_{i}_sc{j}.png")),
        width = 1600, height = 1600)
    par(mar = c(0.5, 0.5, 4, 10))
    # basemap
    plot(states, col = "#cfcfcf", lwd = 3,
         extent = get(glue("ppy{km}")),
         main = glue("{n_name} Trend ({start_year}-{end_year}), fold {i}"),
         cex.main = 3)
    # data
    plot(get(glue("ppy{km}"))[[n]],
         breaks = breaks, col = pal,
         legend = FALSE, add = TRUE)
    # state lines
    plot(states, lwd = 3, add = TRUE)
    # legend
    image.plot(zlim = range(breaks_legend),
               breaks = breaks_legend, col = pal,
               smallplot = c(0.92, 0.94, 0.15, 0.85),
               axis.args = list(at = breaks_legend,
                                labels = scales::percent(breaks_legend, 1),
                                cex.axis = 2.5),
               legend.only = TRUE)
    box(lwd = 3)
    dev.off()
  }


  # > collect and save all the models --------------------------------------------
  models <- list(
    occurrence = list(
      propensity = occ_propensity,
      marginal = occ_marginal,
      treatment = occ_treatment
    ),
    count = list(
      propensity = count_propensity,
      marginal = count_marginal,
      treatment = count_treatment
    )
  )
  if (save_intermediates) saveRDS(models, path(int_dir, glue("models_{i}_sc{j}.rds")))

  gc()

  #} # {end replicates i} ----                                                  # if running locally, not on hpc

  (time_elapsed <- toc()$callback_msg)
  saveRDS(occ_summ, path(int_dir, glue("occ_summ_{i}_sc{j}.rds")))
  saveRDS(count_summ, path(int_dir, glue("count_summ_{i}_sc{j}.rds")))
  saveRDS(vi, path(int_dir, glue("vi_{i}_sc{j}.rds")))

  # > model notes ----------------------------------------------------------------

  models_info <- list(
    experiment_label = experiment_label,
    scenario = j,
    run_date = Sys.Date(),
    season_date_range = season_date_range,
    start_year = start_year,
    end_year = end_year,
    approx_n_train_occ = nrow(train_occ_all),
    approx_n_train_cnt = nrow(train_cnt_all),
    n_trees_rf = n_trees_rf,
    n_trees_cf = n_trees_cf,
    n_threads = n_threads,
    folds = n_folds,
    seed = seeds,
    time_elapsed = time_elapsed,
    propensity_model = propensity_model,
    propensity_package = propensity_package,
    propensity_features = propensity_features,
    marginal_features = marginal_features,
    treatment_features = treatment_features,
    data_notes = data_notes,
    model_notes = model_notes
  )

  sink(path(int_dir, glue("model_info_sc{j}.txt")))
  print(models_info)
  sink()

} # {end scenarios j} ----





