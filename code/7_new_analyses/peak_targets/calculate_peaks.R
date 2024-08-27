# CALCULATES PEAK TARGETS FROM TRAJECTORIES, AND THEN CREATES ENSEMBLE 

# Load libraries
library(arrow)
library(tidyverse)
library(data.table)
library(MMWRweek)
library(RColorBrewer)
library(scales)
library(CombineDistributions)


#### LOAD THE DATA -------------------------------------------------------------
repo_data <- "../covid19-megaround_data/"
folder_path <- "megaround-processed/"
schema <- arrow::schema(
  arrow::field("origin_date", arrow::string()),
  arrow::field("scenario_id", arrow::string()),
  arrow::field("location", arrow::string()),
  arrow::field("target", arrow::string()),
  arrow::field("horizon", double()),
  arrow::field("type", arrow::string()),
  arrow::field("type_id", double()),
  arrow::field("model_name", arrow::string()),
  arrow::field("value", double()),
)

# read in scenario and model key
load(file.path(repo_data,"code/sysdata.rda")) 
scenario_key <- tibble::enframe(scenario2number, name = 'scenario', value = "scenario_id") %>% setDT()
model_key <- tibble::enframe(model2number, name = 'model', value = "model_name") %>% setDT()

# read in data
dc <- arrow::open_dataset(paste0(repo_data, folder_path), 
                          partitioning = "model_name", 
                          format = "parquet", 
                          schema = schema, 
                          factory_options = list(
                            exclude_invalid_files = TRUE))

df_sample <- dplyr::filter(dc, type == "sample", 
                           # location %in% c("06","01"), model_name%in% c("Ensemble","JHU_IDD-CovidSP","UVA-EpiHiper"),
                           target %in% c("inc hosp", "inc death"),
                           horizon <= 52*2) %>%
  dplyr::collect() %>%
  setDT() 

df_sample <- df_sample %>%
  .[, origin_date := lubridate::as_date(origin_date)] %>%
  .[, time_value :=  origin_date + horizon*7-1]

# team_models <- unique(df_sample$model_name)
# for testing
team_models <- c("UVA-adaptive", "MOBS_NEU-GLEAM_COVID")

# make a key where dates are indexed by quarter year/season
date_key_seasons <- data.table::data.table(date = unique(df_sample$time_value)) %>% 
  .[, `:=` (quarter = lubridate::quarter(date, fiscal_start = 3),
            month = lubridate::month(date),
            year = lubridate::year(date))] %>%
  .[, season := ifelse(quarter == 1, "spring",
                       ifelse(quarter == 2, "summer",
                              ifelse(quarter == 3, "fall", "winter")))] %>%
  .[, `:=` (time_frame = ifelse(season == "winter" & month %in% c(1,2), paste0(year - 1, season),paste0(year, season)))] 
#S1 = spring, S2 = summer, S3 = fall, S4 = winter

# make a yearly date key (april to april)
date_key_yearly <- data.table::data.table(date = lubridate::as_date(unique(df_sample$time_value))) %>% 
  .[, `:=` (time_frame = ifelse(date <= lubridate::as_date("2024-04-27"), "2023-24", "2024-25"))] 

# PROBABILITY OF PEAKS ----------
probs <- c(0.01, 0.025, seq(0.05,0.95, 0.05), 0.975, 0.99)

## ~ for whole time period -------
# Calculate peak size 
models_peak_targets <- rbind(do.call(rbind,lapply(team_models, function(i) {
  calculate_peak_size(
    df_sample %>% .[model_name == i], df_all = NULL, i,
    peak_size_target = "peak size hosp",
    peak_group = c("origin_date", "scenario_id", "location", "target",
                   "type_id"),
    quantile_vect = probs
  ) %>% mutate(model_name = i)
})
) %>% prep_peak_size(., probs),
do.call(rbind,lapply(team_models, function(i) {
  calculate_peak_time(
    df_sample %>% .[model_name == i], df_all = NULL, i,
    peak_time_target = "peak time hosp",
    peak_group = c("origin_date", "scenario_id", "location", "target",
                   "type_id")
  ) %>% mutate(model_name = i)
})
) %>% prep_peak_time()
)
## ~ for yearly period -------

models_peak_yearly <- rbind(do.call(rbind,lapply(team_models, function(i) {
  res <- NULL
  for(j in unique(date_key_yearly$time_frame)){
    res_tmp <- calculate_peak_size(
      df_sample %>% 
        .[date_key_yearly, on = .(time_value = date)] %>%
        .[model_name == i & time_frame == j], 
      df_all = NULL, 
      i,
      peak_size_target = "peak size hosp",
      peak_group = c("origin_date", "scenario_id", "location", "target",
                     "type_id"),
      quantile_vect = probs
    ) %>% mutate(model_name = i, time_frame = j)
    res <- rbind(res, res_tmp)
  }
  res
})
) %>% prep_peak_size(., probs),
do.call(rbind,lapply(team_models, function(i) {
  res <- NULL
  for(j in unique(date_key_yearly$time_frame)){
    res_tmp <- calculate_peak_time(
      df_sample %>% 
        .[date_key_yearly, on = .(time_value = date)] %>%
        .[model_name == i & time_frame == j], 
      df_all = NULL, i,
      peak_time_target = "peak time hosp",
      peak_group = c("origin_date", "scenario_id", "location", "target",
                     "type_id")
    ) %>% mutate(model_name = i, time_frame = j)
    res <- rbind(res, res_tmp)
  }
  res
})
) %>% prep_peak_time()
)


## ~ for seasonal period -------

models_peak_seasons <- rbind(do.call(rbind,lapply(team_models, function(i) {
  res <- NULL
  for(j in unique(date_key_seasons$time_frame)){
    res_tmp <- calculate_peak_size(
      df_sample %>% 
        .[date_key_seasons, on = .(time_value = date)] %>%
        .[model_name == i & time_frame == j], 
      df_all = NULL, 
      i,
      peak_size_target = "peak size hosp",
      peak_group = c("origin_date", "scenario_id", "location", "target",
                     "type_id"),
      quantile_vect = probs
    ) %>% mutate(model_name = i, time_frame = j)
    res <- rbind(res, res_tmp)
  }
  res
})
) %>% prep_peak_size(., probs),
do.call(rbind,lapply(team_models, function(i) {
  res <- NULL
  for(j in unique(date_key_seasons$time_frame)){
    res_tmp <- calculate_peak_time(
      df_sample %>% 
        .[date_key_seasons, on = .(time_value = date)] %>%
        .[model_name == i & time_frame == j], 
      df_all = NULL, i,
      peak_time_target = "peak time hosp",
      peak_group = c("origin_date", "scenario_id", "location", "target",
                     "type_id")
    ) %>% mutate(model_name = i, time_frame = j)
    res <- rbind(res, res_tmp)
  }
  res
})
) %>% prep_peak_time()
)

# CALCULATE ENSEMBLES -----------

weight_scheme <- "equal"

# loop through and do all ensembles
for(ens in c("LOP_untrimmed","LOP_trimmed")){
  for(type in c("yearly","seasons")){
    ensemble <- peak_ensemble(
      df_peak = get(paste0("models_peak_",type)),
      weighting_scheme = weight_scheme,
      n_trim = ifelse(ens == "LOP_untrimmed", NA, 2), # trimmed ensemble
      weight_df = NULL,
      ens_group = c("origin_date","scenario_id","location","target","horizon","time_frame"),
      id_var = "model_name"
    )
    arrow::write_parquet(ensemble, 
                         paste0("code/7_new_analyses/peak_targets/ensemble_",ens,"_peak_hosp_",type,".parquet"))
  }
}


# # ADD GROUNDTRUTH ------
# 
# # add peak timing gold standard 
# gold_standard_peak <- gold_standard_data %>%
#   .[date_key, on = .(time_value = date)] %>% na.omit() %>% .[str_detect(target, "inc")] %>%
#   .[, .(peak_value = max(obs), 
#         peak_timing = as.character(time_value[which.max(obs)])), 
#     by = .(geo_value_fullname, fips, target, time_frame)] %>%  
#   .[, peak_timing := as.numeric((as.Date(peak_timing) - as.Date("2023-04-16") + 1) / 7)] %>% # horizon
#   melt(., measure.vars = c("peak_timing","peak_value"), variable.name = "metric", value.name = "obs") 
