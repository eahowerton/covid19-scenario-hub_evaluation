### THIS SCRIPT EVALUATES SMH R17 PEAK TARGETS
library(arrow)
library(tidyverse)
library(data.table)
library(MMWRweek)
library(RColorBrewer)
library(scales)


#### SCORING FUNCTIONS ---------------------------------------------------------
source("code/3_score_projections/scoring_functions.R")

crps <- function(q,v,o,size=1000, rule=2, ties="ordered") {
  x = approx(q,v,runif(size),rule=rule, ties=ties)
  mean(abs(x$y-o[1])) - 1/2*mean(abs(x$y-approx(q,v,runif(size),rule=rule, ties=ties)$y))
}


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

df_sample <- dplyr::filter(dc, type == "sample", target %in% c("inc death", "inc hosp")) %>%
  dplyr::collect()
setDT(df_sample)
df_sample <- df_sample %>%
  .[, origin_date := as.IDate(origin_date)]

# can redo getting the peaks, or just read in the files that Lucie put together to get the peaks
# death, hosp, for 3 time frames
peak_files <- list.files(file.path(repo_data,"visualization/data-visualization/models/round17/peak")) %>% 
  .[grepl("peak_dist", .)] 
df_peak <- sprintf(file.path(repo_data,"visualization/data-visualization/models/round17/peak",
                         peak_files)) %>%
  readr::read_csv() %>%
  setDT() %>%
  .[metric %in% c("peak_value","peak_timing")] %>%
  .[method == "max"] %>% # just do max method for now, since observation data is easy to get max?
  .[model_key, on = .(model_name)] %>%
  .[scenario_key, on = .(scenario_id)] %>%
  na.omit() %>%
  .[,-c("scenario_id", "model_name")] %>%
  rename(model_name = model, scenario_id = scenario)

# truth data (from fluview - check with Lucie)-----
gs_inc_hosp <- setDT(read.csv(paste0(repo_data, "visualization/data-goldstandard/hospitalization.csv"))) %>% 
  .[time_value > "2023-04-16"]
gs_inc_death <-  setDT(read.csv(paste0(repo_data,"visualization/data-goldstandard/fv_death_incidence_num.csv"))) %>%
  .[time_value > "2023-04-16"]
gold_standard_data <- rbindlist(l = list( copy(gs_inc_hosp) %>% 
                                            .[, target := "inc hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, target := "inc death"], 
                                          copy(gs_inc_hosp) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum hosp"], 
                                          copy(gs_inc_death) %>% 
                                            .[, value := cumsum(value), 
                                              by = .(geo_value_fullname, fips)] %>% 
                                            .[, target := "cum death"])) %>%
  .[, time_value := as.IDate(time_value)] %>%
  rename(obs = "value") 
# start with just doing summer and winter (that's what we have data for so far)
gold_standard_peak <- rbindlist(list(gold_standard_data %>% 
                                         .[time_value >= "2023-04-01" & time_value <= "2023-09-30"] %>% 
                                         .[str_detect(target, "inc")] %>%
                                         .[, time_frame := "summer 2023"],
                                     gold_standard_data %>% 
                                         .[time_value >= "2023-10-01" & time_value <= "2024-03-31"] %>% 
                                         .[str_detect(target, "inc")] %>%
                                         .[, time_frame := "winter 2023"],
                                     gold_standard_data %>% 
                                         .[time_value >= "2023-10-01" & time_value <= "2024-09-30"] %>% 
                                         .[str_detect(target, "inc")] %>%
                                         .[, time_frame := "year 2023"])) %>%
  .[, .(peak_value = max(obs), 
        peak_timing = as.character(time_value[which.max(obs)])), 
    by = .(geo_value_fullname, fips, target, time_frame)] %>% 
  melt(., measure.vars = c("peak_timing","peak_value"), variable.name = "metric", value.name = "obs") 

# add observations to df_quantile
df_peak <- df_peak %>% 
  .[gold_standard_peak, on = .(time_frame, target, location_name = geo_value_fullname, metric)] 


#### CALCULATE COVERAGE --------------------------------------------------------
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
# find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile")
cov$quantile = round(cov$quantile,3)
setDT(cov)

### Peak size ------ 
# merge with proj to assign alpha and lwr/upr to each quantile in proj
cov_pksize <- cov[df_peak %>% .[metric == "peak_value"], on = .(quantile = type_id)] %>% 
  .[quantile != 0.5] %>% 
  .[, -c("type")] %>%
  na.omit() %>%
  # reshape to make lwr and upr columnns
  data.table::dcast(model_name + metric + scenario_id + location_name + target + 
                      method + time_frame + fips + alpha + obs ~ variable, value.var = "value")  %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%  # does this work for date and value? probs not
  .[, ":=" (upr = NULL,
            lwr = NULL, 
            obs = NULL)] %>%
  .[, names(.) := lapply(.SD, setattr, "sorted", NULL)] 

### Peak timing ------
cov_pktime <- cov[df_peak %>% .[metric == "peak_timing"], on = .(quantile = type_id)] %>% 
  .[quantile != 0.5] %>% 
  .[, -c("type")] %>%
  na.omit() %>%
  .[, `:=` (value = lubridate::as_date(value),
            obs = lubridate::as_date(obs))] %>%
  # reshape to make lwr and upr columnns
  data.table::dcast(model_name + metric + scenario_id + location_name + target + 
                      method + time_frame + fips + alpha + obs ~ variable, value.var = "value")  %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>%  # does this work for date and value? probs not
  .[, ":=" (upr = NULL,
            lwr = NULL, 
            obs = NULL)] %>%
  .[, names(.) := lapply(.SD, setattr, "sorted", NULL)] 

#### PLOT RESULTS --------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET ------ 
# peak value coverage
cov_pksize %>%
  na.omit() %>%
  .[target == "inc death" &
    #substr(scenario_id,1,1) %in% c("E", "F") &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, target, model_name, scenario_id, time_frame)] %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), rows = vars(time_frame)) + 
  labs(x = "expected coverage", y = "actual coverage") +
  # scale_color_manual(values = c("black", gg_color_hue(8))) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1),
                     label = percent) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank()) 

ggsave("figures/new_analyses/R17_coverage_bymodel_peak_size_death.pdf", width = 10, height = 5)

cov_pktime %>%
  na.omit() %>%
  .[target == "inc hosp" &
      #substr(scenario_id,1,1) %in% c("E", "F") &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, target, model_name, scenario_id, time_frame)] %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), rows = vars(time_frame)) + 
  labs(x = "expected coverage", y = "actual coverage") +
  # scale_color_manual(values = c("black", gg_color_hue(8))) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1),
                     label = percent) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank()) 
ggsave("figures/new_analyses/R17_coverage_bymodel_peak_timing_hosp.pdf", width = 10, height = 5)



# PROBABILITY OF PEAKS ----------

# From the sample values, calculate the peak (just max)
df_sample_test <- dplyr::filter(dc, type == "sample", location == "US", model_name == "JHU_IDD-CovidSP",
                                target %in% c("inc hosp", "inc death"),
                                horizon <= 52*2) %>%
  dplyr::collect() %>%
  setDT() 

df_sample_test <- df_sample_test %>%
  .[, origin_date := lubridate::as_date(origin_date)] %>%
  .[, time_value :=  origin_date + horizon*7-1]

# make a key where dates are indexed by quarter year/season
date_key <- data.table::data.table(date = unique(df_sample_test$time_value)) %>% 
  .[, `:=` (quarter = lubridate::quarter(date, fiscal_start = 3),
            year = lubridate::year(date))] %>%
  .[, season := ifelse(quarter == 1, "spring",
                       ifelse(quarter == 2, "summer",
                              ifelse(quarter == 3, "fall", "winter")))] %>%
  .[, `:=` (time_frame = paste0(year, season))] 
#S1 = spring, S2 = summer, S3 = fall, S4 = winter
  
# same as above calc just taking max
# calculates peak timing and peak value by model, target, scenario and time frame
df_sample_test2 <- df_sample_test %>% 
  .[date_key, on = .(time_value = date)] %>%
  .[, .(peak_value = max(value), 
        peak_timing = as.character(time_value[which.max(value)])), 
    by = .(scenario_id, location, target, type_id, model_name, time_frame)] %>% 
  melt(., measure.vars = c("peak_timing","peak_value"), variable.name = "metric", value.name = "obs") 

# gives the probability of peak in a given week by 
# model, target, scenario and time frame
df_sample_test3 <- df_sample_test %>%
  .[date_key, on = .(time_value = date)] %>%
  .[, `:=` (peak_value = max(value), 
        peak_timing = lubridate::as_date((time_value[which.max(value)]))), 
    by = .(scenario_id, location, target, type_id, model_name, time_frame)] %>%
  .[, is_peak := ifelse(time_value == peak_timing, 1, 0)] %>%
  .[, .(prob_peak = sum(is_peak)/.N), by = .(scenario_id, location, target, model_name, time_value, time_frame)] 

# CRPS for probability of peak timing

