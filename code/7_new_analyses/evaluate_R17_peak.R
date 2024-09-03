### THIS SCRIPT EVALUATES SMH R17 PEAK TARGETS
## CALCULATES WIS FOR PEAK SIZE AND CRPS FOR PEAK TIMING
## COVERAGE FOR BOTH PEAK SIZE AND TIMING
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
# repo_data <- "../covid19-megaround_data/"
# folder_path <- "megaround-processed/"
# schema <- arrow::schema(
#   arrow::field("origin_date", arrow::string()),
#   arrow::field("scenario_id", arrow::string()),
#   arrow::field("location", arrow::string()),
#   arrow::field("target", arrow::string()),
#   arrow::field("horizon", double()),
#   arrow::field("type", arrow::string()),
#   arrow::field("type_id", double()),
#   arrow::field("model_name", arrow::string()),
#   arrow::field("value", double()),
# )

# read in scenario and model key
load(file.path(repo_data,"code/sysdata.rda")) 
scenario_key <- tibble::enframe(scenario2number, name = 'scenario', value = "scenario_id") %>% setDT()
model_key <- tibble::enframe(model2number, name = 'model', value = "model_name") %>% setDT()

# read in data
models_peak_time_seasons <- read_parquet("code/7_new_analyses/peak_targets/models_peak_hosp_seasons.parquet")

# read in ensemble data
time_frame_choice <- "seasons"
all_ens_files <- list.files("code/7_new_analyses/peak_targets/")
ens_files <- all_ens_files[grep("ensemble_LOP",all_ens_files)]
ens_files <- ens_files[grep(time_frame_choice, ens_files)]

ensembles <- rbind(do.call(rbind, lapply(ens_files, function(x) 
  read_parquet(paste0("code/7_new_analyses/peak_targets/",x)) %>% setDT() %>%
    .[, model_name := gsub(paste0("_peak_hosp_",time_frame_choice,".parquet"),"",x)])))

model_output <- read_parquet(paste0("code/7_new_analyses/peak_targets/models_peak_hosp_",time_frame_choice,".parquet")) %>% 
  setDT()

df <- rbindlist(list(ensembles, 
                     model_output %>%
                       .[, value := ifelse(target == "peak time hosp", quantile, value)] %>% 
                       .[, -c("quantile","target_end_date")]),use.names = TRUE)


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

# truth data (from fluview - check with Lucie)-----
gs_inc_hosp <- setDT(read.csv(paste0(repo_data, "visualization/data-goldstandard/hospitalization.csv"))) %>% 
  .[time_value > "2023-04-16"]
gs_inc_death <-  setDT(read.csv(paste0(repo_data,"visualization/data-goldstandard/fv_death_incidence_num.csv"))) %>%
  .[time_value > "2023-04-16"]
gold_standard_data <- rbindlist(l = list( copy(gs_inc_hosp) %>% 
                                            .[, target := "inc hosp"]
                                          # copy(gs_inc_death) %>% 
                                          #   .[, target := "inc death"], 
                                          # copy(gs_inc_hosp) %>% 
                                          #   .[, value := cumsum(value), 
                                          #     by = .(geo_value_fullname, fips)] %>% 
                                          #   .[, target := "cum hosp"], 
                                          # copy(gs_inc_death) %>% 
                                          #   .[, value := cumsum(value), 
                                          #     by = .(geo_value_fullname, fips)] %>% 
                                          #   .[, target := "cum death"])
)) %>%
  .[, time_value := as.IDate(time_value)] %>%
  rename(obs = "value") 

# start with just doing summer and winter (that's what we have data for so far)
gold_standard_peak <- gold_standard_data %>% 
  .[date_key_seasons, on = .(time_value = date)] %>%
  .[, .(`peak size hosp` = max(obs), 
        `peak time hosp` = as.character(time_value[which.max(obs)])), 
    by = .(geo_value_fullname, fips, target, time_frame)] %>% 
  melt(., measure.vars = c("peak time hosp","peak size hosp"), variable.name = "metric", value.name = "obs")  %>%
  na.omit()

# add observations to df_quantile
df_peak <- df %>% 
  .[gold_standard_peak %>% .[, -c("target")], on = .(time_frame, location = fips, target = metric), allow.cartesian = TRUE] 


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
cov_pksize <- cov[df_peak %>% .[target == "peak size hosp"] %>% 
                    .[, type_id := as.numeric(type_id)], 
                  on = .(quantile = type_id)] %>% 
  .[quantile != 0.5] %>% 
  .[, -c("type")] %>%
  .[, obs := as.numeric(obs)] %>%
  # reshape to make lwr and upr columnns
  data.table::dcast(model_name + target + scenario_id + geo_value_fullname + location + 
                      time_frame + alpha + obs ~ variable, value.var = "value")  %>%
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

# PEAK SIZE COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET ------ 
# peak value coverage
cov_pksize %>%
  na.omit() %>%
  .[, time_frame := factor(time_frame, levels = c("2023spring","2023summer","2023fall","2023winter","2024spring"))] %>%
  .[#substr(scenario_id,1,1) %in% c("E", "F") &
      !(model_name %in% c("NCSU-COVSIM", "ensemble_LOP_untrimmed", "ensemble"))] %>%
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

# ggsave("figures/new_analyses/R17_coverage_bymodel_peak_size_hosp.pdf", width = 10, height = 5)


# PEAK SIZE WIS -----------

WIS <- df_peak[target == "peak size hosp" & !is.na(type_id)] %>%
  .[, `:=` (quantile = as.numeric(type_id),
            obs = as.numeric(obs))] %>%
  .[, wis(quantile, value, obs, IS_components = TRUE),
    by = .(target, scenario_id, model_name, location, obs, geo_value_fullname, time_frame)] 

wis_tmp <- WIS %>%
  .[!model_name %in% c("NCSU-COVSIM", "ensemble_LOP_untrimmed")] %>%
  .[, time_frame := factor(time_frame, levels = c("2023spring","2023summer","2023fall","2023winter","2024spring"))] %>%
  # normalize WIS by round, target, location (do not group by scenario or model)
  .[, WIS_SD := sd(WIS, na.rm = TRUE), by = .(target, location, time_frame, geo_value_fullname)] %>%
  .[, WIS_norm := WIS/WIS_SD] %>%
  .[, .(WIS_norm = mean(WIS_norm, na.rm = TRUE)), by = .(target, model_name, scenario_id, time_frame)] %>%
  .[, scenario_letter := substr(scenario_id, 1,1)] 

ggplot(wis_tmp %>% .[scenario_letter == "F"]) + 
  geom_tile(aes(time_frame, model_name, fill = WIS_norm)) 

ggplot(wis_tmp %>% .[scenario_letter == "F"]) + 
  geom_line(aes(time_frame, WIS_norm, group = model_name, colour = model_name)) + theme_minimal()



# PEAK TIMING CRPS ----------
# CRPS(quantile,value,obs)
crps_scores_prep  <- df_peak %>%
  .[target == "peak time hosp"] %>%
  .[, year :=  as.numeric(substr(gsub("EW", "", type_id), 1, 4))] %>%
  .[, week := as.numeric(substr(gsub("EW", "",type_id), 5, 6))] %>%
  .[! geo_value_fullname %in% c("American Samoa", "Puerto Rico", "Virgin Islands")] %>%
  .[, target_end_date :=  MMWRweek2Date(year, week, MMWRday = 7)] %>%
  .[, horizon := as.numeric((as.Date(target_end_date) - as.Date(origin_date) + 1) / 7)] %>%
  # .[, obs2 := ifelse(lubridate::as_date(obs) > target_end_date, 0, 1)]
  .[, obs2 := as.numeric((as.Date(obs) - as.Date(origin_date) + 1) / 7)]
  
# for crps: quantile = value, value = horizon, obs = horizon of observation (?)
crps_scores <- crps_scores_prep %>%
  .[, .(CRPS = crps(value,horizon,obs2)),
    by=.(target, scenario_id, model_name, location, time_frame)] %>% # score across locations 
# by=.(target, scenario_id, model_name, location, time_frame)] 
  # normalize CRPS by round, target, location (do not group by scenario or model)
  .[, CRPS_SD := sd(CRPS, na.rm = TRUE), by = .(target, location, time_frame)] %>% # Q: what to normalize over 
  .[, CRPS_norm := CRPS/CRPS_SD] 

# A <-
  crps_scores %>% 
  .[scenario_id == "F-2023-04-16"] %>%
  .[, time_frame := factor(time_frame, levels = c("2023spring","2023summer","2023fall","2023winter","2024spring"))] %>%
  ggplot() + 
  geom_tile(aes(time_frame, model_name, fill = CRPS_norm)) + 
  theme_minimal()
B <- ggplot() + 
  geom_line(data = gold_standard_data %>% .[geo_value_fullname == "US" & target == "inc hosp"], aes(time_value, obs)) + 
  #add vertical lines to indicate seasons
  geom_vline(xintercept = as.Date(c("2023-03-21", "2023-06-21", "2023-09-21", "2023-12-21","2024-03-21")), linetype = "dashed") +
  theme_minimal() 
# cowplot::plot_grid(A,B, ncol = 1)
