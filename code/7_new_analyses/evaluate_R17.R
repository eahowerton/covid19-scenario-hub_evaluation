### THIS SCRIPT EVALUATES SMH R17 PROJECTIONS
library(arrow)
library(tidyverse)
library(data.table)
library(MMWRweek)
library(RColorBrewer)
library(scales)


#### SCORING FUNCTIONS ---------------------------------------------------------
source("code/3_score_projections/scoring_functions.R")


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

dc <- arrow::open_dataset(paste0(repo_data, folder_path), 
                          partitioning = "model_name", 
                          format = "parquet", 
                          schema = schema, 
                          factory_options = list(
                            exclude_invalid_files = TRUE))


df_quantile <- dplyr::filter(dc, type == "quantile") %>%  #& (target == "inc death" |  target == "inc hosp")
  dplyr::collect()
setDT(df_quantile)
df_quantile <- df_quantile %>%
  .[, origin_date := as.IDate(origin_date)]

# truth data (from fluview - check with Lucie)
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


# add observations to df_quantile
df_quantile <- df_quantile %>% 
  .[, time_value :=  origin_date + horizon*7-1] %>%
  .[gold_standard_data, on = .(time_value, target, location = fips)] %>%
  .[!is.na(origin_date)]


#### CALCULATE COVERAGE --------------------------------------------------------
cov <- data.table(alpha = c(seq(0.1, 0.9, 0.1), 0.95, 0.98))
# find upper and lower intervals for all alpha levels
cov$upr <- cov$alpha/2 + 0.5
cov$lwr <- 1-(cov$alpha/2 + 0.5)
cov <- melt(cov, "alpha", value.name = "quantile")
cov$quantile = round(cov$quantile,3)
setDT(cov)

# merge with proj to assign alpha and lwr/upr to each quantile in proj
cov <- cov[df_quantile, on = .(quantile = type_id)] %>% 
  .[quantile != 0.5] %>% 
  # reshape to make lwr and upr columnns
  data.table::dcast(origin_date + scenario_id + location + target + horizon + time_value  + 
                      model_name + alpha + obs ~ variable, value.var = "value")  %>%
  .[, ":=" (cov = ifelse(obs < upr & obs > lwr, 1, 0))] %>% 
  .[, ":=" (upr = NULL,
            lwr = NULL, 
            obs = NULL)]

#### CALCULATE WIS 00000--------------------------------------------------------

# calculate WIS for all projections
wis <-  df_quantile[, wis(type_id,value,obs,IS_components = TRUE),
                by=.(target, origin_date, time_value, scenario_id, model_name, location, geo_value_fullname, obs)]

#### PLOT RESULTS --------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET
cov %>%
  .[!is.na(cov) &
      #substr(scenario_id,1,1) %in% c("E", "F") &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), rows = vars(target)) + 
  labs(x = "expected coverage", y = "actual coverage") +
  scale_color_manual(values = c("black", gg_color_hue(8))) +
  scale_x_continuous(expand = c(0,0),
                     label = percent) +
  scale_y_continuous(expand = c(0,0),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())

ggsave("figures/new_analyses/R17_coverage_bymodel_overall.pdf", width = 10, height = 5)

# 50% AND 95% COVERAGE ACROSS WEEKS FOR ALL MODELS/TARGETS
cov %>%
  .[!is.na(cov) & 
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble")) &
      #substr(scenario_id,1,1) %in% c("E", "F") &
      alpha %in% c(0.5)] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, time_value, model_name, scenario_id)] %>%
  ggplot(aes(x = time_value, y = cov_summ, color = model_name)) + 
  geom_hline(aes(yintercept = alpha), linetype = "dotted") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(rows = vars(scenario_id), cols = vars(target)) + 
  labs(x = "projection week", y = "actual coverage") +
  scale_color_manual(values = c("black", gg_color_hue(8))) +
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())

ggsave("figures/new_analyses/R17_coverage_bymodel_byweek.pdf", width = 8, height = 8)

## QQ PLOT FOR CUMULATIVE TARGETS IN FINAL WEEK
cov %>%
  .[, max_time_value := max(time_value), by = .(target)] %>%
  .[!is.na(cov) &
      time_value == max_time_value &
      #substr(scenario_id,1,1) %in% c("E", "F") &
      substr(target, 1,3) == "cum" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), rows = vars(target)) + 
  labs(x = "expected coverage", y = "actual coverage") +
  scale_color_manual(values = c("black", gg_color_hue(8))) +
  scale_x_continuous(expand = c(0,0),
                     label = percent) +
  scale_y_continuous(expand = c(0,0),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())
ggsave("figures/new_analyses/R17_coverage_bymodel_cumulative.pdf", width = 8, height = 4)

# WIS OVER TIME BY MODEL AND LOCATION
wis %>%
  .[!is.na(WIS) & 
      target == "inc hosp" &
      substr(scenario_id,1,1) %in% c("E", "F") &
      geo_value_fullname == "Florida" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  ggplot(aes(x = time_value, y = WIS/obs, color = model_name)) + 
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(rows = vars(model_name), cols = vars(scenario_id), scales = "free") + 
  labs(x = "projection week", y = "actual coverage") +
  scale_color_manual(values = c("black", gg_color_hue(8))) +
  scale_x_date(expand = c(0,0)) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())

wis %>%
  .[!is.na(WIS) & 
      target == "inc hosp" &
      geo_value_fullname == "Florida" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, quarter := lubridate::quarter(time_value)] %>%
  ggplot(aes(x = obs, y = WIS, color = as.factor(quarter))) + 
  geom_point() + 
  geom_path() + 
  facet_grid(rows = vars(scenario_id), cols = vars(model_name), scales = "free") + 
  scale_color_brewer(palette = "Set1") + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())

us_latlong <- read.csv("code/7_new_analyses/US_GeoCode.csv")

p1 = gold_standard_data %>%
  .[us_latlong[, c("latitude", "longitude", "Name")], on = .(geo_value_fullname = Name)] %>%
  .[!is.na(time_value)] %>%
  .[, relobs := obs/max(obs), by = c("geo_value_fullname", "target")] %>%
  .[target == "inc hosp"] %>%
  ggplot(aes(x = time_value, y = relobs, color = latitude, group = geo_value_fullname)) + 
  geom_line(alpha = 0.75) + 
  scale_color_distiller(palette = "Blues", direction = 0) + 
  scale_x_date(expand = c(0,0)) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
  
p2 = wis %>%
  .[!is.na(WIS) & 
      target == "inc hosp" &
      substr(scenario_id,1,1) %in% c("E", "F") &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, quarter := lubridate::quarter(time_value)] %>%
  .[, relWIS_rank := rank(WIS/obs), by = .(model_name, location, geo_value_fullname, scenario_id)] %>%
  .[us_latlong[, c("latitude", "longitude", "Name")], on = .(geo_value_fullname = Name)] %>%
  .[!is.na(target)] %>%
  .[, geo_value_fullname := factor(geo_value_fullname)] %>%
  .[, geo_value_fullname := reorder(geo_value_fullname, latitude)] %>%
  ggplot(aes(x = time_value, y = geo_value_fullname, fill = relWIS_rank)) + 
  geom_tile() + 
  facet_grid(cols = vars(scenario_id), rows = vars(model_name), 
             scales = "free") + 
  labs(y = "location (sorted by latitude)") + 
  scale_fill_distiller(palette = "YlOrBr", direction = 0) + 
  scale_x_date(expand = c(0,0)) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none", 
        legend.title = element_blank(),
        panel.grid = element_blank(), 
        strip.background = element_blank())

cowplot::plot_grid(
  cowplot::plot_grid(p1, p1 + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank()), nrow = 1), 
  p2, 
  ncol = 1, rel_heights = c(0.2, 0.8))

ggsave("figures/new_analyses/R17_WIS_bymodel_bylocation.pdf", width = 8, height = 10)

df_quantile %>% 
  .[location == "US" &
      model_name == "Ensemble_LOP" &
      type_id %in% c(0.05, 0.25, 0.5, 0.75, 0.95)] %>%
  .[, quantile := paste0("Q", type_id*100)] %>%
  data.table::dcast(time_value + scenario_id + target ~ quantile, value.var = "value") %>%
  ggplot(aes(x = time_value)) + 
  geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.2) + 
  geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.4) + 
  geom_line(aes(y = Q50)) + 
  geom_point(data = gold_standard_data[fips == "US" & time_value > "2023-04-16"], 
             aes(x = time_value, y = obs)) + 
  facet_grid(rows = vars(target), 
             cols = vars(scenario_id), scales = "free") + 
  theme_bw()
 

df_quantile %>% 
  .[, max_time_value := max(time_value), by = .(target)] %>%
  .[model_name == "Ensemble_LOP" &
      target == "cum hosp" &
      substr(target,1,3) == "cum" & 
      type_id %in% c(0.05, 0.25, 0.5, 0.75, 0.95) & 
      time_value == max_time_value] %>% 
  .[, quantile := paste0("Q", type_id*100)] %>%
  data.table::dcast(time_value + scenario_id + target + geo_value_fullname + obs ~ quantile, value.var = "value") %>%
  .[, scenario_letter := factor(substr(scenario_id,1,1), 
                                levels = rev(LETTERS[1:6]))] %>% 
  ggplot(aes(y = scenario_letter)) + 
  geom_segment(aes(x = Q5, xend = Q95, 
                   yend = scenario_letter), size = 0.2) + 
  geom_segment(aes(x = Q25, xend = Q75, 
                   yend = scenario_letter), size = 0.6) + 
  geom_point(aes(x = Q50), size = 1) +
  geom_text(aes(x = Q50, label = scenario_letter), color = 'white', size = 1) + 
  geom_vline(aes(xintercept = obs), color = "red") + 
  facet_wrap(vars(geo_value_fullname), scales = "free") + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor = element_blank())


### MAIN TEXT FIG OPTIONS ------------------------------------------------------
scenario_labs <- c("Low immune escape", "High immune escape")
names(scenario_labs) <- paste(c("E", "F"), "2023-04-16", sep = "-")

p1 <- cov %>%
  .[!is.na(cov) & 
      target == "inc hosp" &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble")) &
      substr(scenario_id,1,1) %in% c("E", "F") &
      alpha %in% c(0.5, 0.9)] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, time_value, model_name, scenario_id)] %>%
  .[, model_name := factor(model_name, levels = c("JHU_IDD-CovidSP",
                                                  "MOBS_NEU-GLEAM_COVID" , 
                                                  "NotreDame-FRED",
                                                  "UNCC-hierbin",
                                                  "USC-SIkJalpha",
                                                  "UTA-ImmunoSEIRS",
                                                  "UVA-adaptive",
                                                  "UVA-EpiHiper", 
                                                  "Ensemble_LOP"))] %>%
  ggplot(aes(x = time_value, y = cov_summ, 
             alpha = model_name, color = model_name)) + 
  geom_hline(aes(yintercept = alpha), linetype = "dotted") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), 
             rows = vars(alpha), 
             labeller = labeller(scenario_id = scenario_labs), 
             switch = "y") + 
  labs(x = "projection week", y = "actual coverage") +
  scale_alpha_manual(values = c(rep(0.6, 8), 1)) +
  scale_color_manual(values = c(rep("gray", 8), "black")) +
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())


### CHECKS ---
cov %>%
  .[!is.na(cov) &
      #substr(scenario_id,1,1) %in% c("E", "F") &
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble"))] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>% 
  .[substr(scenario_id,1,1) == "E" & 
      alpha == 0.5 & 
      model_name == "Ensemble_LOP" & 
      target == "inc hosp"]

cov %>%
  .[!is.na(cov) & 
      !(model_name %in% c("NCSU-COVSIM", "Ensemble_LOP_untrimmed", "Ensemble")) &
      #substr(scenario_id,1,1) %in% c("E", "F") &
      alpha %in% c(0.5)] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, time_value, model_name, scenario_id)] %>%
  .[substr(scenario_id,1,1) == "E" & 
      alpha == 0.5 & 
      model_name == "Ensemble_LOP" & 
      target == "inc hosp"] %>%
  .[, .(t = mean(cov_summ))]

