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


df_quantile <- dplyr::filter(dc, type == "quantile" & (target == "inc death" |  target == "inc hosp")) %>% 
  dplyr::collect()
setDT(df_quantile)
df_quantile <- df_quantile %>%
  .[, origin_date := as.IDate(origin_date)]

# truth data (from fluview - check with Lucie)
gs_inc_hosp <- read.csv(paste0(repo_data, "visualization/data-goldstandard/hospitalization.csv"))
gs_inc_death <- read.csv(paste0(repo_data,"visualization/data-goldstandard/fv_death_incidence_num.csv"))
gold_standard_data <- rbindlist(l = list(setDT(gs_inc_hosp) %>% 
                                           .[, target := "inc hosp"], 
                                         setDT(gs_inc_death) %>% 
                                           .[, target := "inc death"])) %>%
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

#### PLOT RESULTS --------------------------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# COVERAGE ACROSS LOCATIONS AND WEEKS FOR EACH SCENARIO/TARGET
cov %>%
  .[!is.na(cov) &
      substr(scenario_id,1,1) %in% c("E", "F") &
      model_name != "NCSU-COVSIM"] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, model_name, scenario_id)] %>%
  ggplot(aes(x = alpha, y = cov_summ, color = model_name)) +
  geom_abline(linetype = "dashed") +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), rows = vars(target)) + 
  labs(x = "expected coverage", y = "actual coverage") +
  scale_color_manual(values = c("grey70", "grey45", "black", gg_color_hue(9))) +
  scale_x_continuous(expand = c(0,0),
                     label = percent) +
  scale_y_continuous(expand = c(0,0),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())
ggsave("figures/new_analyses/R17_coverage_bymodel_overall.pdf", width = 8, height = 8)

# 50% AND 95% COVERAGE ACROSS WEEKS FOR ALL MODELS/TARGETS
cov %>%
  .[!is.na(cov) & 
      model_name != "NCSU-COVSIM" &
      substr(scenario_id,1,1) %in% c("E", "F") &
      alpha %in% c(0.5, 0.95)] %>%
  .[, .(cov_summ = sum(cov)/.N), 
    by = .(alpha, origin_date, target, time_value, model_name, scenario_id)] %>%
  ggplot(aes(x = time_value, y = cov_summ, color = model_name)) + 
  geom_hline(aes(yintercept = alpha), size = 1) +
  geom_line() + 
  #geom_point(size = 1) + 
  facet_grid(cols = vars(scenario_id), rows = vars(alpha, target)) + 
  labs(x = "projection week", y = "actual coverage") +
  scale_color_manual(values = c("grey70", "grey45", "black", gg_color_hue(9))) +
  scale_linetype_manual(values = c("dashed", "twodash", "dotdash", rep("solid",9))) +
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),
                     label = percent) +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        panel.grid = element_blank())
ggsave("figures/new_analyses/R17_coverage_bymodel_byweek.pdf", width = 8, height = 8)

