# Functions for calculating peak targets 


#' Calculate missing peak size
#'
#' For rounds with samples output_type, the peak size hosp target is optional.
#' This function test for if peak size hosp is present in a specific team_model
#' submission and if not generate them
#'
#' @param df_team model projection of a team-model (as submitted)
#' @param df_all model projection of a team-model (with samples and quantiles)
#' to append the results, NULL for no data to append
#' @param team_model name of the team-model
#' @param peak_size_target target name as in the submission file
#' @param peak_group character vector of column names to 'group_by' the input
#' data frame to calculate the peak size hospitalization target for each group
#' @param quantile_vect numeric vector of probabilities to produce quantiles
#' corresponding to the given probabilities
#'
#' @details
#' Generates missing target from incident hospitalization trajectories by
#' extracting the max incident hospitalization value for each trajectories and
#' then generates the quantiles distribution of the max valua across all
#' trajectories.
#'
#' @export
#' @importFrom dplyr filter group_by across all_of mutate ungroup bind_rows
#' @importFrom dplyr reframe summarise
calculate_peak_size <- function(
    df_team, df_all, team_model, peak_size_target = "peak size hosp",
    peak_group = c("origin_date", "scenario_id", "location", "target",
                   "type_id"),
    quantile_vect = c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
                      0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
                      0.95, 0.975, 0.99)) {
  
  df_size_target <- dplyr::filter(df_team, target == peak_size_target)
  if (nrow(df_size_target) < 1) {
    df_inc_sample <- dplyr::filter(df_team, type == "sample",
                                   target == "inc hosp")
    if (nrow(df_inc_sample) < 1) {
      warning("Missing inc hosp sample to calculate peak size hosp.")
      df_size <- NULL
    } else {
      message("The submission file from: ", team_model,
              " is missing value for the target: ", peak_size_target,
              ". The missing value will be calculated.")
      df_size <- dplyr::group_by(df_inc_sample, across(all_of(peak_group))) %>%
        dplyr::summarise(max = max(value)) %>%
        dplyr::ungroup() %>%
        dplyr::reframe(value = unname(quantile(max, quantile_vect)),
                       type = "quantile",
                       type_id = quantile_vect,
                       .by = all_of(grep("type", peak_group,
                                         value = TRUE, invert = TRUE))) %>%
        dplyr::mutate(horizon = NA, target = peak_size_target)
    }
    df_all <- rbind(df_all, df_size)
  } else {
    message("The submission file from: ", team_model,
            " contains value for the target: ", peak_size_target)
    df_all <- rbind(df_all, df_size_target)
  }
  
  return(df_all)
}

#' Calculate missing peak time
#'
#' For rounds with samples output_type, the peak time hosp target is optional.
#' This function test for if peak time hosp is present in a specific team_model
#' submission and if not generate them
#'
#' @param df_team model projection of a team-model (as submitted)
#' @param df_all model projection of a team-model (with samples and quantiles)
#' to append the results, NULL for no data to append
#' @param team_model name of the team-model
#' @param peak_time_target target name as in the submission file
#' @param peak_group character vector of column names to 'group_by' the input
#' data frame to calculate the peak time hospitalization target for each group
#'
#' @details
#' Generates missing target from incident hospitalization trajectories by
#' calculating cumulative probability of peak across the complete time series.
#' For each trajectories, the horizon of the maximum incident hospitalization
#' is extracted.
#' For each horizon: calculate the average of how many trajectories have peaked
#' on this week (for example: 45 trajectories on 100 trajectories have peaked on
#' horizon 4) (probability of peak).
#' Calculate the cumulative probability of peak across trajectories:
#' following previous example cumulative on horizon 4 = (45/100 + cumulative
#' probability on horizon 3)
#'
#' @export
#' @importFrom dplyr filter group_by across all_of mutate ungroup bind_rows
#' @importFrom dplyr arrange select group_split summarise distinct
#' @importFrom MMWRweek MMWRweek
calculate_peak_time <- function(
    df_team, df_all, team_model, peak_time_target = "peak time hosp",
    peak_group = c("origin_date", "scenario_id", "location", "target",
                   "type_id")) {
  
  df_time_target <- dplyr::filter(df_team, target == peak_time_target)
  if (nrow(df_time_target) < 1) {
    df_inc_sample <- dplyr::filter(df_team, type == "sample",
                                   target == "inc hosp") %>%
      dplyr::mutate(type_id = as.numeric(type_id),
                    origin_date = as.Date(origin_date))
    if (nrow(df_inc_sample) < 1) {
      warning("Missing inc hosp sample to calculate peak time hosp.")
      df_all <- df_all
    } else {
      message("The submission file from: ", team_model,
              " is missing value for the target: ", peak_time_target,
              ". The missing value will be calculated.")
      
      df_time <- dplyr::group_by(df_inc_sample, across(all_of(peak_group))) %>%
        mutate(sel = ifelse(max(value) == value, horizon, NA)) %>%
        mutate(sel2 = ifelse(all(value == max(value)), 0, 1)) %>%
        filter(sel2 != 0) %>%
        filter(!is.na(sel)) %>%
        arrange(horizon) %>%
        mutate(sel3 = ifelse(duplicated(type_id), 0, 1)) %>%
        filter(sel3 == 1) %>%
        select(-contains("sel")) %>%
        ungroup()
      
      lst_time <-
        dplyr::group_split(df_time,
                           dplyr::across(all_of(grep("type", peak_group,
                                                     value = TRUE,
                                                     invert = TRUE))))
      n_sample <- df_inc_sample %>%
        dplyr::summarise(n = n(), .by = grep("type|value",
                                             colnames(df_inc_sample),
                                             value = TRUE, invert = TRUE)) %>%
        .$n %>%
        unique()
      peak_time <- lapply(lst_time, function(dft) {
        df_epitime <- NULL
        if (nrow(dft) != n_sample) {
          print(paste0("No. Peaks != ", n_sample))
          print(nrow(dft))
        }
        for (i in 1:max(df_inc_sample$horizon, na.rm = TRUE)) {
          peak_prob <- nrow(dplyr::filter(dft, horizon == i)) / nrow(dft)
          if (!is.null(df_epitime)) {
            peak_cum <- filter(df_epitime, horizon == i - 1) %>% .$value
            peak_cum <- peak_cum + peak_prob
          }   else {
            peak_cum <- peak_prob
          }
          if (peak_cum >= 1) peak_cum <- 1
          date <- MMWRweek::MMWRweek(unique(dft$origin_date) + (i * 7) - 1)
          if (nchar(date$MMWRweek) < 2) {
            date <- paste0("EW", date$MMWRyear, "0", date$MMWRweek)
          } else {
            date <- paste0("EW", date$MMWRyear, date$MMWRweek)
          }
          df_epi <- distinct(dft[, grep("type", peak_group, value = TRUE,
                                        invert = TRUE)]) %>%
            mutate(horizon = i,
                   type = "cdf",
                   type_id = date,
                   value = peak_cum,
                   target = peak_time_target)
          df_epitime <- rbind(df_epi, df_epitime)
        }
        df_epitime$horizon <- NA
        return(df_epitime)
      }) %>%
        dplyr::bind_rows()
      # df_all <- dplyr::mutate(df_all, origin_date = as.Date(origin_date))
      df_all <- rbind(df_all, peak_time)
    }
  } else {
    message("The submission file from: ", team_model,
            " contains value for the target: ", peak_time_target)
    df_all <- rbind(df_all, df_time_target)
  }
  return(df_all)
}
#' Calculate LOP Ensemble(s) for peak targets
#'
#' Calculate LOP Ensemble (trimmed or untrimmed, with weight or not) for both
#' peak size hosp and peak time hosp by using
#' `CombineDistributions::aggregate_cdfs()`. The input file should match the
#' expected input format of this function.
#'
#' @param df_peak data frame containing the peak target information to calculate
#'  ensemble on, the data frame should be in the expected format
#' @param weighting_scheme string to indicate how to weight in the aggregate (
#'  see function documentation)
#' @param n_trim integer denoting the number of models to trim (
#'  see function documentation). If NA, no trim applied
#' @param weight_df data frame containing the weight information by team-model,
#'  should have 2 columns: `model_id` and `weight`.
#' @param ens_group character vector of group of column on which to calculate
#'  the ensemble (for example, aggregate by scenario, target, horizon, location,
#'  and age group for the FLU projections)
#' @param id_var string containing the name of the column that identifies
#' unique cdfs (parameter from `CombineDistributions::aggregate_cdfs()`.
#' By default `model_id`.
#'
#' @details
#' The function returns a data frame in the SMH standard format
#'
#'
#' @importFrom dplyr filter arrange mutate select
#' @importFrom CombineDistributions aggregate_cdfs
#' @export
peak_ensemble <- function(df_peak, weighting_scheme, n_trim, weight_df,
                          ens_group, id_var = "model_id") {
  
  df_peak_time <- dplyr::filter(df_peak, grepl("time", target)) %>%
    dplyr::arrange(target_end_date)
  df_peak_size <- dplyr::filter(df_peak, grepl("size", target))
  if (is.na(n_trim)) {
    ret_val <- unique(df_peak_time$target_end_date)
    df_lop_peak_time <-
      CombineDistributions::aggregate_cdfs(df_peak_time, id_var = id_var,
                                           group_by = ens_group, method = "LOP",
                                           weighting_scheme = weighting_scheme,
                                           weights = weight_df,
                                           ret_quantiles = NA,
                                           reorder_quantiles = FALSE,
                                           ret_values = ret_val) %>%
      dplyr::mutate(target_end_date = value, value = quantile,
                    type_id = make_epiweek(target_end_date),
                    type = "cdf", horizon = NA) %>%
      dplyr::select(-quantile, -target_end_date)
    
    ret_quant <- unique(df_peak_size$quantile)
    df_lop_peak_size <-
      CombineDistributions::aggregate_cdfs(df_peak_size, id_var = id_var,
                                           group_by = ens_group, method = "LOP",
                                           weighting_scheme = weighting_scheme,
                                           weights = weight_df,
                                           ret_quantiles = ret_quant,
                                           ret_values = NA)  %>%
      dplyr::mutate(horizon = NA, type = "quantile",
                    type_id = quantile) %>%
      dplyr::select(-quantile)
    
  } else {
    ret_val <- unique(df_peak_time$target_end_date)
    df_lop_peak_time <-
      CombineDistributions::aggregate_cdfs(df_peak_time, id_var = id_var,
                                           group_by = ens_group, method = "LOP",
                                           weighting_scheme = weighting_scheme,
                                           n_trim = n_trim, ret_quantiles = NA,
                                           reorder_quantiles = FALSE,
                                           ret_values = ret_val) %>%
      dplyr::mutate(target_end_date = value, value = quantile,
                    type_id = make_epiweek(target_end_date),
                    type = "cdf", horizon = NA) %>%
      dplyr::select(-quantile, -target_end_date)
    
    ret_quant <- unique(df_peak_size$quantile)
    df_lop_peak_size <-
      CombineDistributions::aggregate_cdfs(df_peak_size, id_var = id_var,
                                           group_by = ens_group, method = "LOP",
                                           weighting_scheme = weighting_scheme,
                                           n_trim = n_trim,
                                           ret_quantiles = ret_quant,
                                           ret_values = NA) %>%
      dplyr::mutate(horizon = NA, type = "quantile",
                    type_id = quantile) %>%
      dplyr::select(-quantile)
  }
  df_lop_peak <- rbind(df_lop_peak_size, df_lop_peak_time)
  return(df_lop_peak)
}


#' Format Peak time submission for ensemble calculation
#'
#' Re-format the model output to the expected format for peak time ensemble
#' calculation:
#' - transform the column `value` as `quantile`
#' - calculated `target_end_date` column date extracted from `output_type_id`
#'  column
#' - transform the `value` column as numeric version of the date extracted
#'  from `output_type_id` column
#'
#' @param df date frame, model output containing "peak time" target information
#'
#' @details The output data frame of this function, will only contains
#' `"peak time"` target information. Any additional information is removed from
#' the inputted data frame
#'
#' @noRd
#' @importFrom MMWRweek MMWRweek2Date
#' @importFrom dplyr filter mutate arrange
prep_peak_time <- function(df) {
  df_time <- dplyr::filter(df, grepl("peak time", target)) %>%
    dplyr::mutate(quantile = value,
                  target_end_date =
                    MMWRweek2Date(as.numeric(substr(gsub("EW", "",
                                                         type_id), 1,
                                                    4)),
                                  as.numeric(substr(gsub("EW", "",
                                                         type_id), 5,
                                                    6)), 7),
                  value =
                    MMWRweek2Date(as.numeric(substr(gsub("EW", "",
                                                         type_id), 1,
                                                    4)),
                                  as.numeric(substr(gsub("EW", "",
                                                         type_id), 5,
                                                    6)), 7)) %>%
    dplyr::mutate(value = as.numeric(value)) %>%
    dplyr::arrange(target_end_date)
  return(df_time)
}


#' Format Peak size submission for ensemble calculation
#'
#' Re-format the model output to the expected format for peak size ensemble
#' calculation:
#' - filtered `output_type_id` to include the value inputted in the
#'  `list_quantile` parameter
#' - added the column `quantile` as an equivalent to the `output_type_id` column
#' - added `target_end_date` column as a `NA` Date format
#'
#' @param df date frame, model output containing "peak size" target information
#' @param list_quantiles numeric vector of probabilities, to select
#'  quantiles of interest in the submission files
#'
#' @details The output data frame of this function, will only contains
#' `"peak size"` target information. Any additional information is removed from
#' the inputted data frame
#'
#' @noRd
#' @importFrom dplyr filter mutate arrange
prep_peak_size <- function(df, list_quantiles) {
  df_size <- dplyr::filter(df, grepl("peak size", target),
                           !is.na(type_id)) %>%
    dplyr::arrange(type_id) %>%
    dplyr::filter(type_id %in% list_quantiles) %>%
    dplyr::mutate(quantile = type_id,
                  target_end_date = as.Date(NA))
  return(df_size)
}

#' 
#' #' Adapt peak target information for ensemble calculation
#' #'
#' #' Re-frame the input data frame and generate missing variable in the expected
#' #' format to be able to calculate LOP ensembles on peak targets.
#' #'
#' #' @param path_model path to the folder containing the submission files to use
#' #'  to calculate the peak ensemble
#' #' @param list_quantiles numeric vector of probabilities, to select
#' #'  quantiles of interest in the submission files (peak size only)
#' #' @param peak_target character vector, peak name target to include in the
#' #' ensemble calculation
#' #' @param partition character vector indicating if the submission file is
#' #'  partitioned and if so, which field (or column) names correspond to the path
#' #'  segments. By default, `NULL` (no partition).
#' #' @param schema arrow schema object, use to load the partition files.
#' #' By default, `NULL`
#' #' @param orig_date origin date used to identify a specific round, should be
#' #' in the "YYYY-MM-DD" format. By default, `NULL`.
#' #' @param round_path path the a csv with round information. `NULL` by default
#' #' @param round_number numeric to select data for a specific round number,
#' #'  round number should be written  like "1", "5", ... and the parameter
#' #'  `round_path` should be filled too. If `NULL` (default), all round selected
#' #' @param excl_model character string, regex of team-model name to exclude from
#' #' the file, if `NULL` (default) no exclusion
#' #'
#' #' #' @details
#' #' ## Load the data
#' #'
#' #' The function contains two possible ways to load the data.
#' #'
#' #' 1. By using the `round_path, round_number` parameters to identify the files
#' #' to load.
#' #'
#' #' 2. (Recommended for partitioned files) By using the `partition`, `schema`,
#' #' and `orig_date` to identify which data to load. Used the
#' #' `arrow::open_dataset()` function
#' #'
#' #' As the parameters are set to `NULL` by default, a method would need to be
#' #' chosen and set accordingly.
#' #'
#' #' @importFrom dplyr mutate filter arrange bind_rows
#' #' @importFrom MMWRweek MMWRweek2Date
#' #' @noRd
#' peak_prep <- function(path_model, list_quantiles, peak_target, partition = NULL,
#'                       schema = NULL, orig_date = NULL, round_number = NULL,
#'                       round_path = NULL, excl_model = NULL) {
#'   if (is.null(partition)) {
#'     pttrn <- ".{4}-.{2}-.{2}-|.csv|.zip|.gz|.pqt|.parquet"
#'     # List of files to read per round
#'     name_peak <- path_files(path_model, ensemble = FALSE,
#'                             round_path = round_path, round = round_number)
#'     # Files exclusion if necessary
#'     if (!is.null(excl_model)) {
#'       name_peak <- grep(excl_model, name_peak, value = TRUE, invert = TRUE)
#'     }
#'     # Re-Frame data and calculate missing variable
#'     df_peak <- lapply(name_peak, function(x) {
#'       df <- read_files(x) %>%
#'         dplyr::filter(target %in% peak_target) %>%
#'         dplyr::mutate(origin_date = as.Date(origin_date))
#'       # Peak time
#'       if (nrow(df) > 0) {
#'         df_time <- prep_peak_time(df)
#'         df_size <- prep_peak_size(df, list_quantiles)
#'         df_tot <- rbind(df_time, df_size) %>%
#'           dplyr::mutate(model_id = gsub(pttrn, "", basename(x)))
#'       } else {
#'         df_tot <- NULL
#'       }
#'       return(df_tot)
#'     }) %>%
#'       dplyr::bind_rows()
#'   } else {
#'     df_model <- hubData::connect_model_output(path_model,
#'                                               partition_names = partition,
#'                                               file_format = "parquet",
#'                                               schema = schema) %>%
#'       dplyr::filter(!grepl("(e|E)nsemble", model_id), target %in% peak_target,
#'                     origin_date == as.Date(orig_date)) %>%
#'       hubData::collect_hub() %>%
#'       dplyr::rename(model_id = model_id)
#'     if (!is.null(excl_model))
#'       df_model <- dplyr::filter(df_model, !grepl(excl_model, model_id))
#'     df_time <- prep_peak_time(df_model)
#'     df_size <- prep_peak_size(df_model, list_quantiles)
#'     df_peak <- rbind(df_time, df_size)
#'   }
#'   df_peak <- dplyr::mutate(df_peak,
#'                            location = as.character(ifelse(nchar(location) == 1,
#'                                                           paste0("0", location),
#'                                                           location)))
#'   return(df_peak)
#' }

#' Calculate the epiweek in "EWYYYWW" format
#'
#' For a specific date or a vector of date, calculate the epiweek in a
#' "EWYYYYWW" format
#'
#' @param date date or vector of date to transform
#'
#' @export
#' @importFrom MMWRweek MMWRweek
make_epiweek <- function(date) {
  epi_date <- MMWRweek::MMWRweek(date)
  epi_year <- epi_date[["MMWRyear"]]
  epi_week <- as.character(epi_date[["MMWRweek"]])
  epi_week <- ifelse(nchar(epi_week) < 2, paste0(0, epi_week), epi_week)
  epi_info <- paste0("EW", epi_year, epi_week)
  return(epi_info)
}