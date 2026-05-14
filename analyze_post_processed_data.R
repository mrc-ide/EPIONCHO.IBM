library(dplyr)
library(ggplot2)

# plot all age values. If you want to change the age group of what is plotted, just edit the filters in the function
create_time_plot <- function(r_df, desc) {
  return(
    ggplot() +
      geom_line(
        data=r_df %>% filter((age_start == 0  | age_start == 5) & age_end == 80),
        aes(
          x=year,
          y=median_val
        )
      ) +
      geom_ribbon(
        data=r_df %>% filter((age_start == 0  | age_start == 5) & age_end == 80),
        aes(
          x=year,
          ymin=lower,
          ymax=upper
        ),
        alpha=0.3
      )
      facet_grid(Ke~ABR) +
      theme_bw() +
      labs(
        x = "Year",
        y = desc
      )
  )
}

# plot age grouped values. If you want to change the age groups plotted, just edit the filters in the function
create_compare_plot_age_groups <- function(r_df, desc, years=c(50, 95)) {
  x_breaks = r_df %>% filter(!((age_start == 0 | age_start == 5) & age_end == 80)) %>% distinct(age_start, age_end) %>% rowMeans() %>% floor()
  return(
    ggplot() +
      geom_line(
        data=r_df %>% filter(year %in% years & !((age_start == 0  | age_start == 5) & age_end == 80)),
        aes(
          x=(age_start + age_end) / 2,
          y=mean_val
        )
      ) +
      geom_ribbon(
        data=r_df %>% filter(year %in% years & !((age_start == 0  | age_start == 5) & age_end == 80)),
        aes(
          x=(age_start + age_end) / 2,
          ymin=lower,
          ymax=upper
        ),
        alpha=0.3
      ) +
      facet_grid(Ke~ABR+year, nrow=3) +
      theme_bw() +
      labs(
        y = desc
      ) +
      scale_x_continuous("Age Groups", breaks=x_breaks, limits=c(0, 81))
  )
}

create_baseline_measure_across_ke_abr <- function(r_df, desc, year_of_baseline = 90) {
  r_df <- r_df %>% 
    mutate(
      Ke = factor(
        Ke,
        levels = as.character(sort(unique(Ke)))
      )
    )
  return(
    ggplot() +
      geom_line(
        data=r_df %>% filter((age_start == 0  | age_start == 5) & age_end == 80 & year == year_of_baseline),
        aes(
          x=ABR,
          y=median_val,
          color=Ke
        )
      ) +
      geom_ribbon(
        data=r_df %>% filter((age_start == 0  | age_start == 5) & age_end == 80 & year == year_of_baseline),
        aes(
          x=ABR,
          ymin=lower,
          ymax=upper,
          fill=Ke
        ),
        alpha=0.3
      ) +
      theme_bw() +
      labs(
        x = "ABR",
        y = desc,
        colour = "Ke",
        fill = "Ke"
      )
  )
}

mutate_r_data <- function(df, measure_name, measure_delimiter = "_", scale_val = 100, rename_measure = "") {
  start_time = Sys.time();
  new_label <- if (measure_delimiter != "") {
    gsub(measure_delimiter, "-", measure_name, fixed = TRUE)
  } else {
    measure_name
  }
  if (rename_measure == "") {
    rename_measure = new_label
  }
  
  # Convert to data table for quicker operations
  dt <- data.table::as.data.table(df)
  cols <- names(dt)
  if (measure_delimiter != "") {
    pattern <- paste0("^", measure_name, "(?=[0-9])")
    matches <- grepl(pattern, cols, perl = TRUE)
    cols[matches] <- sub(pattern, paste0(new_label, "_"), cols[matches], perl = TRUE)
    cols[cols == measure_name] <- new_label
    data.table::setnames(dt, names(dt), cols)
  }
  
  if (new_label %in% names(dt)) {
    data.table::setnames(dt, new_label, paste0(new_label, "_0_80"))
  }
  
  target_cols <- grep(paste0("^", new_label, "_"), names(dt), value = TRUE)
  
  long <- data.table::melt(
    dt,
    measure.vars    = target_cols,
    variable.name   = "age_groups",
    value.name      = "prevalence",
    variable.factor = FALSE
  )
  
  long[, c("age_start", "age_end") := {
    parts <- regmatches(age_groups, regexec("_(\\d+)_(\\d+)$", age_groups))
    age_s <- as.integer(vapply(parts, function(x) if (length(x) >= 3) x[[2]] else NA_character_, character(1)))
    age_e <- as.integer(vapply(parts, function(x) if (length(x) >= 3) x[[3]] else NA_character_, character(1)))
    list(age_s, age_e)
  }]
  
  result <- long[, .(
    lower = as.numeric(quantile(prevalence, 0.025)) * scale_val,
    upper = as.numeric(quantile(prevalence, 0.975)) * scale_val,
    median_val = median(prevalence) * scale_val,
    mean_val = mean(prevalence) * scale_val,
    measure = rename_measure
  ), by = .(ABR, Ke, year, age_start, age_end)]
  cat(paste0(round(difftime(Sys.time(), start_time, units = 'secs'), 2), " secs"))
  return(as.data.frame(result) %>% select(ABR, Ke, year, age_start, age_end, measure, everything()))
}

mutate_all_r_data <- function(r_data, all_ke_values, all_abr_values) {
  print("Starting Processing for all values: ")
  tmp_r_processed_data_list <- list()
  ind <- 1
  for (abr_val in abr_vals) {
    cat(paste0("Start | abr ", abr_val))
    for (key in names(all_r_measures)) {
      cat(paste0(" | ", key, "..."))
      delim = all_r_measures[[key]]
      if (key == "prev")
        tmp_r_processed_data_list[[ind]] <- mutate_r_data(r_data$mf_prev_df %>% filter_vals(key, all_ke_values, c(abr_val)), key, delim, scale_val = 100, rename_measure=r_measure_to_cpp_measure_map[[key]])
      if (key == "mean.mf.per.snip")
        tmp_r_processed_data_list[[ind]] <- mutate_r_data(r_data$mf_intensity_df %>% filter_vals(key, all_ke_values, c(abr_val)), key, delim, scale_val = 1, rename_measure=r_measure_to_cpp_measure_map[[key]])
      if (key == "ov16_seroprevalence_combined_seroreversion")
        tmp_r_processed_data_list[[ind]] <- mutate_r_data(r_data$ov16_trends_df %>% filter_vals(key, all_ke_values, c(abr_val)), key, delim, scale_val = 100, rename_measure=r_measure_to_cpp_measure_map[[key]])
      if (key == "ov16_seroprevalence_combined_seroreversion_adj")
        tmp_r_processed_data_list[[ind]] <- mutate_r_data(r_data$ov16_adj_trends_df %>% filter_vals(key, all_ke_values, c(abr_val)), key, delim, scale_val = 100, rename_measure=r_measure_to_cpp_measure_map[[key]])
      if (grepl("_prev", key))
        tmp_r_processed_data_list[[ind]] <- mutate_r_data(r_data$morbidity_df %>% filter_vals(key, all_ke_values, c(abr_val)), key, delim, scale_val = 100, rename_measure=r_measure_to_cpp_measure_map[[key]])
      if (grepl("worm_burden", key))
        tmp_r_processed_data_list[[ind]] <- mutate_r_data(r_data$worm_df %>% filter_vals(key, all_ke_values, c(abr_val)), key, delim, scale_val = 1, rename_measure=r_measure_to_cpp_measure_map[[key]])
      ind = ind + 1
    }
    cat(" done\n")
  }
  return(bind_rows(tmp_r_processed_data_list))
}

# Comment out measures you don't need to reduce the run time of the summary code
all_r_measures = list(
  "prev"="",
  "mean.mf.per.snip"="\\.",
  "ov16_seroprevalence_combined_seroreversion"="_",
  "ov16_seroprevalence_combined_seroreversion_adj"="_",
  "SevereItch_prev"="_",
  "RSD_prev"="_",
  "Atrophy_prev"="_",
  "HG_prev"="_",
  "Depig_prev"="_",
  "Blindness_prev"="_",
  "OAE_prev"="_",
  "male_worm_burden"="_",
  "fertile_female_worm_burden"="_",
  "infertile_female_worm_burden"="_"
)

filter_vals <- function(df, key="", k_E = c(0.2), abr = c(1000)) {
  return (
    df %>% 
      as.data.frame() %>% 
      filter(Ke %in% k_E & ABR %in% abr) %>%
      select(Ke, ABR, year, run_num, starts_with(key))
  )
}

r_data <- readRDS("path/to/processed/data.RDS")
all_kEs <- unique(as.data.frame(r_data$atp_df)$Ke)
abr_vals <- sort(unique(as.data.frame(r_data$atp_df)$ABR))

# Rather than summarising data all abrs and kes, you can change the function below
# to only summarise smaller subsets, which should take a lot less time than everything at once.
all_r_processed_data <- mutate_all_r_data(r_data, all_kEs, abr_vals)

# Save the summarised data. Since it takes a long time, you don't want to sit and wait for the processing again
# instead you can just load the file
saveRDS(all_r_processed_data, "path/to/summarised/data.RDS")


for (abr_val in abr_vals) {
  for (m in unique(all_r_processed_data$measure)) {
    print(create_compare_plot(all_r_processed_data %>% filter(ABR == abr_val & measure == m), m))
    print(create_compare_plot_age_groups(all_r_processed_data %>% filter(ABR == abr_val & measure == m), m))
  }
}

for (m in unique(all_r_processed_data$measure)) {
  print(create_baseline_measure_across_ke_abr(all_r_processed_data %>% filter(measure == m), m))
}