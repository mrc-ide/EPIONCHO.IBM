process_multiple_runs <- function(files='', outputs_per_year = 4, ov16_indiv = FALSE, morbidity_runs = FALSE, verbose=TRUE, ov16_indiv_locations=c(1), ov16_indiv_location_names=c("baseline"), min_iter=0, max_iter=100000) {
  allOutputs <- data.frame()
  fileToUse <- files

  i <- 1
  total_files <- length(list.files(fileToUse))
  start_time <- Sys.time() 
  mf_prev_df <- NA
  atp_df <- NA
  mf_intensity_df <- NA
  morbidity_df <- NA
  worm_df <- NA
  oae_incidence_df <- NA
  ov16_trends_df <- NA
  ov16_adj_trends_df <- NA
  ov16_indiv_df <- NA
  for (file in list.files(fileToUse)) {
    if(verbose) {
      if (total_files > 10 & (i %% floor(total_files/10)) == 0) {
        print(paste("Time Elapsed:", Sys.time()-start_time, ":", i / (total_files)))
        gc()
      }
    }
    iter = -1
    tryCatch(
      {
        # Assume that the iter number is always the last value before '.rds'
        # I.e output_...._17.rds, 17 is the iter number
        all_vals <- strsplit(file, "_")[[1]]
        iter <- as.numeric(strsplit(all_vals[length(all_vals)], "\\.")[[1]][1])
        print(paste0("i: ", i, " Iter: ", iter))
        if (iter > max_iter | iter < min_iter) {
          next
        }
        tmpRDSData <- readRDS(paste(fileToUse, file, sep=""))
      },
      error = function(e) {
        message(paste("Error occurred while reading:", fileToUse))
      }
    )

    
    print(file)
    timesteps_in_year = 1/tmpRDSData$timestep_used
    # selector <- which(tmpRDSData$years %% 1 == 0)
    total_timesteps <- 1:(length(tmpRDSData$years))
    year_position <- total_timesteps %% timesteps_in_year
    year_position[year_position == 0] <- timesteps_in_year
    output_locs_year <- floor(seq(1, timesteps_in_year, length.out=(outputs_per_year+1))[2:(outputs_per_year+1)])
    selector <- total_timesteps[year_position %in% output_locs_year]

    num_vals <- length(selector)
    
    # Used for Ov16 individual matrix processing
    ov16_indiv_matrix <- tmpRDSData$ov16_indiv_matrix
    age <- ov16_indiv_matrix[, ov16_indiv_locations[1]]
    num_individuals <- length(age)
    ov16_indiv_df_cols <- ncol(ov16_indiv_matrix) + 3
    if(i == 1) {
      mf_prev_df <- matrix(ncol=(ncol(tmpRDSData$all_mf_prevalence_age_grouped)*1)+4, nrow=total_files*num_vals)
      colnames(mf_prev_df) <- c(
        colnames(tmpRDSData$all_mf_prevalence_age_grouped), 
        # paste0("30-", colnames(tmpRDSData$all_mf_prevalence_age_grouped)), - TODO: need to update model
        "Ke", "ABR", "year", 'run_num'
      )

      atp_df <- matrix(ncol=11, nrow=total_files*num_vals)
      colnames(atp_df) <- c("l3_per_blackfly", "l3_prevalence_per_blackfly", "recorded_abr", "atp", "Ke", "ABR", "h", "m", "beta", "year", 'run_num')

      mf_intensity_df <- matrix(ncol=(ncol(tmpRDSData$all_mf_intensity_age_grouped)*1)+4, nrow=total_files*num_vals)
      colnames(mf_intensity_df) <- c(
        colnames(tmpRDSData$all_mf_intensity_age_grouped), 
        # paste0("30-", colnames(tmpRDSData$all_mf_intensity_age_grouped)), - TODO: need to update model
        "Ke", "ABR", "year", 'run_num'
      )

      worm_df <- matrix(ncol=(ncol(tmpRDSData$worm_burden_outputs))+4, nrow=total_files*num_vals)
      colnames(worm_df) <- c(
        colnames(tmpRDSData$worm_burden_outputs),
        "Ke", "ABR", "year", 'run_num'
      )
      
      if (morbidity_runs) {
        morbidity_df <- matrix(ncol=ncol(tmpRDSData$all_morbidity_prevalence_outputs)+4, nrow=total_files*num_vals)
      colnames(morbidity_df) <- c(colnames(tmpRDSData$all_morbidity_prevalence_outputs), "Ke", "ABR", "year", 'run_num')

        oae_incidence_df <- matrix(ncol=ncol(tmpRDSData$oae_incidence_outputs)+4, nrow=total_files*num_vals)
      colnames(oae_incidence_df) <- c(colnames(tmpRDSData$oae_incidence_outputs), "Ke", "ABR", "year", 'run_num')
      }
      
      ov16_trends_df <- matrix(ncol=ncol(tmpRDSData$ov16_timetrend_outputs)+4, nrow=total_files*num_vals)
      colnames(ov16_trends_df) <- c(colnames(tmpRDSData$ov16_timetrend_outputs), "Ke", "ABR", "year", 'run_num')
      
      ov16_adj_trends_df <- matrix(ncol=ncol(tmpRDSData$ov16_timetrend_outputs_adj)+4, nrow=total_files*num_vals)
      colnames(ov16_adj_trends_df) <- c(colnames(tmpRDSData$ov16_timetrend_outputs_adj), "Ke", "ABR", "year", 'run_num')
      
      if (ov16_indiv) {
        ov16_indiv_df <- data.frame(matrix(ncol=ov16_indiv_df_cols, nrow=num_individuals*total_files))
        ov16_indiv_df_final_col_names <- c()
        ov16_indiv_df_names <- c("age", "sex", "mf_status", "mf_intens_raw", "ov16_status", "serorevert_fast")
        ov16_indiv_index <- 1
        for (loc in ov16_indiv_locations) {
          ov16_indiv_df_final_col_names <- c(
            ov16_indiv_df_final_col_names,
            paste(ov16_indiv_location_names[loc], ov16_indiv_df_names, sep="_")
          )
        }
        ov16_indiv_df_final_col_names <- c(ov16_indiv_df_final_col_names, "ABR", "Ke", "run_num")
        colnames(ov16_indiv_df) <- ov16_indiv_df_final_col_names
      }
    }
    
    kE <- -1
    if (any("Ke" %in% names(tmpRDSData))) {
      kE <- tmpRDSData$Ke
    }

    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-10)] <- tmpRDSData$L3[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-9)] <- tmpRDSData$blackfly_l3_prevalence[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-8)] <- tmpRDSData$ABR_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-7)] <- tmpRDSData$L3[selector] * tmpRDSData$ABR_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-6)] <- rep(kE, length(selector))
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-5)] <- tmpRDSData$ABR
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-4)] <- NaN# tmpRDSData$h_recorded[selector] - TODO: need to update model
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-3)] <- NaN# tmpRDSData$m_recorded[selector] - TODO: need to update model
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-2)] <- NaN# tmpRDSData$beta_recorded[selector] - TODO: need to update model
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-1)] <- tmpRDSData$years[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df))] <- rep(iter, length(selector))


    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(tmpRDSData$all_mf_prevalence_age_grouped))] <- tmpRDSData$all_mf_prevalence_age_grouped[selector,]
    # mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(tmpRDSData$all_mf_prevalence_age_grouped)+1):(ncol(mf_prev_df)-4)] <- tmpRDSData$all_mf_prevalence_30_age_grouped[selector,]
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)-3] <- rep(kE, length(selector))
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)-2] <- tmpRDSData$ABR
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)-1] <- tmpRDSData$years[selector]
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)] <- rep(iter, length(selector))
    
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(tmpRDSData$all_mf_intensity_age_grouped))] <- tmpRDSData$all_mf_intensity_age_grouped[selector,]
    # mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(tmpRDSData$all_mf_intensity_age_grouped)+1):(ncol(mf_intensity_df)-4)] <- tmpRDSData$all_mf_intensity_30_age_grouped[selector,]
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)-3] <- rep(kE, length(selector))
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)-2] <- tmpRDSData$ABR
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)-1] <- tmpRDSData$years[selector]
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)] <- rep(iter, length(selector))

    worm_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(tmpRDSData$worm_burden_outputs))] <- tmpRDSData$worm_burden_outputs[selector,]
    worm_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(worm_df)-3] <- rep(kE, length(selector))
    worm_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(worm_df)-2] <- tmpRDSData$ABR
    worm_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(worm_df)-1] <- tmpRDSData$years[selector]
    worm_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(worm_df)] <- rep(iter, length(selector))
    
    if (morbidity_runs) {
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(morbidity_df)-4)] <-tmpRDSData$all_morbidity_prevalence_outputs[selector,]
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)-3] <- rep(kE, length(selector))
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)-2] <- tmpRDSData$ABR
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)-1] <- tmpRDSData$years[selector]
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)] <- rep(iter, length(selector))
  
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(oae_incidence_df)-4)] <- tmpRDSData$oae_incidence_outputs[selector,]
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)-3] <- rep(kE, length(selector))
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)-2] <- tmpRDSData$ABR
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)-1] <- tmpRDSData$years[selector]
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)] <- rep(iter, length(selector))
    }
    

    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(ov16_trends_df)-4)] <- tmpRDSData$ov16_timetrend_outputs[selector,]
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)-3] <- rep(kE, length(selector))
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)-2] <- tmpRDSData$ABR
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)-1] <- tmpRDSData$years[selector]
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)] <- rep(iter, length(selector))
    
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(ov16_adj_trends_df)-4)] <- tmpRDSData$ov16_timetrend_outputs_adj[selector,]
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)-3] <- rep(kE, length(selector))
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)-2] <- tmpRDSData$ABR
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)-1] <- tmpRDSData$years[selector]
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)] <- rep(iter, length(selector))
    

    if (ov16_indiv) {
      start_index <- 1+num_individuals*(i-1)
      end_index <- num_individuals*i
      ov16_indiv_df[start_index:end_index,1:(ov16_indiv_df_cols-3)] <- ov16_indiv_matrix
      ov16_indiv_df[start_index:end_index,ov16_indiv_df_cols-2] <- tmpRDSData$ABR
      ov16_indiv_df[start_index:end_index,ov16_indiv_df_cols-1] <- kE
      ov16_indiv_df[start_index:end_index,ov16_indiv_df_cols] <- iter
    }

    i <- i + 1

  }

  if (ov16_indiv) {
    ov16_indiv_df <- as.data.frame(ov16_indiv_df)
    ov16_indiv_df["age_groups"] <- ifelse(
      ceiling(ov16_indiv_df$age*5)/5 == 0, 
      0,
      ifelse(ov16_indiv_df$age <= 75, ceiling(ov16_indiv_df$age/5)*5 - 2.5, 77.5)
    )
  }
  all_return <- list(mf_prev_df, mf_intensity_df, worm_df, morbidity_df, oae_incidence_df, ov16_trends_df, ov16_adj_trends_df, ov16_indiv_df, atp_df)
  names(all_return) <- c("mf_prev_df", "mf_intensity_df", "worm_df", "morbidity_df", "oae_incidence_df", "ov16_trends_df", "ov16_adj_trends_df", "ov16_indiv_df", "atp_df")
  return(all_return)
}

outputs_folder = Sys.getenv("OUTPUT_FOLDER")

processed_file_path = Sys.getenv("PROCESS_DATA_PATH")
print(outputs_folder)
print(processed_file_path)

saveRDS(process_multiple_runs(file=paste0(outputs_folder, "/"), morbidity_runs=TRUE), processed_file_path)