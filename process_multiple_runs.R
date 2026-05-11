process_multiple_runs <- function (files='', ov16_indiv = FALSE, morbidity_runs = FALSE, verbose=TRUE, ov16_indiv_locations=c(1, 2, 3, 4, 5)) {
  allOutputs <- data.frame()
  fileToUse <- files

  i <- 1
  total_files <- length(list.files(fileToUse))
  start_time <- Sys.time() 
  mf_prev_df <- NA
  mf_intensity_df <- NA
  morbidity_df <- NA
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
    tryCatch(
      {
        tmpRDSData <- readRDS(paste(fileToUse, file,sep=""))
      },
      error = function(e) {
        message(paste("Error occurred while reading:", fileToUse))
      }
    )
    
    print(file)
    selector <- which(tmpRDSData$years %% 1 == 0)
    num_vals <- length(selector)
    
    # Used for Ov16 individual matrix processing
    ov16_indiv_matrix <- tmpRDSData$ov16_indiv_matrix
    age <- ov16_indiv_matrix[, ov16_indiv_locations[1]]
    num_individuals <- length(age)
    ov16_indiv_df_cols <- ncol(ov16_indiv_matrix) + 3
    if(i == 1) {
      mf_prev_df <- matrix(ncol=ncol(tmpRDSData$all_mf_prevalence_age_grouped)+4, nrow=total_files*num_vals)
      colnames(mf_prev_df) <- c(colnames(tmpRDSData$all_mf_prevalence_age_grouped), "Ke", "ABR", "year", 'run_num')

      atp_df <- matrix(ncol=10, nrow=total_files*num_vals)
      colnames(atp_df) <- c("l3_per_blackfly", "recorded_abr", "atp", "Ke", "ABR", "h", "m", "beta", "year", 'run_num')

      mf_intensity_df <- matrix(ncol=ncol(tmpRDSData$all_mf_intensity_age_grouped)+4, nrow=total_files*num_vals)
      colnames(mf_intensity_df) <- c(colnames(tmpRDSData$all_mf_intensity_age_grouped), "Ke", "ABR", "year", 'run_num')
      
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
        ov16_indiv_df_col_names <- c()
        ov16_loc_names <- c("pre_mda", "post_vc_1", "post_vc_2", "start_biannual_mda", "survey_time")
        ov16_indiv_df_cols <- c("age", "sex", "mf_status", "mf_intens_raw", "ov16_status_no_seroreversion", "ov16_status_instant_serorevert", "ov16_status_finite_seroreversion")
        ov16_indiv_index <- 1
        for (loc in ov16_indiv_locations) {
          ov16_indiv_df_col_names <- c(
            ov16_indiv_df_col_names,
            paste(ov16_loc_names[loc], ov16_indiv_df_cols, sep="_")
          )
        }
        ov16_indiv_df_col_names <- c(ov16_indiv_df_col_names, "ABR", "Ke", "run_num")
        colnames(ov16_indiv_df) <- ov16_indiv_df_col_names
      }
    }
    
    kE <- -1
    if (any("Ke" %in% names(tmpRDSData))) {
      kE <- tmpRDSData$Ke
    }

    iter <- as.numeric(strsplit(strsplit(file, "_")[[1]][5], "\\.")[[1]][1])
    print(paste0("i: ", i, " Iter: ", iter))

    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-9)] <- tmpRDSData$L3[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-8)] <- tmpRDSData$ABR_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-7)] <- tmpRDSData$L3[selector] * tmpRDSData$ABR_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-6)] <- rep(kE, length(selector))
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-5)] <- tmpRDSData$ABR
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-4)] <- tmpRDSData$h_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-3)] <- tmpRDSData$m_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-2)] <- tmpRDSData$beta_recorded[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df)-1)] <- tmpRDSData$years[selector]
    atp_df[(1+(num_vals*(i-1))):(i*num_vals),(ncol(atp_df))] <- rep(iter, length(selector))


    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(mf_prev_df)-4)] <- tmpRDSData$all_mf_prevalence_age_grouped[selector,]
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)-3] <- rep(kE, length(selector))
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)-2] <- tmpRDSData$ABR_recorded[selector]
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)-1] <- tmpRDSData$years[selector]
    mf_prev_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_prev_df)] <- rep(iter, length(selector))
    
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(mf_intensity_df)-4)] <- tmpRDSData$all_mf_intensity_age_grouped[selector,]
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)-3] <- rep(kE, length(selector))
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)-2] <- tmpRDSData$ABR_recorded[selector]
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)-1] <- tmpRDSData$years[selector]
    mf_intensity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(mf_intensity_df)] <- rep(iter, length(selector))
    
    if (morbidity_runs) {
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(morbidity_df)-4)] <-tmpRDSData$all_morbidity_prevalence_outputs[selector,]
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)-3] <- rep(kE, length(selector))
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)-2] <- tmpRDSData$ABR_recorded[selector]
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)-1] <- tmpRDSData$years[selector]
      morbidity_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(morbidity_df)] <- rep(iter, length(selector))
  
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(oae_incidence_df)-4)] <- tmpRDSData$oae_incidence_outputs[selector,]
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)-3] <- rep(kE, length(selector))
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)-2] <- tmpRDSData$ABR_recorded[selector]
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)-1] <- tmpRDSData$years[selector]
      oae_incidence_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(oae_incidence_df)] <- rep(iter, length(selector))
    }
    

    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(ov16_trends_df)-4)] <- tmpRDSData$ov16_timetrend_outputs[selector,]
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)-3] <- rep(kE, length(selector))
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)-2] <- tmpRDSData$ABR_recorded[selector]
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)-1] <- tmpRDSData$years[selector]
    ov16_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_trends_df)] <- rep(iter, length(selector))
    
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),1:(ncol(ov16_adj_trends_df)-4)] <- tmpRDSData$ov16_timetrend_outputs_adj[selector,]
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)-3] <- rep(kE, length(selector))
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)-2] <- tmpRDSData$ABR_recorded[selector]
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)-1] <- tmpRDSData$years[selector]
    ov16_adj_trends_df[(1+(num_vals*(i-1))):(i*num_vals),ncol(ov16_adj_trends_df)] <- rep(iter, length(selector))
    

    if (ov16_indiv) {
      start_index <- 1+num_individuals*(i-1)
      end_index <- num_individuals*i
      # list_of_ov16_indiv_outputs <- list()
      # loc_index <- 0
      # for (loc in ov16_indiv_locations) {
      #   list_of_ov16_indiv_outputs[[(7*loc_index)+1]] <- ov16_indiv_matrix
      #   list_of_ov16_indiv_outputs[[(7*loc_index)+2]] <- ifelse(ov16_indiv_matrix[,(loc*7)-5] == 1, "Male", "Female")
      #   list_of_ov16_indiv_outputs[[(7*loc_index)+3]] <- ov16_indiv_matrix[,(loc*7)-1]
      #   list_of_ov16_indiv_outputs[[(7*loc_index)+4]] <- ov16_indiv_matrix[,(loc*7)]
      #   list_of_ov16_indiv_outputs[[(7*loc_index)+5]] <- ov16_indiv_matrix[,(loc*7)-2]
      #   loc_index + 1
      # }
      ov16_indiv_df[start_index:end_index,1:(ov16_indiv_df_cols-3)] <- ov16_indiv_matrix#list_of_ov16_indiv_outputs
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
  all_return <- list(mf_prev_df, mf_intensity_df, morbidity_df, oae_incidence_df, ov16_trends_df, ov16_adj_trends_df, ov16_indiv_df, atp_df)
  names(all_return) <- c("mf_prev_df", "mf_intensity_df", "morbidity_df", "oae_incidence_df", "ov16_trends_df", "ov16_adj_trends_df", "ov16_indiv_df", "atp_df")
  return(all_return)
}

saveRDS(process_multiple_runs(file="path/to/outputs/"), "path/to/save/file.RDS")