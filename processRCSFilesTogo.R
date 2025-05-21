library(dplyr)
processRCSFiles <- function (files='/', which_timestep=2, verbose=TRUE, onlyCalcMFP=FALSE, useSerorevert="none", onlyCalcOv16Trends=FALSE) {
  allOutputs <- data.frame()
  fileToUse <- files

  i <- 1
  mda_vals <- c()
  abr_vals <- c()
  mf_prev_vals <- c()
  total_files <- length(list.files(fileToUse))
  start_time <- Sys.time() 
  mf_prev_df <- data.frame()
  ov16_trend_df <- matrix(ncol=5, nrow=0)
  ov16_serorevert_trend_df <- matrix(ncol=5, nrow=0)
  for (file in list.files(fileToUse)) {
    if(verbose) {
      if ((i %% floor(total_files/10)) == 0) {
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
    
    file_parts <- unlist(strsplit(file, "_"))
    run_num <- as.numeric(sub("\\.rds$", "", file_parts[length(file_parts)]))

    mf_prev_vals <- c(mf_prev_vals, tmpRDSData$mf_prev[(119/(1/366))])
    mda_vals <- c(mda_vals, tmpRDSData$MDA)
    abr_vals <- c(abr_vals, tmpRDSData$ABR)
    vctr.ctrl.val = tmpRDSData$vctr.ctrl.eff
    if(onlyCalcOv16Trends) {
      sero_prev_vals <- tmpRDSData$ov16_seroprevalence[c(1,seq(from=183, to=length(tmpRDSData$ov16_seroprevalence), by=183))]
      serorev_prev_vals <- tmpRDSData$ov16_seroprevalence_sero[c(1,seq(from=183, to=length(tmpRDSData$ov16_seroprevalence_sero), by=183))]
      total_sero_vals <- length(sero_prev_vals)
      if(i == 1) {
        ov16_trend_df <- matrix(ncol=5, nrow=total_files*total_sero_vals)
        colnames(ov16_trend_df) <- c("ABR", "vctr.ctrl.eff", "Ke", "run_num", "ov16_vals")
        ov16_serorevert_trend_df <- matrix(ncol=5, nrow=total_files*total_sero_vals)
        colnames(ov16_serorevert_trend_df) <- c("ABR", "vctr.ctrl.eff", "Ke", "run_num", "ov16_vals")
      }
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),1] <- rep(tmpRDSData$ABR, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),2] <- rep(vctr.ctrl.val, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),3] <- rep(tmpRDSData$Ke, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),4] <- rep(run_num, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),5] <- sero_prev_vals

      ov16_serorevert_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),1] <- rep(tmpRDSData$ABR, total_sero_vals)
      ov16_serorevert_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),2] <- rep(vctr.ctrl.val, total_sero_vals)
      ov16_serorevert_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),3] <- rep(tmpRDSData$Ke, total_sero_vals)
      ov16_serorevert_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),4] <- rep(run_num, total_sero_vals)
      ov16_serorevert_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),5] <- serorev_prev_vals
      i <- i + 1
      next
    }
    if(onlyCalcMFP) {
      curr_mf_prev_vals <- tmpRDSData$mf_prev[c(1,seq(from=183, to=length(tmpRDSData$mf_prev), by=183))]
      total_mf_prev_vals <- length(curr_mf_prev_vals)
      if(i == 1) {
        mf_prev_df <- data.frame(matrix(ncol=5, nrow=total_files*total_mf_prev_vals))
        colnames(mf_prev_df) <- c("ABR", "Ke", "vctr.ctrl.eff", "run_num", "mf_prev")
      
      }
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),1] <- rep(tmpRDSData$ABR, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),2] <- rep(tmpRDSData$Ke, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),3] <- rep(vctr.ctrl.val, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),4] <- rep(run_num, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),5] <- curr_mf_prev_vals
      i <- i + 1
      next
    }
    
    matrix_to_use <- tmpRDSData$ov16_seropositive_matrix
    if(useSerorevert == "finite") {
      matrix_to_use <- tmpRDSData$ov16_seropositive_matrix_serorevert
    } else if (useSerorevert == "instant") {
      matrix_to_use <- tmpRDSData$ov16_seropositive_matrix_serorevert_instant
    }
    
    age <- matrix_to_use[,which_timestep*9-8]
    sex <- ifelse(matrix_to_use[,which_timestep*9-7] == 1, "Male", "Female")
    mf_prev <- matrix_to_use[,which_timestep*9-6]
    ov16_seropos <- matrix_to_use[,which_timestep*9-5]
    ov_l3 <- matrix_to_use[,which_timestep*9-4]
    ov_l4 <- matrix_to_use[,which_timestep*9-3]
    ov_mating_no_mf <- matrix_to_use[,which_timestep*9-2]
    ov_mating_detectable_mf <- matrix_to_use[,which_timestep*9-1]
    ov_mating_any_mf <- matrix_to_use[,which_timestep*9]

    tmpNumRows <- length(age)
    if(i == 1) {
      if(useSerorevert != "none") {
        allOutputs <- data.frame(matrix(ncol=14, nrow=tmpNumRows*total_files))
        colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "Ke", "vctr.ctrl.eff", "sero_type", "run_num")
      } else {
        allOutputs <- data.frame(matrix(ncol=13, nrow=tmpNumRows*total_files))
        colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "Ke", "vctr.ctrl.eff", "run_num")
      }
    }

    startIndex <- 1+tmpNumRows*(i-1)
    endIndex <- tmpNumRows*i
    if(useSerorevert != "none") {
        allOutputs[startIndex:endIndex,-14] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, rep(tmpRDSData$ABR, tmpNumRows), rep(tmpRDSData$Ke, tmpNumRows), rep(vctr.ctrl.val, tmpNumRows),rep("no_infection", tmpNumRows))
        allOutputs[startIndex:endIndex,14] <- run_num
    } else {
        allOutputs[startIndex:endIndex,-13] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, rep(tmpRDSData$ABR, tmpNumRows), rep(tmpRDSData$Ke, tmpNumRows), rep(vctr.ctrl.val, tmpNumRows))
        allOutputs[startIndex:endIndex,13] <- run_num
    }

    i <- i + 1
  }

  print(paste("Overall Pre Treatment MF Prevalence:", mean(mf_prev_vals)))
  print(paste("Avg MDA Years:", mean(mda_vals)))

  if(onlyCalcOv16Trends) {
    ov16Return <- list(ov16_trend_df, ov16_serorevert_trend_df)
    names(ov16Return) <- c("ov16_trend_df", "ov16_serorevert_trend_df")
    return(ov16Return)
  }

  if(onlyCalcMFP) {
    mfpReturn <- list(mf_prev_df, mf_prev_vals, mda_vals, abr_vals)
    names(mfpReturn) <- c("mf_prev_df", 'mf_prev', 'mda_vals', 'abr_vals')
    return(mfpReturn)
  }

  allOutputs <- allOutputs %>% mutate(age_groups = case_when(
                                                age <= 30 ~ ceiling(age/1)*1,
                                                age <= 75 ~ ceiling(age/5)*5,
                                                TRUE ~ 80
                                              ),
                                        age_groups_old = case_when(
                                                age < 5 ~ ceiling(age/1)*1,
                                                age < 11 ~ 7.5,
                                                age < 16 ~ 13,
                                                age < 21 ~ 18,
                                                age < 26 ~ 23,
                                                age < 31 ~ 28,
                                                age < 36 ~ 33,
                                                age < 41 ~ 38,
                                                age < 46 ~ 48,
                                                age < 51 ~ 53,
                                                age < 56 ~ 58,
                                                age < 61 ~ 63,
                                                TRUE ~ 73
                                              ),
                                        age_groups_lsr = case_when(
                                              age < 11 ~ 7.5,
                                              age < 16 ~ 13,
                                              age < 21 ~ 18,
                                              age < 26 ~ 23,
                                              age < 31 ~ 28,
                                              age < 36 ~ 33,
                                              age < 41 ~ 38,
                                              age < 46 ~ 48,
                                              age < 51 ~ 53,
                                              age < 56 ~ 58,
                                              age < 61 ~ 63,
                                              TRUE ~ 73
                                          )
                                      )
  returnVal <- list(allOutputs, mf_prev_vals, abr_vals, mda_vals)
  names(returnVal) <- c("allOutputs", "mf_prev", "abr_vals", "mda_vals")
  return(returnVal)
}

saveRDS(processRCSFiles(onlyCalcMFP = TRUE), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_mfp_trends_data.RDS")

saveRDS(processRCSFiles(onlyCalcOv16Trends=TRUE), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_ov16_trends.RDS")

saveRDS(processRCSFiles(onlyCalcMFP = TRUE), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_mfp_data.RDS")

saveRDS(processRCSFiles(which_timestep=2), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_100_mount_data_start_2015.RDS")

saveRDS(processRCSFiles(which_timestep=2, useSerorevert="finite"), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_finite_seroreversion_data_start_2015.RDS")

saveRDS(processRCSFiles(which_timestep=2, useSerorevert="instant"), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_instant_seroreversion_data_start_2015.RDS")

saveRDS(processRCSFiles(which_timestep=3), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_100_mount_data_mid_2015.RDS")

saveRDS(processRCSFiles(which_timestep=3, useSerorevert="finite"), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_finite_seroreversion_data_mid_2015.RDS")

saveRDS(processRCSFiles(which_timestep=3, useSerorevert="instant"), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_instant_seroreversion_data_mid_2015.RDS")

saveRDS(processRCSFiles(which_timestep=4), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_100_mount_data_end_2015.RDS")

saveRDS(processRCSFiles(which_timestep=4, useSerorevert="finite"), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_finite_seroreversion_data_end_2015.RDS")

saveRDS(processRCSFiles(which_timestep=4, useSerorevert="instant"), "ov16_test_togo/oti_agg_data_final_worm_revert/togo_oti_instant_seroreversion_data_end_2015.RDS")
