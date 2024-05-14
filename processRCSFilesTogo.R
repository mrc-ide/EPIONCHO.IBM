library(dplyr)
processRCSFiles <- function (files='/rds/general/user/ar722/home/ov16_test/ov16_output_bassar_mean_age/', which_timestep=2, verbose=TRUE, onlyCalcMFP=FALSE, useSerorevert=FALSE, onlyCalcOv16Trends=FALSE) {
  allOutputs <- data.frame()
  fileToUse <- files#paste("data/", files, sep="")

  i <- 1
  mda_vals <- c()
  abr_vals <- c()
  mf_prev_vals <- c()
  total_files <- length(list.files(fileToUse))
  start_time <- Sys.time() 
  mf_prev_df <- data.frame()
  ov16_trend_df <- matrix(ncol=5, nrow=0)
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
    
    mf_prev_vals <- c(mf_prev_vals, tmpRDSData$mf_prev[(119/(1/366))])
    mda_vals <- c(mda_vals, tmpRDSData$MDA)
    abr_vals <- c(abr_vals, tmpRDSData$ABR)
    vctr.ctrl.val = tmpRDSData$vctr.ctrl.eff
    if(onlyCalcOv16Trends) {
      sero_prev_vals <- tmpRDSData$ov16_seroprevalence[c(1,seq(from=183, to=length(tmpRDSData$ov16_seroprevalence), by=183))]
      total_sero_vals <- length(sero_prev_vals)
      if(i == 1) {
        ov16_trend_df <- matrix(ncol=5, nrow=total_files*total_sero_vals)
        colnames(ov16_trend_df) <- c("ABR", "vctr.ctrl.eff", "Ke", "run_num", "ov16_vals")
      }
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),1] <- rep(tmpRDSData$ABR, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),2] <- rep(vctr.ctrl.val, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),3] <- rep(tmpRDSData$Ke, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),4] <- rep(i,total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),5] <- sero_prev_vals
      i <- i + 1
      next
    }
    if(onlyCalcMFP) {
      if(i == 1) {
        mf_prev_df <- data.frame(matrix(ncol=7, nrow=total_files))
        colnames(mf_prev_df) <- c("ABR", "Ke", "vctr.ctrl.eff", "run_num", "mf_prev_start_2015", "mf_prev_mid_2015", "mf_prev_after_2015")
      }
      mf_prev_df[i,] <- list(tmpRDSData$ABR, tmpRDSData$Ke, vctr.ctrl.val, i, tmpRDSData$mf_prev[(119/(1/366) - 1)], tmpRDSData$mf_prev[(119/(1/366) + 0.5/(1/366))], tmpRDSData$mf_prev[(119/(1/366) + 1/(1/366))])
      i <- i + 1
      next
    }
    
    matrix_to_use <- tmpRDSData$ov16_seropositive_matrix
    if(useSerorevert) {
      matrix_to_use <- tmpRDSData$ov16_seropositive_matrix_serorevert
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
      if(useSerorevert) {
        allOutputs <- data.frame(matrix(ncol=14, nrow=tmpNumRows*total_files))
        colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "Ke", "vctr.ctrl.eff", "sero_type", "run_num")
      } else {
        allOutputs <- data.frame(matrix(ncol=13, nrow=tmpNumRows*total_files))
        colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "Ke", "vctr.ctrl.eff", "run_num")
      }
    }

    startIndex <- 1+tmpNumRows*(i-1)
    endIndex <- tmpNumRows*i
    if(useSerorevert) {
        allOutputs[startIndex:endIndex,-14] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, rep(tmpRDSData$ABR, tmpNumRows), rep(tmpRDSData$Ke, tmpNumRows), rep(vctr.ctrl.val, tmpNumRows),rep("no_infection", tmpNumRows))
        allOutputs[startIndex:endIndex,14] <- i
    } else {
        allOutputs[startIndex:endIndex,-13] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, rep(tmpRDSData$ABR, tmpNumRows), rep(tmpRDSData$Ke, tmpNumRows), rep(vctr.ctrl.val, tmpNumRows))
        allOutputs[startIndex:endIndex,13] <- i
    }

    i <- i + 1
  }

  print(paste("Overall Pre Treatment MF Prevalence:", mean(mf_prev_vals)))
  print(paste("Avg MDA Years:", mean(mda_vals)))

  if(onlyCalcOv16Trends) {
    ov16Return <- list(ov16_trend_df)
    names(ov16Return) <- c("ov16_trend_df")
    return(ov16Return)
  }

  if(onlyCalcMFP) {
    mfpReturn <- list(mf_prev_df, mf_prev_vals, mda_vals, abr_vals)
    names(mfpReturn) <- c("mf_prev_df", 'mf_prev', 'mda_vals', 'abr_vals')
    return(mfpReturn)
  }

  #abr_vs_mf_prev <- ggplot() + geom_boxplot(aes(x=factor(abr_vals), y=mf_prev_vals)) +
  #  scale_y_continuous(breaks=seq(0, 0.3, 0.05), limits = c(0, 0.3))

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
  returnVal <- list(allOutputs, mf_prev_vals, abr_vals, mda_vals)#, abr_vs_mf_prev)
  names(returnVal) <- c("allOutputs", "mf_prev", "abr_vals", "mda_vals")#, "abr_vs_mf_prev")
  return(returnVal)
}

saveRDS(processRCSFiles(onlyCalcMFP = TRUE), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_mfp_data.RDS")

saveRDS(processRCSFiles(which_timestep=2), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_100_mount_data_start_2015.RDS")

saveRDS(processRCSFiles(which_timestep=2, useSerorevert=TRUE), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_seroreversion_data_start_2015.RDS")

saveRDS(processRCSFiles(which_timestep=3), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_100_mount_data_mid_2015.RDS")

saveRDS(processRCSFiles(which_timestep=3, useSerorevert=TRUE), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_seroreversion_data_mid_2015.RDS")

saveRDS(processRCSFiles(which_timestep=4), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_100_mount_data_end_2015.RDS")

saveRDS(processRCSFiles(which_timestep=4, useSerorevert=TRUE), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_seroreversion_data_end_2015.RDS")

saveRDS(processRCSFiles(onlyCalcOv16Trends=TRUE), "/rds/general/user/ar722/home/ov16_test/bassar_agg_data/togo_bassar_ov16_trends.RDS")
