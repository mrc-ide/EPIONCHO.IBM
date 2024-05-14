library(dplyr)
processRCSFiles <- function (files='/rds/general/user/ar722/home/ov16_test/ov16_output/', verbose=TRUE, onlyCalcMFP=FALSE, useSerorevert=FALSE, onlyCalcOv16Trends=FALSE) {
  allOutputs <- data.frame()
  fileToUse <- files#paste("data/", files, sep="")

  i <- 1
  mda_vals <- c()
  abr_vals <- c()
  mf_prev_vals <- c()
  total_files <- length(list.files(fileToUse))
  start_time <- Sys.time() 
  mf_prev_df <- matrix(ncol=4, nrow=0)#data.frame()
  ov16_trend_df <- matrix(ncol=4, nrow=0)
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
    
    mf_prev_vals <- c(mf_prev_vals, tail(tmpRDSData$mf_prev)[1])
    mda_vals <- c(mda_vals, tmpRDSData$MDA)
    abr_vals <- c(abr_vals, tmpRDSData$ABR)
    if(onlyCalcOv16Trends) {
      sero_prev_vals <- tmpRDSData$ov16_seroprevalence[c(1,seq(from=183, to=length(tmpRDSData$ov16_seroprevalence), by=183))]
      total_sero_vals <- length(sero_prev_vals)
      if(i == 1) {
        ov16_trend_df <- matrix(ncol=4, nrow=total_files*length(sero_prev_vals))
        colnames(ov16_trend_df) <- c("ABR", "Ke", "run_num", "ov16_vals")
      }
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),1] <- rep(tmpRDSData$ABR, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),2] <- rep(tmpRDSData$Ke, total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),3] <- rep(i,total_sero_vals)
      ov16_trend_df[(1+(total_sero_vals*(i-1))):(i*total_sero_vals),4] <- sero_prev_vals
      i <- i + 1
      next
    }
    if(onlyCalcMFP) {
      curr_mf_prev_vals <- tmpRDSData$mf_prev[c(1,seq(from=366, to=length(tmpRDSData$mf_prev), by=366))]
      total_mf_prev_vals <- length(curr_mf_prev_vals)
      if(i == 1) {
        mf_prev_df <- matrix(ncol=4, nrow=total_files*total_mf_prev_vals)
        colnames(mf_prev_df) <- c("ABR", "Ke", "run_num", "mf_prev")
      }
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),1] <- rep(tmpRDSData$ABR, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),2] <- rep(tmpRDSData$Ke, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),3] <- rep(i, total_mf_prev_vals)
      mf_prev_df[(1+(total_mf_prev_vals*(i-1))):(i*total_mf_prev_vals),4] <- curr_mf_prev_vals
      i <- i + 1
      next
    }
    
    matrix_to_use <- tmpRDSData$ov16_seropositive_matrix
    if(useSerorevert) {
      matrix_to_use <- tmpRDSData$ov16_seropositive_matrix_serorevert
    }
    
    age <- matrix_to_use[,1]
    sex <- ifelse(matrix_to_use[,2] == 1, "Male", "Female")
    mf_prev <- matrix_to_use[,3]
    ov16_seropos <- matrix_to_use[,4]
    ov_l3 <- matrix_to_use[,5]
    ov_l4 <- matrix_to_use[,6]
    ov_mating_no_mf <- matrix_to_use[,7]
    ov_mating_detectable_mf <- matrix_to_use[,8]
    ov_mating_any_mf <- matrix_to_use[,9]

    tmpNumRows <- length(age)
    if(i == 1) {
      if(useSerorevert) {
        allOutputs <- data.frame(matrix(ncol=13, nrow=tmpNumRows*total_files))
        colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "Ke", "sero_type", "run_num")
      } else {
        allOutputs <- data.frame(matrix(ncol=12, nrow=tmpNumRows*total_files))
        colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "Ke", "run_num")
      }
    }

    startIndex <- 1+tmpNumRows*(i-1)
    endIndex <- tmpNumRows*i
    if(useSerorevert) {
        allOutputs[startIndex:endIndex,-13] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, rep(tmpRDSData$ABR, tmpNumRows), rep(tmpRDSData$Ke, tmpNumRows), rep(tmpRDSData$sero_type, tmpNumRows))
        allOutputs[startIndex:endIndex,13] <- i
    } else {
        allOutputs[startIndex:endIndex,-12] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, rep(tmpRDSData$ABR, tmpNumRows), rep(tmpRDSData$Ke, tmpNumRows))
        allOutputs[startIndex:endIndex,12] <- i
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

  allOutputs <- allOutputs %>% mutate(age_groups_old = case_when(
                                                age <= 30 ~ ceiling(age/1)*1,
                                                age <= 75 ~ ceiling(age/5)*5,
                                                TRUE ~ 80
                                              ),
                                      age_groups = case_when(
                                                ceiling(age*5)/5 == 0 ~ 0,
                                                age <= 75 ~ ceiling(age/5)*5 - 2.5,
                                                TRUE ~ 77.5
                                              ),
                                      age_groups_floor = case_when(
                                                age <= 75 ~ ceiling(age/5)*5,
                                                TRUE ~ 80
                                              )
                                      )
  returnVal <- list(allOutputs, mf_prev_vals, abr_vals, mda_vals)#, abr_vs_mf_prev)
  names(returnVal) <- c("allOutputs", "mf_prev", "abr_vals", "mda_vals")#, "abr_vs_mf_prev")
  return(returnVal)
}

saveRDS(processRCSFiles(onlyCalcMFP = TRUE), "/rds/general/user/ar722/home/ov16_test/gabon_agg_data/gabon_mfp_abr_all_age_data.RDS")#gabon_mfp_data.RDS")

#saveRDS(processRCSFiles(onlyCalcOv16Trends=TRUE), "/rds/general/user/ar722/home/ov16_test/gabon_agg_data/gabon_ov16_trends.RDS")

#saveRDS(processRCSFiles(), "/rds/general/user/ar722/home/ov16_test/gabon_agg_data/gabon_mfp_abr_data.RDS")#gabon_100_mount_data.RDS")

#saveRDS(processRCSFiles(useSerorevert=TRUE), "/rds/general/user/ar722/home/ov16_test/gabon_agg_data/gabon_seroreversion_data.RDS")