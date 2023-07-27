library(dplyr)
library(ggplot2)
library(tidyr)

processRCSFiles <- function () {
  allOutputs <- data.frame()
  #colnames(allOutputs) <- c("age", "sex", "ov16_pos", "mf_prev", "ABR", "run_num")#, "age_pre", "sex_pre", "ov16_pos_pre", "mf_prev_pre", "age_after_year", "sex_after_year", "ov16_pos_after_year", "mf_prev_after_year", "run_num")

  files <- c('ov16_gabon_3/')
  fileToUse <- paste("data/", files[1], sep="")

  i <- 1
  mda_vals <- c()
  abr_vals <- c()
  mf_prev_vals <- c()
  total_files <- length(list.files(fileToUse))
  for (file in list.files(fileToUse)) {
    print(i)
    tmpRDSData <- readRDS(paste(fileToUse, file,sep=""))
    age <- tmpRDSData$ov16_seropositive_matrix[,10]#tmpRDSData$all_infection_burdens[,2]

    sex <- ifelse(tmpRDSData$ov16_seropositive_matrix[,11] == 1, "Male", "Female")#ifelse(tmpRDSData$all_infection_burdens[,3]==1, "Male", "Female")

    mf_prev <- tmpRDSData$ov16_seropositive_matrix[,12]#tmpRDSData$mf_indv_prevalence

    ov16_seropos <- tmpRDSData$ov16_seropositive_matrix[,13]#tmpRDSData$ov16_seropositive

    ov_l3 <- tmpRDSData$ov16_seropositive_matrix[,14]
    ov_l4 <- tmpRDSData$ov16_seropositive_matrix[,15]
    ov_mating_no_mf <- tmpRDSData$ov16_seropositive_matrix[,16]
    ov_mating_detectable_mf <- tmpRDSData$ov16_seropositive_matrix[,17]
    ov_mating_any_mf <- tmpRDSData$ov16_seropositive_matrix[,18]

    mda_vals <- c(mda_vals, tmpRDSData$MDA)
    abr_vals <- c(abr_vals, tmpRDSData$ABR)
    mf_prev_vals <- c(mf_prev_vals, tail(tmpRDSData$mf_prev)[1])

    tmpNumRows <- length(age)
    if(i == 1) {
      allOutputs <- data.frame(matrix(ncol=11, nrow=tmpNumRows*total_files))
      colnames(allOutputs) <- c("age", "sex", "ov16_pos", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "mf_prev", "ABR", "run_num")#, "age_pre", "sex_pre", "ov16_pos_pre", "mf_prev_pre", "age_after_year", "sex_after_year", "ov16_pos_after_year", "mf_prev_after_year", "run_num")
    }

    startIndex <- 1+tmpNumRows*(i-1)
    endIndex <- tmpNumRows*i
    allOutputs[startIndex:endIndex,-11] <- list(age, sex, ov16_seropos, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, mf_prev, tmpRDSData$ABR)#, age_pre, sex_pre, ov16_seropos_pre, mf_prev_pre, age_after_year, sex_after_year, ov_after_year, mf_prev_after_year)
    allOutputs[startIndex:endIndex,11] <- i

    i <- i + 1
  }

  print(paste("Overall Pre Treatment MF Prevalence:", mean(mf_prev_vals)))
  print(paste("Avg MDA Years:", mean(mda_vals)))
  par(mfrow=c(1,2))
  hist(mda_vals)
  print(paste("Avg ABR Value:", mean(abr_vals)))
  hist(abr_vals)
  par(mfrow=c(1,1))

  abr_vs_mf_prev <- ggplot() + geom_boxplot(aes(x=factor(abr_vals), y=mf_prev_vals)) +
    scale_y_continuous(breaks=seq(0, 0.3, 0.05), limits = c(0, 0.3))

  returnVal <- list(allOutputs, mf_prev_vals, abr_vals, mda_vals, abr_vs_mf_prev)
  names(returnVal) <- c("allOutputs", "mf_prev", "abr_vals", "mda_vals", "abr_vs_mf_prev")
  return(returnVal)
}

readGabonData <- function() {
  fileName <- "../Ov16 Data/MW_ML_GAB_data.csv"
  data <- read.csv(fileName, check.names = FALSE)
  data <- data[,1:(length(data)-16)]
  cols_to_keep <- c(2, 8, 9, 10, 14, 16, 17, 20, 21, 33, 34)
  data <- data[,cols_to_keep]
  data <- data %>% filter(Country == "Gabon") %>% mutate(
    ov16_seropos = ifelse(get("Ov16 result") == 1, 1, 0),
    mf_pos = ifelse(Oncho_MF_mean_by_snip > 0, 1, 0),
    sex = ifelse(get("Sex (1=M, 2=F)") == 1, 'Male', 'Female'),
    age_groups = case_when(
      Age <= 30 ~ round(Age/1)*1,
      Age <= 75 ~ round(Age/5)*5,
      TRUE ~ 80
    )
  ) %>% filter(!is.na(ov16_seropos))

  print(paste("Overall Gabon MF Prevalence:", mean(data$mf_pos)))
  print(paste("Overall Gabon Ov16 SeroPrevalence:", mean(data$ov16_seropos)))

  returnVal <- list(data, mean(data$mf_pos))
  names(returnVal) <- c("data", "mf_prev")
  return(returnVal)
}

readGabonDataRDTElisa <- function() {
  fileName <- "../Ov16 Data/Ov16 Gabon ELISA.csv"
  data <- read.csv(fileName) %>% select(1:14) %>% drop_na()
  colnames(data) <- c('id', 'ov16_status_rdt', 'mf_status', 'sample_id', 'inferred_conc', 'region', 'district_name', 'village_name', 'ov16_elisa_plate_num', 'sample_id_ind_id_mismatch', 'ov16_result_elisa', 'ov16_status_elisa', 'age', 'age_group_orig')
  # view duplicate individual id vals, remove them as we don't know which one is correct
  dupIds <- as.numeric(names(which(table(data$id) > 1)))
  data[data$id %in% dupIds,]
  data <- data %>% filter(!(id %in% dupIds)) %>% mutate(
    age_groups = case_when(
      age <= 30 ~ round(age/5)*5,
      age <= 75 ~ round(age/5)*5,
      TRUE ~ 80
    ),
    ov16_status_elisa = case_when(
      ov16_result_elisa == 'Invalid' ~ as.numeric(NaN),
      TRUE ~ as.numeric(ov16_status_elisa)
    )
  )
  print(paste("Total Individuals:", length(unique(data$id))))
  print(paste("Overall Gabon MF Prevalence:", mean(data$mf_status)))
  print(paste("Overall Gabon Ov16 ELISA SeroPrevalence:", mean(data$ov16_status_elisa, na.rm=TRUE)))
  print(paste("Overall Gabon Ov16 RDT SeroPrevalence:", mean(data$ov16_status_rdt)))

  return(data)
}

processActualGabonData <- function(data, justRDT=TRUE, useSex=TRUE) {
  if(justRDT) {
    groupByCols <- c('age_groups')
    if(useSex) {
      groupByCols <- c('age_groups', 'sex')
    }
    data <- data %>% group_by(!!!syms(groupByCols)) %>% summarise(
      ov16_prev = mean(ov16_seropos),
      mf_prev = mean(mf_pos)
    ) %>% as.data.frame() %>% drop_na()

    fit <- lm(ov16_prev ~ log(age_groups+1), data=data)
    data[, c('yhat', 'lb', 'ub')] <- predict(fit, list(age_groups=data$age_groups), interval='prediction')
    # ggplot() + geom_point(aes(x=age_groups, y=ov16_prev), color="red", data=summarize_gabon_data_2) +
    #   geom_line(aes(x=age_groups, y=ov16_prev), color="red", data=summarize_gabon_data_2) +
    #   geom_line(aes(x=age_groups, y=yhat), data=summarize_gabon_data_2)
  } else {
    data <- data %>% group_by(age_groups) %>% summarise(
      ov16_prev_rdt = mean(ov16_status_rdt),
      ov16_prev_elisa = mean(ov16_status_elisa, na.rm=TRUE),
      mf_prev = mean(mf_status)
    ) %>% as.data.frame()

    fit <- lm(ov16_prev_rdt ~ log(age_groups+1), data=data)
    data[, c('yhat_rdt', 'lb_rdt', 'ub_rdt')] <- predict(fit, list(age_groups=data$age_groups), interval='prediction')

    fit <- lm(ov16_prev_elisa ~ log(age_groups+1), data=data)
    data[, c('yhat_elisa', 'lb_elisa', 'ub_elisa')] <- predict(fit, list(age_groups=data$age_groups), interval='prediction')

  }

  return(data)
}

calcSensSpecSeroPrev <- function(run_seropos_data, sens=0.43, spec=0.9998) {
  indv <- length(run_seropos_data)
  prob <- runif(indv)
  new_seropos_data<-rep(0, indv)
  pos <- which(run_seropos_data==1)
  neg <- which(run_seropos_data==0)

  new_seropos_data[pos] <- as.numeric(prob[pos] < sens)
  new_seropos_data[neg] <- as.numeric(prob[neg] > spec)
  return(new_seropos_data)
}

summarizeSimulationData <- function(data, sensSpec=c(1,1), groupBySex=TRUE) {
  df <- data
  sens=sensSpec[1]; spec=sensSpec[2]
  hypothesisNames <- list("Any Worm", "L3 Exposure", "L4-L5 moult", "Mating Worm with any MF", "Mating Worm with detectable MF", "Mating worm")
  if(sens != 1 | spec != 1) {
    df <- df %>% dplyr::group_by(run_num, age, sex, ABR, age_groups) %>%
      dplyr::summarise(ov16_pos=calcSensSpecSeroPrev(ov16_pos, sens, spec),
                       ov16_pos_l3=calcSensSpecSeroPrev(ov16_pos_l3, sens, spec),
                       ov16_pos_l4=calcSensSpecSeroPrev(ov16_pos_l4, sens, spec),
                       ov16_pos_mating_no_mf=calcSensSpecSeroPrev(ov16_pos_mating_no_mf, sens, spec),
                       ov16_pos_mating_detectable_mf=calcSensSpecSeroPrev(ov16_pos_mating_detectable_mf, sens, spec),
                       ov16_pos_mating_any_mf=calcSensSpecSeroPrev(ov16_pos_mating_any_mf, sens, spec),
                       mf_prev=mean(mf_prev)) %>% as.data.frame()
  }
  groupByCols1 <- c("run_num", "age_groups")
  groupByCols2 <- c("age_groups", "Hypothesis")
  hypNamesLoc <- 3:8
  if(groupBySex) {
    groupByCols1 <- c("run_num", "age_groups", "sex")
    groupByCols2 <- c("age_groups", "sex", "Hypothesis")
    hypNamesLoc <- 4:9
  }
  df <- df %>%
    dplyr::group_by(!!!syms(groupByCols1)) %>%
    dplyr::summarise(ov16_any_worm_prev=mean(ov16_pos),
                     ov16_l3_prev=mean(ov16_pos_l3),
                     ov16_l4_prev=mean(ov16_pos_l4),
                     ov16_mating_no_mf_prev=mean(ov16_pos_mating_no_mf),
                     ov16_mating_detectable_mf_prev=mean(ov16_pos_mating_detectable_mf),
                     ov16_mating_any_mf_prev=mean(ov16_pos_mating_any_mf),
                     mf_prev=mean(mf_prev)) %>% as.data.frame()
  colnames(df)[all_of(hypNamesLoc)] <- hypothesisNames
  df <- df %>%
    as.data.frame() %>% pivot_longer(cols=hypNamesLoc, names_to="Hypothesis", values_to = "ov16_prev")

  df <- df %>% dplyr::group_by(!!!syms(groupByCols2)) %>%
    dplyr::summarise(ov16_prev=mean(ov16_prev), mf_prev=mean(mf_prev)) %>% as.data.frame() %>%
    mutate(
      mf_prev = case_when(
        age_groups == 0 & mf_prev > 0 ~ 0,
        TRUE ~ mf_prev
      )
    )

  return(df)
}

makePlots <- function(df1, df2, df3, groupBySex=FALSE, useYhat=FALSE, type='') {
  type <- ifelse(type=='', type, paste('_', type, sep=""))
  plot1 <- ggplot() + geom_smooth(aes(x=age_groups, y=mf_prev), color="blue", data=df1, se=FALSE) +
    geom_line(aes(x=age_groups, y=mf_prev), color="black", data=df2) +
    geom_line(aes(x=age_groups, y=mf_prev), color="red", data=df3)
  if(groupBySex) {
    plot1 <- ggplot() + geom_smooth(aes(x=age_groups, y=mf_prev, linetype=sex), color="blue", data=df1, se=FALSE) +
      geom_line(aes(x=age_groups, y=mf_prev, linetype=sex), color="black", data=df2) +
      geom_line(aes(x=age_groups, y=mf_prev, linetype=sex), color="red", data=df3)
  }
  plot2 <- ggplot()
  plot3 <- ggplot()
  ov16_var <- paste('ov16_prev', type, sep="")
  yhat_var <- paste('yhat', type, sep="")
  lb_var <- paste('lb', type, sep="")
  ub_var <- paste('ub', type, sep="")
  if(useYhat) {
    plot2 <- plot2 + geom_line(aes(x=age_groups, y=get(yhat_var), color="Actual Data", linetype="Actual Data"), data=df1)
    plot3 <- plot3 + geom_line(aes(x=age_groups, y=get(yhat_var), color="Actual Data", linetype="Actual Data"), data=df1)
  } else {
    plot2 <- plot2 + geom_smooth(aes(x=age_groups, y=get(ov16_var), color="Actual Data", linetype="Actual Data"), data=df1, se=FALSE)
    plot3 <- plot3 + geom_smooth(aes(x=age_groups, y=get(ov16_var), color="Actual Data", linetype="Actual Data"), data=df1, se=FALSE)
  }
  plot2 <- plot2 +
    geom_point(aes(x=age_groups, y=get(ov16_var)), color="black", data=df1) +
    geom_line(aes(x=age_groups, y=ov16_prev, color=Hypothesis, linetype=Hypothesis), linewidth=1.1, data=df2) +
    #scale_color_manual(name="Hypotheses", values=c("black", "red", "green", "blue", "pink", "cyan", "yellow")) +
    scale_color_manual(name="Hypotheses", values=c("black", "red", "green", "blue", "red", "green", "blue")) +
    scale_linetype_manual(name="Hypotheses", values=c("solid", "dashed", "dashed", "dashed", "dotted", "dotted", "dotted")) +
    ggtitle("Simulated Data 100% Sens and Spec") +
    xlab("Age") +
    ylab("Ov16 Seroprevalence (%)")

  plot3 <- plot3 +
    geom_point(aes(x=age_groups, y=get(ov16_var)), color="black", data=df1) +
    geom_line(aes(x=age_groups, y=ov16_prev, color=Hypothesis, linetype=Hypothesis), linewidth=1.1, data=df3) +
    #scale_color_manual(name="Hypotheses", values=c("black", "red", "green", "blue", "pink", "cyan", "yellow")) +
    scale_color_manual(name="Hypotheses", values=c("black", "red", "green", "blue", "red", "green", "blue")) +
    scale_linetype_manual(name="Hypotheses", values=c("solid", "dashed", "dashed", "dashed", "dotted", "dotted", "dotted")) +
    ggtitle("Simulated Data Adjusted for 43% Sensitivity and 99.98% Specificity") +
    xlab("Age") +
    ylab("Ov16 Seroprevalence (%)")

  returnVals <- list(plot1, plot2, plot3)
  names(returnVals) <- list("mf", "full", "edited")
  return(returnVals)
}



actualGabonData <- readGabonData()

summarize_gabon_data <- processActualGabonData(actualGabonData$data)

tmpDf <- summarizeSimulationData(summarizeSimulatedGabon)
# values from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6169192/
tmpDf2 <- summarizeSimulationData(summarizeSimulatedGabon, c(0.43,0.9998))

p_males <- makePlots(summarize_gabon_data[summarize_gabon_data$sex == "Male",], tmpDf[tmpDf$sex == "Male",], tmpDf2[tmpDf2$sex == "Male",])
p_males

p_females <- makePlots(summarize_gabon_data[summarize_gabon_data$sex == "Female",], tmpDf[tmpDf$sex == "Female",], tmpDf2[tmpDf2$sex == "Female",])
p_females

# Not split by sex

summarize_gabon_data_2 <- processActualGabonData(actualGabonData$data, useSex = FALSE)


tmpDf3 <- summarizeSimulationData(summarizeSimulatedGabon, groupBySex=FALSE)
# values from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6169192/
tmpDf4 <- summarizeSimulationData(summarizeSimulatedGabon, c(0.43,0.9998), groupBySex=FALSE)

p_comb <- makePlots(summarize_gabon_data_2, tmpDf3, tmpDf4, useYhat=FALSE)
p_comb


actualGabonData_RDTELISA<-readGabonDataRDTElisa()
summarized_gabon_rdt_elisa <- processActualGabonData(t, justRDT=FALSE)

makePlots(summarized_gabon_rdt_elisa, tmpDf3, tmpDf4, useYhat=TRUE, type='elisa')

#sens = TP / TP + FN
#spec = TN / TN + FP
#summarizeSimulatedGabon %>% group_by(run_num, age, sex, ABR, age_groups) %>% summarise(ov16_pos_edited=calcSensSpecSeroPrev(ov16_pos))


# tmpDf <- summarizeSimulatedGabon %>% dplyr::group_by(run_num, age_groups) %>%
#                                                                   dplyr::summarise(ov16_any_worm_prev=mean(ov16_pos), ov16_l3_prev=mean(ov16_pos_l3),
#                                                                                                ov16_l4_prev=mean(ov16_pos_l4), ov16_mating_no_mf_prev=mean(ov16_pos_mating_no_mf),
#                                                                                                ov16_mating_detectable_mf_prev=mean(ov16_pos_mating_detectable_mf),ov16_mating_any_mf_prev=mean(ov16_pos_mating_any_mf),
#                                                                                                mf_prev=mean(mf_prev)) %>% as.data.frame() %>% pivot_longer(cols=3:8, names_to="Hypothesis", values_to = "ov16_prev")
# tmpDf <- tmpDf %>% dplyr::group_by(age_groups, Hypothesis) %>% dplyr::summarise(ov16_prev=mean(ov16_prev), mf_prev=mean(mf_prev)) %>% as.data.frame() %>% mutate(
#   mf_prev = case_when(
#     age_groups == 0 & mf_prev > 0 ~ 0,
#     TRUE ~ mf_prev
#   )
# )
