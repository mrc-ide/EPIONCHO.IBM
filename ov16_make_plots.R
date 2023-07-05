# fit to gabon - no treatment
# Mali - bakoye - pre-control, history of intervention (20-24 years) - RDT losing sens
# TOGO - have MF data - identify the villages, find pre-control, history of intervention, see what we have when we look at ov16
# change ABR to match MFP
# make a presentation
# Sampling biting
# sigmoidicity
# catalytic model

library(dplyr)
library(ggplot2)

processRCSFiles <- function () {
  #allOutputs <- data.frame(matrix(ncol=25))
  #colnames(allOutputs) <- c("age", "sex", "ov16_pos", "mf_prev", "age_pre", "sex_pre", "ov16_pos_pre", "mf_prev_pre", "age3", "age4", "age5", "age6", "sex3", "sex4", "sex5", "sex6",  "ov16_pos3", "ov16_pos4", "ov16_pos5", "ov16_pos6", "mf_prev3", "mf_prev4", "mf_prev5", "mf_prev6", "run_num")

  allOutputs <- data.frame()#data.frame(matrix(ncol=19))
  #colnames(allOutputs) <- c("age", "sex", "ov16_pos", "mf_prev", "age_pre", "sex_pre", "ov16_pos_pre", "mf_prev_pre", "ov_l3", "ov_l4", "ov_mating_no_mf", "ov_mating_detectable_mf", "ov_mating_any_mf", "ov_l3_pre", "ov_l4_pre", "ov_mating_no_mf_pre", "ov_mating_detectable_mf_pre", "ov_mating_any_mf_pre", "run_num")

  files <- c('ov16_all_hypotheses_3/', 'ov16_output_any_worm_gamm_7/', 'ov16_output_any_worm_gamm_6/', 'ov16_output_any_worm_44/', 'ov16_output_l3/', 'ov16_output_l3_2/', 'ov16_output_l3_3/')
  fileToUse <- paste("data/", files[1], sep="")

  i <- 1
  mda_vals <- c()
  abr_vals <- c()
  mf_prev_vals <- c()
  total_files <- length(list.files(fileToUse))
  for (file in list.files(fileToUse)) {
    print(i)
    #c(treat.start-1, treat.stop, treat.stop+(1/DT)+1)
    tmpRDSData <- readRDS(paste(fileToUse, file,sep=""))
    age <- tmpRDSData$ov16_seropositive_matrix[,10]

    age_pre <- tmpRDSData$ov16_seropositive_matrix[,1]

    sex <- ifelse(tmpRDSData$ov16_seropositive_matrix[,11] == 1, "Male", "Female")

    sex_pre <- ifelse(tmpRDSData$ov16_seropositive_matrix[,2]==1, "Male", "Female")

    mf_prev <- tmpRDSData$ov16_seropositive_matrix[,12]

    mf_prev_pre <- tmpRDSData$ov16_seropositive_matrix[,3]

    ov16_seropos <- tmpRDSData$ov16_seropositive_matrix[,13]

    ov16_seropos_pre <- tmpRDSData$ov16_seropositive_matrix[,4]

    # age3 <- tmpRDSData$ov16_seropositive_matrix[,7]
    # age4 <- tmpRDSData$ov16_seropositive_matrix[,10]
    # age5 <- tmpRDSData$ov16_seropositive_matrix[,13]
    # age6 <- tmpRDSData$ov16_seropositive_matrix[,16]
    #
    # sex3 <- ifelse(tmpRDSData$ov16_seropositive_matrix[,8]==1, 'Male', 'Female')
    # sex4 <- ifelse(tmpRDSData$ov16_seropositive_matrix[,11]==1, 'Male', 'Female')
    # sex5 <- ifelse(tmpRDSData$ov16_seropositive_matrix[,14]==1, 'Male', 'Female')
    # sex6 <- ifelse(tmpRDSData$ov16_seropositive_matrix[,17]==1, 'Male', 'Female')
    #
    # ov3 <- tmpRDSData$ov16_seropositive_matrix[,9]
    # ov4 <- tmpRDSData$ov16_seropositive_matrix[,12]
    # ov5 <- tmpRDSData$ov16_seropositive_matrix[,15]
    # ov6 <- tmpRDSData$ov16_seropositive_matrix[,18]
    #
    # mf_prev3 <- tmpRDSData$mf_indv_prev_matrix[,9]
    # mf_prev4 <- tmpRDSData$mf_indv_prev_matrix[,12]
    # mf_prev5 <- tmpRDSData$mf_indv_prev_matrix[,15]
    # mf_prev6 <- tmpRDSData$mf_indv_prev_matrix[,18]

    ov_l3 <- tmpRDSData$ov16_seropositive_matrix[,14]
    ov_l4 <- tmpRDSData$ov16_seropositive_matrix[,15]
    ov_mating_no_mf <- tmpRDSData$ov16_seropositive_matrix[,16]
    ov_mating_detectable_mf <- tmpRDSData$ov16_seropositive_matrix[,17]
    ov_mating_any_mf <- tmpRDSData$ov16_seropositive_matrix[,18]

    ov_l3_pre <- tmpRDSData$ov16_seropositive_matrix[,5]
    ov_l4_pre <- tmpRDSData$ov16_seropositive_matrix[,6]
    ov_mating_no_mf_pre <- tmpRDSData$ov16_seropositive_matrix[,7]
    ov_mating_detectable_mf_pre <- tmpRDSData$ov16_seropositive_matrix[,8]
    ov_mating_any_mf_pre <- tmpRDSData$ov16_seropositive_matrix[,9]



    tmpNumRows <- length(age)
    if(i == 1) {
      allOutputs <- data.frame(matrix(ncol=19, nrow=tmpNumRows*total_files))
      colnames(allOutputs) <- c("age", "sex", "ov16_pos", "mf_prev", "age_pre", "sex_pre", "ov16_pos_pre", "mf_prev_pre", "ov16_pos_l3", "ov16_pos_l4", "ov16_pos_mating_no_mf", "ov16_pos_mating_detectable_mf", "ov16_pos_mating_any_mf", "ov16_pos_l3_pre", "ov16_pos_l4_pre", "ov16_pos_mating_no_mf_pre", "ov16_pos_mating_detectable_mf_pre", "ov16_pos_mating_any_mf_pre", "run_num")
    }

    startIndex <- 1+tmpNumRows*(i-1)
    endIndex <- tmpNumRows*i
    # allOutputs[startIndex:endIndex,-25] <- list(age, sex, ov16_seropos, mf_prev, age_pre, sex_pre, ov16_seropos_pre, mf_prev_pre, age3, age4, age5, age6, sex3, sex4, sex5, sex6, ov3, ov4, ov5, ov6, mf_prev3, mf_prev4, mf_prev5, mf_prev6)
    # allOutputs[startIndex:endIndex,25] <- i
    allOutputs[startIndex:endIndex,-19] <- list(age, sex, ov16_seropos, mf_prev, age_pre, sex_pre, ov16_seropos_pre, mf_prev_pre, ov_l3, ov_l4, ov_mating_no_mf, ov_mating_detectable_mf, ov_mating_any_mf, ov_l3_pre, ov_l4_pre, ov_mating_no_mf_pre, ov_mating_detectable_mf_pre, ov_mating_any_mf_pre)
    allOutputs[startIndex:endIndex,19] <- i

    mda_vals <- c(mda_vals, tmpRDSData$MDA)

    abr_vals <- c(abr_vals, tmpRDSData$ABR)

    mf_prev_vals <- c(mf_prev_vals, tmpRDSData$mf_prev[36600])
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
     scale_y_continuous(breaks=seq(0.6, 0.7, 0.01), limits = c(0.6, 0.7))
  returnVals <- list("allOutputs"=allOutputs, "abr_vs_mf_graph"=abr_vs_mf_prev, "mf_prev_vals"=mf_prev_vals)
  return(returnVals)
}

tmpAllOutputs <- processRCSFiles()

tmpAllOutputs$abr_vs_mf_graph


# viz

age_var <- 'age'
sex_var <- 'sex'
ov_var <- 'ov16_pos_l4'
ov_pre_var <- 'ov16_pos_l4_pre'
mf_var <- 'mf_prev'

allOutputs <- tmpAllOutputs$allOutputs %>% mutate(age_groups = #cut(get(age_var), c(seq(0, 30, 0.5), seq(35, 80, 5)), include.lowest = TRUE, labels=c(seq(0.5, 30, 0.5), seq(35, 80, 5))),
                                         case_when(
                                       get(age_var) <= 30 ~ ceiling(get(age_var)/1)*1,
                                       get(age_var) <= 75 ~ ceiling(get(age_var)/5)*5,
                                       TRUE ~ 80
                                       ),
                                       age_groups_pre = #cut(get(age_var), c(seq(0, 30, 0.5), seq(35, 80, 5)), include.lowest = TRUE, labels=c(seq(0.5, 30, 0.5), seq(35, 80, 5)))
                                                            case_when(
                                         age_pre <= 30 ~ ceiling(age_pre/1)*1,
                                         age_pre <= 75 ~ ceiling(age_pre/5)*5,
                                         TRUE ~ 80
                                       )
                                       )
allOutputs$age_groups <- as.numeric(as.character(allOutputs$age_groups))
allOutputs$age_groups_pre <- as.numeric(as.character(allOutputs$age_groups_pre))


tmpDf <- allOutputs %>% dplyr::group_by(run_num, age_groups, get(sex_var)) %>% dplyr::summarise(ov16_prev=mean(get(ov_var)), mf_prev=mean(get(mf_var))) %>% as.data.frame()
colnames(tmpDf) <- c('run_num', 'age_groups', 'sex', 'ov16_prev', 'mf_prev')
tmpDf <- tmpDf %>% dplyr::group_by(age_groups, sex) %>% dplyr::summarise(ov16_prev=mean(ov16_prev), mf_prev=mean(mf_prev)) %>% as.data.frame()
tmpDf[1,4] <- 0
#tmpDf <- allOutputs %>% dplyr::group_by(age_groups, get(sex_var)) %>% dplyr::summarise(ov16_prev=mean(get(ov_var)), mf_prev=mean(get(mf_var))) %>% as.data.frame()
#colnames(tmpDf) <- c('age_groups', 'sex', 'ov16_prev', 'mf_prev')


tmpDf2 <- allOutputs %>% dplyr::group_by(run_num, age_groups_pre, sex_pre) %>% dplyr::summarise(ov16_prev_pre=mean(get(ov_pre_var)), mf_prev_pre=mean(mf_prev_pre)) %>% as.data.frame() %>%
  dplyr::group_by(age_groups_pre, sex_pre) %>% dplyr::summarise(ov16_prev_pre=mean(ov16_prev_pre), mf_prev_pre=mean(mf_prev_pre)) %>% as.data.frame()
tmpDf2[1:2,4] <- 0

#tmpDf2 <- allOutputs %>% dplyr::group_by(age_groups_pre, sex_pre) %>% dplyr::summarise(ov16_prev_pre=mean(ov16_pos_pre), mf_prev_pre=mean(mf_prev_pre)) %>% as.data.frame()


ov16_graph <- ggplot() +
  geom_line(aes(x=age_groups_pre, y=ov16_prev_pre*100, color="Pre Treatment", linetype=sex_pre), linewidth=1.2, data=tmpDf2) +
  geom_line(aes(x=age_groups, y=ov16_prev*100, color='Post Treatment', linetype=sex), linewidth=1.2, data=tmpDf) +
  xlab("Age") +
  ylab("OV16 Seroprevalence (%)") +
  scale_x_continuous(breaks=c(seq(0,10,5), seq(20, 80, 10))) +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.x = element_line( linewidth=.15, color=rgb(0, 0, 0, 20, maxColorValue=255)),
    panel.grid.major.y = element_line( linewidth=.15, color=rgb(0, 0, 0, 20, maxColorValue=255)),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    text=element_text(size=15),
    axis.text = element_text(size=20),
    axis.title= element_text(size=15)) +
  scale_y_continuous(breaks=seq(0,101,20), limits=c(0, 100)) +
  scale_linetype_manual(values=c("dashed", "dotted")) +
  ggtitle("Ov16 Seroprevalence in presence of Mating Worm Pair") +
  scale_color_manual(values=c("red", "black"))

ov16_graph

#ggsave("ov16_any_adult_worm_graph_gamma_dot4.png", ov16_graph, width=5000, height = 4000, units="px", dpi=600)
ggsave("ov16_mating_worm_gamma_dot4.png", ov16_graph, width=5000, height = 4000, units="px", dpi=600)

mf_prev_graph <- ggplot()  +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.x = element_line( linewidth=.1, color=rgb(0, 0, 0, 20, maxColorValue=255)),
    panel.grid.major.y = element_line( linewidth=.1, color=rgb(0, 0, 0, 20, maxColorValue=255)),
    panel.background = element_rect(fill = 'white', colour = 'black')
  ) +
  geom_line(aes(x=age_groups_pre, y=mf_prev_pre*100, color='Pre Treatment', linetype=sex_pre), linewidth=1.2, data=tmpDf2) +
  geom_line(aes(x=age_groups, y=mf_prev*100, color='Post Treatment', linetype=sex), linewidth=1.2, data=tmpDf) +
  scale_x_continuous(seq(0, 81, 5), name="Age") +
  scale_y_continuous(seq(0, 101, 20), limits=c(0, 100), name="mf prevalence (%)") +
  scale_linetype_manual(values=c("dashed", "dotted")) +
  scale_color_manual(values=c("red", "black"))
mf_prev_graph

mf_nested <- ggplot() +
  geom_line(aes(x=age_groups, y=mf_prev*100, color='Post Treatment', linetype=sex),  linewidth=1.2, data=tmpDf) +
  scale_x_continuous(seq(0, 81, 5), name="") +
  scale_y_continuous(seq(0, 5, 1), limits=c(0, 5), name="") +
  scale_linetype_manual(values=c("dashed", "dotted")) +
  scale_color_manual(values=c("red", "black")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.x = element_line( linewidth=.1, color=rgb(0, 0, 0, 20, maxColorValue=255)),
    panel.grid.major.y = element_line( linewidth=.1, color=rgb(0, 0, 0, 20, maxColorValue=255)),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    legend.position='none',
    axis.text = element_text(size=7)
  )

mf_prev_graph <- mf_prev_graph +
  geom_rect(aes(xmin=14, xmax=70, ymin=15, ymax=50), alpha=0.3, color="black") +
  annotation_custom(ggplotGrob(mf_nested), xmin = 10, xmax = 70, ymin = 10, ymax = 50) +
  ggtitle("MF Prevalence vs Age")

mf_prev_graph


ggsave("mf_prev_graph_dot4.png", mf_prev_graph, width=3500, height = 2000, units="px", dpi=600)


grid_graph <- grid.arrange(mf_prev_graph, ov16_graph, ncol=2)
grid_graph

ggsave("both_graphs.png", grid_graph, width=7000, height = 2000, units="px", dpi=600)

# ---- Togo and Gabon

readGabonData <- function() {
  fileName <- "../Ov16 Data/MW_ML_GAB_data.csv"
  data <- read.csv(fileName, check.names = FALSE)
  data <- data[,1:(length(data)-16)]
  cols_to_keep <- c(2, 8, 9, 10, 14, 16, 17, 20, 21, 33, 34)
  data <- data[,cols_to_keep]
  print(colnames(data))
  data <- data %>% filter(Country == "Gabon") %>% mutate(
    ov16_seropos = ifelse(get("Ov16 result") == 1, 1, 0),
    mf_pos = ifelse(Oncho_MF_mean_by_snip > 0, 1, 0),
    sex = ifelse(get("Sex (1=M, 2=F)") == 1, 'Male', 'Female'),
    age_groups = case_when(
      Age <= 30 ~ ceiling(Age/1)*1,
      Age <= 70 ~ ceiling(Age/10)*10,
      TRUE ~ 80
    )
  ) %>% filter(!is.na(ov16_seropos))
  print(length(unique(data$IndividualNo)))

  return(data)
}

gabonData <- readGabonData() %>% group_by(age_groups, sex) %>% summarise(
  ov16_prev = mean(ov16_seropos),
  mf_prev = mean(mf_pos)
)

gabonData %>% ggplot() + geom_line(aes(x=age_groups, y=mf_prev*100, color=sex)) +
  ylim(0, 100)

gabonData %>% ggplot() + geom_point(aes(x=age_groups, y=ov16_prev*100, color=sex)) +
  geom_smooth(aes(x=age_groups, y=ov16_prev*100, color=sex), se=FALSE) +
  ylim(0, 50) +
  ylab("Ov16 Seroprevalence %") +
  xlab("Age Groups") +
  ggtitle("Ov16 Seroprevalence vs Age")

togoDataFunc <- function() {
  # 2 years each side except for 7.5 (2.5y gap) and 70.5 (9.5y gap)
  ageGroups <- c(7.5, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 70.5)
  male_pop <- c(96, 106, 31, 32, 58, 45, 42, 45, 37, 32, 29, 42)
  female_pop <- c(76, 106, 64, 86, 137, 90, 63, 44, 53, 30, 32, 38)
  mf_pos_male <- c(1, 1, 0, 2, 9, 6, 6, 4, 3, 2, 1, 1)
  mf_pos_female <- c(3, 3, 1, 4, 7, 7, 8, 3, 3, 0, 3, 4)
  data <- data.frame(age_groups=ageGroups, male_pop=male_pop, female_pop=female_pop, mf_pos_male=mf_pos_male, mf_pos_female=mf_pos_female) %>%
    mutate(mf_prev_male = mf_pos_male/male_pop,
           mf_prev_female=mf_pos_female/female_pop)
  return(data)
}

togoOv16DataFunc <- function() {
  # 2 years each side except for 7.5 (2.5y gap) and 75.5 (4.5y gap)
  ageGroups <- c(7.5, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 63, 68, 75.5)
  participants <- c(87, 76, 41, 44, 75, 58, 40, 27, 43, 24, 30, 13, 16, 2)
  ov16_pos <- c(13, 12, 7, 18, 31, 22, 20, 11, 15, 15, 17, 6, 9, 2)

  ageGroups <- c(0, 7.5, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 63, 73)
  participants <- c(1, 87, 76, 41, 44, 75, 58, 40, 27, 43, 24, 30, 13, 18)
  ov16_pos <- c(0, 13, 12, 7, 18, 31, 22, 20, 11, 15, 15, 17, 6, 11)
  data <- data.frame(age_groups=ageGroups, participants=participants, ov16_seropos=ov16_pos) %>%
    mutate(ov16_seroprev = ov16_pos/participants)
  return(data)
}

togoData <- togoDataFunc()

togoData %>% ggplot() + geom_line(aes(x=age_groups, y=mf_prev_male*100)) +
  geom_line(aes(x=age_groups, y=mf_prev_female*100), color="red") +
  ylim(0, 50) +
  xlab("Age") +
  ylab("MFP")

togoOv16Data <- togoOv16DataFunc()

togoOv16Data %>% ggplot() + geom_point(aes(x=age_groups, y=ov16_seroprev*100)) +
  geom_smooth(aes(x=age_groups, y=ov16_seroprev*100), color='red') +
  ylim(0, 100) +
  xlab("Age") +
  ylab("Ov16 Seroprevalence")
