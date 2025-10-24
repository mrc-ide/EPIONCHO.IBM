#####Minimal square errors#####
setwd("set to location of datafile") 
#download the dataset with the village level microfilarial prevalence of Togo, and the simulations from EPIONCHO-IBM with vector control and MDA
#The analysis is conducted per region of Togo,  Special Intervention Zone (SIZ) status and if the village baseline endemicity is known/unknown.
#open the following libraries
library(ggpubr)
library(patchwork) 
library(dplyr)

###SAVANES SIZ WITHOUT BASELINE####
#unknown
SavanesYesEND1 <- subset(SavanesYes5, VILNOM %in% c("FAREO","N'KPE 2", "DIMONGUE","GNANGBANDI","GBEMBA-BAS","LEGBANDO","NAGBAKOU","TOUTIONGA")) #removed SADORI, as it belongs to hyper at baseline
SavanesYesEND2 <- subset(SavanesYes5, VILNOM %in% c("MANTCHE","KPINKPARPAK","NATOUNDJENGA", "NALOGBANDI","TCHITCHILINGA","TOGOU"))
SavanesYesEND3 <- subset(SavanesYes5, VILNOM %in% c("POPORKOU","NAMBOSSI","KOULAGNIERE","NATOUNKPARGOU","NATOUNDJENGA","NASSIELE","TCHRI","KOUKOUMBOU","MOUKAGA","KPATIBORI","BONSOUGOU", "BOUTCHAKOU","DJANDJATIE","KPINTIDJOUAGA","NABOLI","PANCERYS","PANSIERI","SIMBO","SOUGTANGOU","TCHOUNTCHONGA","TCHRI","YIYINGOU"))
SavanesYesEND4 <- subset(SavanesYes5, VILNOM %in% c("TCHRI","KOUKOUMBOU"))
SavanesYesEND1$ENDE <- "Hypoendemic trends"
SavanesYesEND2$ENDE <- "Mesoendemic trends"
SavanesYesEND3N <- subset(SavanesYesEND3, VILNOM %in% c("LOKPANO","NATOUNKPARGOU","KOULAGNIERE","NATOUNDJENGA","NASSIELE","MOUKAGA","KPATIBORI","BONSOUGOU", "BOUTCHAKOU","TCHRI","KOUKOUMBOU","DJANDJATIE","KPINTIDJOUAGA","NABOLI","PANCERYS","PANSIERI","SIMBO","SOUGTANGOU","TCHOUNTCHONGA","YIYINGOU"))
SavanesYesEND3N$ENDE <- "Hyperendemic trends"

#now the known endemicity levels
SavanesYes1B #hypo 100% VC
PANGA_dataset <- subset(SavanesYes2, VILNOM == "PANGA") #meso
# Remove "PANGA" from the original dataset
SavanesYes2100 <- subset(SavanesYes2, VILNOM != "PANGA") #meso 100% VC

#HYPOENDEMIC
list1 <- lapply(means.examplesSSOB, `[[`, 1)
list11 <- lapply(means.examplesSSRB, `[[`, 1)
list111 <- lapply(means.examplesSSPB, `[[`, 1)
#meso
list2 <- lapply(means.examplesSSOB, `[[`, 2)
list22 <- lapply(means.examplesSSRB, `[[`, 2)
list222 <- lapply(means.examplesSSPB, `[[`, 2)
#hyper
list3 <- lapply(means.examplesSSOB, `[[`, 3)
list33 <- lapply(means.examplesSSRB, `[[`, 3)
list333 <- lapply(means.examplesSSPB, `[[`, 3)
#holo
list4 <- lapply(means.examplesSSOB, `[[`, 4)
list44 <- lapply(means.examplesSSRB, `[[`, 4)
list444 <- lapply(means.examplesSSPB, `[[`, 4)
#hypo with 100% VC
list5 <- lapply(means.examplesSSO100, `[[`, 1)
list55 <- lapply(means.examplesSSR100, `[[`, 1)
list555 <- lapply(means.examplesSSP100, `[[`, 1)
#meso with 100% VC
list6 <- lapply(means.examplesSSO100, `[[`, 2)
list66 <- lapply(means.examplesSSR100, `[[`, 2)
list666 <- lapply(means.examplesSSP100, `[[`, 2)
#hyper with 100% VC
list7 <- lapply(means.examplesSSO100, `[[`, 3)
list77 <- lapply(means.examplesSSR100, `[[`, 3)
list777 <- lapply(means.examplesSSP100, `[[`, 3)
#holo with 100% VC
list8 <- lapply(means.examplesSSO100, `[[`, 4)
list88 <- lapply(means.examplesSSR100, `[[`, 4)
list888 <- lapply(means.examplesSSP100, `[[`, 4)
# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year-1) #minus a year for Savanes
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]]),
  hypoendemic100 = list(list5[[2]], list55[[2]], list555[[2]]),
  mesoendemic100 = list(list6[[2]], list66[[2]], list666[[2]]),
  hyperendemic100 = list(list7[[2]], list77[[2]], list777[[2]]),
  holoendemic100 = list(list8[[2]], list88[[2]], list888[[2]])
)

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = SavanesYesEND1,
  mesoendemic = SavanesYesEND2,
  hyperendemic = SavanesYesEND3,
  holoendemic = SavanesYesEND4
)

# Observed known baseline
known_datasets <- list(
  hypoendemic = SavanesYes1B,#no time constraint for hypoendemicity, as there are no recent surveys
  mesoendemic = SavanesYes2100[SavanesYes2100$Survey_year >= 1989, ], #PANGA_dataset or SavanesYes2100
  hypoendemic100 = SavanesYes1B,#no time constraint for hypoendemicity, as there are no recent surveys
  mesoendemic100 = SavanesYes2100[SavanesYes2100$Survey_year >= 1989, ]#PANGA_dataset or SavanesYes2100
)

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)
#rerum mesoendemic with Panga datest (  mesoendemic1 = PANGA_dataset) which follows mesoendemicity instead of SavanesYes2100 that follows mesoendemicity with 100% VC efficacy




###SAVANES no-SIZ WITHOUT BASELINE####
#unknown
#hypoendemic and mesoendemic without biannual (T?ne and Kpendjal prefecture)
SavanesNoEND1 <- subset(SavanesNo5, VILNOM %in% c("LOUGOU","TINNOGO")) 
SavanesNoEND1$ENDE <- "Hypoendemic trends"
#now the known endemicity levels
SavanesNo2 #meso
SavanesHyperNo <- SavanesNo3
# Extract the "SAMOMONI" row from SavanesNo5
samomoni_row <- SavanesNo5 %>% filter(VILNOM == "SAMOMONI")
# Append the "SAMOMONI" row to the new dataset
SavanesHyperNo <- bind_rows(SavanesHyperNo, samomoni_row)
#SavanesHyperNo <- data.frame(samomoni_row)
# Create a new row based on "SAMOMONI" row with specified changes
new_samomoni_row <- samomoni_row %>%
  mutate(Survey_year = 1974, 
         PREV_cr = 60.2, 
         PREV_st = NA, 
         CILower = 51.1, 
         CIUpper = 68.6)
# Append the new row to the dataset
SavanesHyperNo <- bind_rows(SavanesHyperNo, new_samomoni_row)
SavanesHyperNo <- SavanesHyperNo %>%
  mutate(ENDE = ifelse(VILNOM == "SAMOMONI", "Hyperendemic", ENDE))
#HYPOENDEMIC
list1 <- lapply(means.examplesSNSO, `[[`, 1)
list11 <- lapply(means.examplesSNSR, `[[`, 1)
list111 <- lapply(means.examplesSNSP, `[[`, 1)
#meso
list2 <- lapply(means.examplesSNSO, `[[`, 2)
list22 <- lapply(means.examplesSNSR, `[[`, 2)
list222 <- lapply(means.examplesSNSP, `[[`, 2)
#hyper
list3 <- lapply(means.examplesSNSO, `[[`, 3)
list33 <- lapply(means.examplesSNSR, `[[`, 3)
list333 <- lapply(means.examplesSNSP, `[[`, 3)
#holo
list4 <- lapply(means.examplesSNSO, `[[`, 4)
list44 <- lapply(means.examplesSNSR, `[[`, 4)
list444 <- lapply(means.examplesSNSP, `[[`, 4)
#hypo with 100% VC
list5 <- lapply(means.examplesSNSO100, `[[`, 1)
list55 <- lapply(means.examplesSNSR100, `[[`, 1)
list555 <- lapply(means.examplesSNSP100, `[[`, 1)
#meso with 100% VC
list6 <- lapply(means.examplesSNSO100, `[[`, 2)
list66 <- lapply(means.examplesSNSR100, `[[`, 2)
list666 <- lapply(means.examplesSNSP100, `[[`, 2)
#hyper with 100% VC
list7 <- lapply(means.examplesSNSO100, `[[`, 3)
list77 <- lapply(means.examplesSNSR100, `[[`, 3)
list777 <- lapply(means.examplesSNSP100, `[[`, 3)
#holo with 100% VC
list8 <- lapply(means.examplesSNSO100, `[[`, 4)
list88 <- lapply(means.examplesSNSR100, `[[`, 4)
list888 <- lapply(means.examplesSNSP100, `[[`, 4)

# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year-1) #minus a year for Savanes
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]]),
  hypoendemic100 = list(list5[[2]], list55[[2]], list555[[2]]),
  mesoendemic100 = list(list6[[2]], list66[[2]], list666[[2]]),
  hyperendemic100 = list(list7[[2]], list77[[2]], list777[[2]]),
  holoendemic100 = list(list8[[2]], list88[[2]], list888[[2]])
)

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = SavanesNoEND1
)

# Observed known baseline
known_datasets <- list(
  mesoendemic = SavanesNo2[SavanesNo2$Survey_year >= 2000, ],
  hyperendemic = SavanesHyperNo[SavanesHyperNo$Survey_year >= 2000, ],
  mesoendemic100 = SavanesNo2[SavanesNo2$Survey_year >= 2000, ],
  hyperendemic100 = SavanesHyperNo[SavanesHyperNo$Survey_year >= 2000, ]
)

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)
#rerun the simulation for hyperendemicity100 VC, but instead o f-1 in "round(obs_data$Survey_year-1)", do time +3



###KARA WITHOUT BASELINE###
KaraMerged <- full_join(KaraYes0, KaraNo0)
#now the believed endemicity levels
KaraYesEND1 <- subset(KaraMerged, VILNOM %in% c("ADETEYO","AGBANG2","AGBANSODA","BAOULINSE","BOUGABOU","BOULARE","BAWLESI","KANDJO","KAWA BASSAR II","KOUDJOUKADA","KOUNKOUMBOULE","KPAGBAZIBIYO","KPELOUWAI","KPETAB","LANGA","LEZIYO","NANDOUNGBALE","OURO-GAODE (DAKO)","PESSIDE-ANCIEN","PIYADE","POYO","TCHALOUDE","TCHOLOKOUDE","TOUMBOUA","TOUNDOUNON")) 
KaraYesEND2 <- subset(KaraMerged, VILNOM %in% c("ABOUDA","AGBARADA","BOUNOH","BOWINDO","HALALOMOU (FILANDI)","HOUNDE","KARBONGOU","KASSOU","KAWA","KONFOUH + DIAB","OTI-VILLAGE & BIDJAB","TCHABOUA","TCHIRKPENI(KATCHAMBA)","TCHITCHIKPOLA","KPANGBASSIBIYO","POWAI","KASSI (LANDA)"))
KaraYesEND3 <- subset(KaraMerged, VILNOM %in% c("DANDJESSI","DJAMDE KAWA","HOURTA","KATCHA-KONKOMBA","KAWA-BASSAR","KISSAFO","KOULWERE","KOUTANTAGOU","KOUTANTAGOU &TAPOUNT","KPABTE","MADJATOM","POSSAO","SABOUNDI","SAKPONE","SEKOU-BAS","WARTEMA","ZONE MARAICHERE","WELOUDE (KPAYABOW)", "KADJOL II"))
KaraYesEND4 <- subset(KaraMerged, VILNOM %in% c("GOULBI","KOFFI-FERME", "KOUTOUGOU SOLLA","KPANTIIYAGOU","NARITA/PESSIDE","SOLA","TCHAKASSOU","WASSI","WASSITE","TCHITCHIRA FERME II","AHO-LAO","SIKAN", "PESSIDE FERME+WASITE","TOUGUEL"))
KaraYes1U <- KaraYes1
villages_to_add <- KaraYes2[KaraYes2$VILNOM %in% c("KPESSIDE", "LEON", "ANIMA"), ]
# Add these rows to KaraYes3
KaraYes3U <- rbind(KaraYes3, villages_to_add)
# Remove these villages from the original dataset
KaraYes2U <- KaraYes2[!KaraYes2$VILNOM %in% c("KPESSIDE", "LEON", "ANIMA"), ]
#add Kemini survey in 2014 (and remove it form Centrale SIZ no baseline)
villages_to_add <- CentraleNo0[CentraleNo0$VILNOM %in% c("KEMENI"), ]
# Add these rows to PlateauxNo3
KaraYes2U <- rbind(villages_to_add,KaraYes2U)
KaraYes2U$ENDE <- "Mesoendemic"
#pass Titira from hyper to holoendemic
villages_to_add <- KaraYes2[KaraYes2$VILNOM %in% c("KPESSIDE", "LEON", "ANIMA"), ]
# Add these rows to KaraYes3
KaraYes3U <- KaraYes3[!KaraYes3$VILNOM %in% c("TITIRA"), ]
KaraYes3U <- rbind(KaraYes3U, villages_to_add)
KaraYes3U$ENDE <- "Hyperendemic"
KaraYes4U <- KaraYes4
#HYPOENDEMIC
list1 <- lapply(means.examplesKSOB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesKSRB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesKSPB, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list2 <- lapply(means.examplesKSOB, `[[`, 2)
list22 <- lapply(means.examplesKSRB, `[[`, 2)
list222 <- lapply(means.examplesKSPB, `[[`, 2)
list3 <- lapply(means.examplesKNSOB, `[[`, 3)
list33 <- lapply(means.examplesKNSRB, `[[`, 3)
list333 <- lapply(means.examplesKNSPB, `[[`, 3)
list4 <- lapply(means.examplesKSOB, `[[`, 4)
list44 <- lapply(means.examplesKSRB, `[[`, 4)
list444 <- lapply(means.examplesKSPB, `[[`, 4)

# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]])
)

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = KaraYesEND1,
  mesoendemic = KaraYesEND2,
  hyperendemic = KaraYesEND3,
  holoendemic = KaraYesEND4
)

# Observed known baseline
known_datasets <- list(
  hypoendemic = KaraYes1U[KaraYes1U$Survey_year >= 1990, ],
  mesoendemic = KaraYes2U[KaraYes2U$Survey_year >= 1990, ],
  hyperendemic = KaraYes3U[KaraYes3U$Survey_year >= 1990, ],
  holoendemic = KaraYes4U[KaraYes4U$Survey_year >= 1990, ]
)

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)




###Centrale SIZ WITHOUT BASELINE####
#unknown
CentraleYesEND1 <- subset(CentraleYes0, VILNOM %in% c("GNEZIME")) 
CentraleYesEND3 <- subset(CentraleYes0, VILNOM %in% c("AGBAMASSOUMOU","DANTCHESSI","TCHATOU KOURA","TCHIDAO","NABOUN-KOURA","MOUSSOUKOUDJOU","TCHAKPISSI","TCHETCHEKOU"))
CentraleYesEND4 <- subset(CentraleYes0, VILNOM %in% c("BANDA","ASSAWOH-KOURA","BATTO","KOIDA"))

#now the known endemicity levels
#pass Bouzalo from non-SIZ to SIZ (it is next to Mo in tchaoudjo)
villages_to_add <- CentraleNo3[CentraleNo3$VILNOM %in% c("BOUZALO "), ]
# Add these rows to PlateauxNo3
CentraleYes3B <- rbind(villages_to_add,CentraleYes3B)
#add missing surveys of bouzalo from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9064110/
# Step 1: Identify and copy the row with "BOUZALO" in the "VILNOM" column
bouzalo_row <- CentraleYes3B[CentraleYes3B$VILNOM == "BOUZALO ", ]
# Step 2: Create four new rows based on the copied row
new_rows <- do.call(rbind, replicate(4, bouzalo_row, simplify = FALSE))
# Step 3: Update the specified columns with new data
new_rows$BASNOM <- "Mo"
new_rows$DISNOM <- "Tchaoudjo"
new_rows$Survey_year <- c(1985.246, 1989.165, 1990.821, 1993.165)
new_rows$PREV_cr <- c(63.7, 68.1, 52.9, 7.5)
new_rows$CILower <- c(56.8, 59.6, 43.6, 5.3)
new_rows$CIUpper <- c(70.1, 75.5, 62.1, 10.5)
# Step 4: Bind the new rows to the original dataset
CentraleYes3B <- rbind(CentraleYes3B, new_rows)

#HYPOENDEMIC
list1 <- lapply(means.examplesKSOB, `[[`, 1)
list11 <- lapply(means.examplesKSRB, `[[`, 1)
list111 <- lapply(means.examplesKSPB, `[[`, 1)
#meso
list2 <- lapply(means.examplesKSOB, `[[`, 2)
list22 <- lapply(means.examplesKSRB, `[[`, 2)
list222 <- lapply(means.examplesKSPB, `[[`, 2)
#hyper
list3 <- lapply(means.examplesKSOB, `[[`, 3)
list33 <- lapply(means.examplesKSRB, `[[`, 3)
list333 <- lapply(means.examplesKSPB, `[[`, 3)
#holo
list4 <- lapply(means.examplesKSOB, `[[`, 4)
list44 <- lapply(means.examplesKSRB, `[[`, 4)
list444 <- lapply(means.examplesKSPB, `[[`, 4)

# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year-1) #minus a year for Savanes
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]]))

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = CentraleYesEND1,
  hyperendemic = CentraleYesEND3,
  holoendemic = CentraleYesEND4
)

# Observed known baseline
known_datasets <- list(
  hyperendemic = CentraleYes3B[CentraleYes3B$Survey_year >= 1995, ])

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)




###Centrale non-SIZ WITHOUT BASELINE####
#unknown
CentraleNoEND1 <- subset(CentraleNo0, VILNOM %in% c("BALERIDE","DIKPELEOU/DJAOULLA","KPAKPARASSOU+NGOBO","LAMA WERE-LAOUDA","MOTCHOKPLI","N'KENGBE","TCHEMBERI","YOUROUROU")) 
CentraleNoEND2 <- subset(CentraleNo0, VILNOM %in% c("AKAWOLO","ATCHAVE","BLOU-ELAVAGNON","KPAMBOURE","KPEIDA","OKOU-KOPE","OUDJOMBOI","PAGALA-BOUZIYA","PANLAO","SADA-MONO","SOUKOUNDE","TALABA","YOVO-KOPE"))
CentraleNoEND3 <- subset(CentraleNo0, VILNOM %in% c("AGBANDI-MONO","KATCHALIKADI","OGOUDA & SOMBO","TAKADE","YELOUM BAGNAN"))

#known
CentraleNo1B <- subset(CentraleNo1, DISNOM %in% c("SOTOUBOUA", "TCHAOUDJO") & MesoYes != 1)
CentraleNo1A <- subset(CentraleNo1, DISNOM %in% c("TCHAMBA", "BLITTA") & MesoYes != 1)
CentraleNo1A<- subset(CentraleNo1, MesoYes != 1)
CentraleNo2B <- subset(CentraleNo2, DISNOM %in% c("SOTOUBOUA", "TCHAOUDJO"))
CentraleNo2A <- subset(CentraleNo2, DISNOM %in% c("TCHAMBA", "BLITTA"))
CentraleNo2A<-CentraleNo2
villages_to_add <- E8[E8$VILNOM %in% c("SOMIEDA-LAOUDE MONO"), ]
# Add these rows to PlateauxNo3
CentraleNo2A <- rbind(villages_to_add,CentraleNo2A)
CentraleNo2A$ENDE <- "Mesoendemic"
#now the known endemicity levels
CentraleNo3B <- subset(CentraleNo3, DISNOM %in% c("SOTOUBOUA", "TCHAOUDJO")) #not enough to do biannual
CentraleNo3A <- subset(CentraleNo3, DISNOM %in% c("TCHAMBA", "BLITTA"))
CentraleNo3A<-CentraleNo3
CentraleNo3A$ENDE <- "Hyperendemic"
#Bouzalo is SIZ, remove it
CentraleNo3A <- subset(CentraleNo3A, !(VILNOM == "BOUZALO " )) #KPODJI is hyperendemic

#HYPOENDEMIC
list1 <- lapply(means.examplesKNSO, `[[`, 1)
list11 <- lapply(means.examplesKNSR, `[[`, 1)
list111 <- lapply(means.examplesKNSP, `[[`, 1)
#meso
list2 <- lapply(means.examplesKNSO, `[[`, 2)
list22 <- lapply(means.examplesKNSR, `[[`, 2)
list222 <- lapply(means.examplesKNSP, `[[`, 2)
#hyper
list3 <- lapply(means.examplesKNSO, `[[`, 3)
list33 <- lapply(means.examplesKNSR, `[[`, 3)
list333 <- lapply(means.examplesKNSP, `[[`, 3)
#holo
list4 <- lapply(means.examplesKNSO, `[[`, 4)
list44 <- lapply(means.examplesKNSR, `[[`, 4)
list444 <- lapply(means.examplesKNSP, `[[`, 4)

# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year-1) #minus a year for Savanes
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]]))

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = CentraleNoEND1,
  mesoendemic = CentraleNoEND2,
  hyperendemic = CentraleNoEND3
)

# Observed known baseline
known_datasets <- list(
  hypoendemic = CentraleNo1A[CentraleNo1A$Survey_year >= 1995, ],
  mesoendemic = CentraleNo2A[CentraleNo2A$Survey_year >= 1995, ],
  hyperendemic = CentraleNo3A[CentraleNo3A$Survey_year >= 1995, ])

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)





###Plateaux  WITHOUT BASELINE####
#unknown
PlateauxNoBEND <- subset(PlateauxNo0, DISNOM %in% c("DANYI", "HAHO","OGOU","AMOU") | VILNOM %in% c("KPATI COPE","ATINKPASSA","TSOKPLE")) #not enough to do biannual
PlateauxNoAEND <- subset(PlateauxNo0, DISNOM %in% c("AGOU", "AKÉBOU", "MOYEN-MONO", "EST-MONO","KLOTO","WAWA","KPÉLÉ","Wawa","ANIÉ"))
#annual
#Dont account AMOUGODO and KONTA+AGBATEHI and DAYES-DODZI and IGBOM?DJI, that is already accounted for in the plots with recorded endemicity
PlateauxNoEND1A <- subset(PlateauxNoAEND, VILNOM %in% c("AGOKPLAME","AGOUDOUVOU","ALÉ KOPÉ","BLOUDOKOPE","DEVELEBE","GAVO KOSSI","KPOVENOU","LETSOUKOPE","MANGOTIGOMÉ","NOUGBESSOU KOPE","WOGLOKOPE","ZIONOU","ZOGBEGAN-OGA","TOYIGBO")) 
PlateauxNoEND2A <- subset(PlateauxNoAEND, VILNOM %in% c("AHLON DZINDZI","KOUMASSE","TOME"))
PlateauxNoEND3A <- subset(PlateauxNoAEND, VILNOM %in% c("ANANIVIKODZI","ATEWE-ZONGO","GUIN KOPE","KLO-MAYONDI","KPIME-SEVA","NYIVE","ODOMI ABRA","SUKUL-KPODJI","TUTU ZIONOU"))
#PlateauxNoEND4A <- subset(PlateauxNoAEND, VILNOM %in% c())
#biannual
#Dont account G. AKEME, as it is already accounted as hyperendemic
PlateauxNoEND1B <- subset(PlateauxNoBEND, VILNOM %in% c("AGBOROU KOPE","AGNAM","AKPAKA","ANYAM-KOPE","ASSANTE","BAGAOU","EBAFEI-KOPE","GOTHA ADJA","GOTHA KABYE","GROKOPÉ","HAHONOU","HAOUSSA KPÉDJI","KOME","KPÉDJI","KPÉTÉ MAVA","MEDZE","KPELE KOPE")) 
PlateauxNoEND2B <- subset(PlateauxNoBEND, VILNOM %in% c("AMOU AKPEKPE","AMOUZOUKOPE(DJEMENI)","DJAKPO","OTCHANARI","ATALAKPOTA","TCHÉKÉLÉ","HOUNO-KOPE","FAWUKPE","MAYABA-KOPE"))
PlateauxNoEND3B <- subset(PlateauxNoBEND, VILNOM %in% c("HETRE","GLIVE","ILEKOHAN","DENOU BUMUEBI","MOBA KOPÉ","S. OUTOUALA","TANAGO","AMOUTCHOU","AMOUTO","TOIGBO","ATINKPASSA","EGBOWOU AMOU","IGBOWOU-AMOU","GLELOU+OMOUVA","KPATI COPE","PIDINA","TSOKPLE"))


#now the known endemicity levels
PlateauxNo1 <- PlateauxNo1 %>%
  mutate(DISNOM = ifelse(VILNOM == "GBAGBADJAKOU I", "ANIÉ", DISNOM))
villages_to_add <- PlateauxNo1[PlateauxNo1$VILNOM %in% c("AMOUTA"), ]
# Add these rows to PlateauxNo3
PlateauxNo2 <- rbind(villages_to_add,PlateauxNo2)
PlateauxNo2$ENDE <- "Mesoendemic"
PlateauxNo2 <- PlateauxNo2 %>%
  mutate(PREV_cr = ifelse(VILNOM == "AMOUTA" & PREV_cr=="36.54", PREV_st, PREV_cr))
PlateauxNo2 <- PlateauxNo2 %>%
  mutate(CILower = ifelse(VILNOM == "AMOUTA" & PREV_cr=="47.2", 37.4, CILower))
PlateauxNo2 <- PlateauxNo2 %>%
  mutate(CIUpper = ifelse(VILNOM == "AMOUTA" & PREV_cr=="47.2", 57.2, CIUpper))
PlateauxNo1 <- subset(PlateauxNo1, !(VILNOM == "AMOUTA")) #Remove an unnusual value of village Konta
villages_to_add <- PlateauxNo1[PlateauxNo1$VILNOM %in% c("KEMEDISSO"), ]
# Add these rows to PlateauxNo3
PlateauxNo2 <- rbind(villages_to_add,PlateauxNo2)
PlateauxNo2$ENDE <- "Mesoendemic"
PlateauxNo2 <- PlateauxNo2 %>%
  mutate(PREV_cr = ifelse(VILNOM == "KEMEDISSO" & PREV_cr=="39.05", PREV_st, PREV_cr))
PlateauxNo2 <- PlateauxNo2 %>%
  mutate(CILower = ifelse(VILNOM == "KEMEDISSO" & PREV_cr=="40.60", 31.3, CILower))
PlateauxNo2 <- PlateauxNo2 %>%
  mutate(CIUpper = ifelse(VILNOM == "KEMEDISSO" & PREV_cr=="40.60", 50.6, CIUpper))
PlateauxNo1 <- subset(PlateauxNo1, !(VILNOM == "KEMEDISSO")) #Remove an unnusual value of village Konta
villages_to_add <- PlateauxNo2[PlateauxNo2$VILNOM %in% c("ATOME"), ]
# Add these rows to PlateauxNo3
PlateauxNo3 <- rbind( villages_to_add,PlateauxNo3)
PlateauxNo3$ENDE <- "Hyperendemic"
PlateauxNo3 <- PlateauxNo3 %>%
  mutate(PREV_cr = ifelse(VILNOM == "ATOME" & PREV_cr=="59.76", PREV_st, PREV_cr))
PlateauxNo3 <- PlateauxNo3 %>%
  mutate(CILower = ifelse(VILNOM == "ATOME" & PREV_cr=="65.4", 61.3, CILower))
PlateauxNo3 <- PlateauxNo3 %>%
  mutate(CIUpper = ifelse(VILNOM == "ATOME" & PREV_cr=="65.4", 69.3, CIUpper))
PlateauxNo2 <- subset(PlateauxNo2, !(VILNOM == "ATOME")) #Remove an unnusual value of village Konta
PlateauxNo1B <- subset(PlateauxNo1, DISNOM %in% c("DANYI", "HAHO","OGOU","AMOU"))
PlateauxNo1A <- subset(PlateauxNo1, DISNOM %in% c("AGOU", "AKÉBOU","MOYEN-MONO", "EST-MONO","KLOTO","WAWA","ANIÉ"))
PlateauxNo2B <- subset(PlateauxNo2, DISNOM %in% c("AMOU","DANYI", "HAHO","OGOU")) 
PlateauxNo2A <- subset(PlateauxNo2, DISNOM %in% c("AGOU", "AKÉBOU","MOYEN-MONO", "EST-MONO","KLOTO","WAWA", "Kpele","ANIÉ"))
#PlateauxNo2A <- subset(PlateauxNo2A, !(VILNOM == "KONTA" & POS_BASE == 85)) #Remove an unnusual value of village Konta
PlateauxNo2A <- PlateauxNo2A %>% #instead, correct year from 2002 to 1992
  mutate(Survey_year = ifelse(VILNOM == "KONTA"& POS_BASE == 85, 1992, Survey_year))
PlateauxNo2A <- subset(PlateauxNo2A, !(VILNOM == "KPODJI")) 
#oBE AND Obetodji are the same. add the latter to the dataset
library(dplyr)
# Extract the "SAMOMONI" row from SavanesNo5
obe_row <- E8 %>% filter(VILNOM == "PYACOPE(OBETODJI)")
# Append the "PYACOPE(OBETODJI)" row to the new dataset
PlateauxNo2A <- bind_rows(PlateauxNo2A, obe_row)
#SavanesHyperNo <- data.frame(samomoni_row)
PlateauxNo2A <- PlateauxNo2A %>%
  mutate(ENDE = ifelse(VILNOM == "PYACOPE(OBETODJI)", "Mesoendemic", ENDE))
PlateauxNo2B <- subset(PlateauxNo2B, !(VILNOM == "ADJABOULOUKOUKOPE" & CECI == 9999  &  POS_BASE == 21)) #Remove an unnusual value, same as the survey 1997.715
PlateauxNo2B <- subset(PlateauxNo2B, !(VILNOM == "KPODJI" )) #KPODJI is hyperendemic
PlateauxNo3 <- PlateauxNo3 %>%
  mutate(DISNOM = ifelse(VILNOM == "KAMALO-KOPE", "ANIÉ", DISNOM))
PlateauxNo3 <- PlateauxNo3 %>%
  mutate(DISNOM = ifelse(VILNOM == "AROUKAKOPE", "OGOU", DISNOM))
PlateauxNo3 <- PlateauxNo3 %>%
  mutate(DISNOM = ifelse(VILNOM == "KONIGBO", "ANIÉ", DISNOM))
PlateauxNo3B <- subset(PlateauxNo3, DISNOM %in% c("AMOU","DANYI","HAHO","OGOU")) 
PlateauxNo3A <- subset(PlateauxNo3, DISNOM %in% c("AGOU", "AKÉBOU", "MOYEN-MONO", "EST-MONO","KLOTO","WAWA", "Kpele", "ANIÉ"))
#confirm no villages were left out --> checked
PlateauxNo3C <- subset(PlateauxNo3, !(DISNOM == "AMOU"|DISNOM =="DANYI"|DISNOM == "HAHO"|DISNOM =="OGOU"|DISNOM =="AGOU"|DISNOM == "AKÉBOU"|DISNOM == "ANIÉ"|DISNOM == "MOYEN-MONO"|DISNOM == "EST-MONO"|DISNOM =="KLOTO"|DISNOM =="WAWA"|DISNOM == "Kpele")) #Remove an unnusual value of village Konta
#pass KPODJI from meso to hyperendemic
villages_to_add <- PlateauxNo2[PlateauxNo2$VILNOM %in% c("KPODJI"), ]
# Add these rows to PlateauxNo3
PlateauxNo3B <- rbind(villages_to_add,PlateauxNo3B)
PlateauxNo3B$ENDE <- "Hyperendemic"
#Pass DAYES-DODZI from no endemicity to hyper
villages_to_add <- PlateauxNo0[PlateauxNo0$VILNOM %in% c("DAYES-DODZI"), ]
# Add these rows to PlateauxNo3
PlateauxNo3A <- rbind(villages_to_add,PlateauxNo3A)
PlateauxNo3A$ENDE <- "Hyperendemic"
villages_to_add <- E8[E8$VILNOM %in% c("KESSIBO-DZODZI (DAYI DODJI)"), ]
# Add these rows to PlateauxNo3
PlateauxNo3A <- rbind(villages_to_add,PlateauxNo3A)
PlateauxNo3A$ENDE <- "Hyperendemic"
#pass ALAMASSOu from holoendemic to hyperendemic
villages_to_add <- PlateauxNo4[PlateauxNo4$VILNOM %in% c("ALAMASSOU"), ]
# Add these rows to PlateauxNo3
PlateauxNo3B <- rbind(villages_to_add,PlateauxNo3B)
PlateauxNo3B$ENDE <- "Hyperendemic"
#add new survey of Dayes Djodzi
Dayes_row <- E8 %>% filter(VILNOM == "DAYES-DODZI")
new_Dayes <- Dayes_row %>%
  mutate(Survey_year = 2000, 
         PREV_cr = 2.2, 
         PREV_st = NA, 
         CILower = 0.7, 
         CIUpper = 5.8)
# Append the new row to the dataset
PlateauxNo3A <- bind_rows(PlateauxNo3A, new_Dayes)
PlateauxNo3A <- PlateauxNo3A %>%
  mutate(ENDE = ifelse(VILNOM == "DAYES-DODZI", "Hyperendemic", ENDE))

#HYPOENDEMIC
list1 <- lapply(means.examplesPO, `[[`, 1)
list11 <- lapply(means.examplesPR, `[[`, 1)
list111 <- lapply(means.examplesPP, `[[`, 1)
#meso
list2 <- lapply(means.examplesPO, `[[`, 2)
list22 <- lapply(means.examplesPR, `[[`, 2)
list222 <- lapply(means.examplesPP, `[[`, 2)
#hyper
list3 <- lapply(means.examplesPO, `[[`, 3)
list33 <- lapply(means.examplesPR, `[[`, 3)
list333 <- lapply(means.examplesPP, `[[`, 3)
#holo
list4 <- lapply(means.examplesPO, `[[`, 4)
list44 <- lapply(means.examplesPR, `[[`, 4)
list444 <- lapply(means.examplesPP, `[[`, 4)
#hypo with biannual
list5 <- lapply(means.examplesPOB, `[[`, 1)
list55 <- lapply(means.examplesPRB, `[[`, 1)
list555 <- lapply(means.examplesPPB, `[[`, 1)
#meso with biannual
list6 <- lapply(means.examplesPOB, `[[`, 2)
list66 <- lapply(means.examplesPRB, `[[`, 2)
list666 <- lapply(means.examplesPPB, `[[`, 2)
#hyper with biannual
list7 <- lapply(means.examplesPOB, `[[`, 3)
list77 <- lapply(means.examplesPRB, `[[`, 3)
list777 <- lapply(means.examplesPPB, `[[`, 3)
#holo with biannual
list8 <- lapply(means.examplesPOB, `[[`, 4)
list88 <- lapply(means.examplesPRB, `[[`, 4)
list888 <- lapply(means.examplesPPB, `[[`, 4)

# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year-1) #minus a year for Savanes
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]]),
  hypoendemicbiannual = list(list5[[2]], list55[[2]], list555[[2]]),
  mesoendemicbiannual = list(list6[[2]], list66[[2]], list666[[2]]),
  hyperendemicbiannual = list(list7[[2]], list77[[2]], list777[[2]]),
  holoendemicbiannual = list(list8[[2]], list88[[2]], list888[[2]])
)

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = PlateauxNoEND1A,
  Mesoendemic = PlateauxNoEND2A,
  Hyperendemic = PlateauxNoEND3A,
  hypoendemicbiannual = PlateauxNoEND1B,
  mesoendemicbiannual = PlateauxNoEND2B,
  hyperendemicbiannual = PlateauxNoEND3B
)

# Observed known baseline
known_datasets <- list(
  hypoendemic = PlateauxNo1A[PlateauxNo1A$Survey_year >= 1992, ],
  mesoendemic = PlateauxNo2A[PlateauxNo2A$Survey_year >= 1992, ],
  hyperendemic = PlateauxNo3A[PlateauxNo3A$Survey_year >= 1992, ],
  hypoendemicbiannual = PlateauxNo1B[PlateauxNo1B$Survey_year >= 1992, ],
  mesoendemicbiannual = PlateauxNo2B[PlateauxNo2B$Survey_year >= 1992, ],
  hyperendemicbiannual = PlateauxNo3B[PlateauxNo3B$Survey_year >= 1992, ]
)

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)
#rerun the simulation for hyperendemicity100 VC, but instead o f-1 in "round(obs_data$Survey_year-1)", do time +3






###Maritime WITHOUT BASELINE###
#now the believed endemicity levels
MaritimeNoEND1 <- subset(MaritimeNo0, VILNOM %in% c("ADIKPE","AGOMENOU","AGOTIME","AGOTO","AKATI ZOGBE","AKE-KONDJI","ALOKPA","ATIKPATAFO","AVEGODOE","BATOE","DEKPO","DROUGBOKOPE","ESSE KOLEVE","FRANGADOUA","HAHO-KPODJI","KPEHO","MOUSSOUHOE","NOUSSOUKOPE","SAKPA-KPENSI","TOFA-KOPE","TOKPLI","TOKPLI (ZOUME)","TOVE","VOULE")) 
MaritimeNoEND2 <- subset(MaritimeNo0, VILNOM %in% c("AFOMONOU","AFOKONOU","GBANDIDI","GOGOKONDJI"))
MaritimeNoEND3 <- subset(MaritimeNo0, VILNOM %in% c("AFANGADJI","KAYIDO","DJREKPON","TOGBA","DZREKPON","MAWUSSOU","LAKATA-KONDJI"))
villages_to_add <- E8[E8$VILNOM %in% c("KONTA+AGBATEHI"), ]
# Add these rows to PlateauxNo3
MaritimeNoEND3 <- rbind(villages_to_add,MaritimeNoEND3)
MaritimeNoEND3$ENDE <- "Hyperendemic"

#now the known endemicity levels
MaritimeNo1
MaritimeNo3

#HYPOENDEMIC
list1 <- lapply(means.examplesMO, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list11 <- lapply(means.examplesMR, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list111 <- lapply(means.examplesMP, `[[`, 1) #A list with 4 numeric columns (plus the fifth). We just need to select the second one, that corresponds to the prevalence
list2 <- lapply(means.examplesMO, `[[`, 2)
list22 <- lapply(means.examplesMR, `[[`, 2)
list222 <- lapply(means.examplesMP, `[[`, 2)
list3 <- lapply(means.examplesMO, `[[`, 3)
list33 <- lapply(means.examplesMR, `[[`, 3)
list333 <- lapply(means.examplesMP, `[[`, 3)
list4 <- lapply(means.examplesMO, `[[`, 4)
list44 <- lapply(means.examplesMR, `[[`, 4)
list444 <- lapply(means.examplesMP, `[[`, 4)

# Function to compute MSE, as before
compute_MSE <- function(obs_data, sim_data, obs_label, sim_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  
  for (i in 1:3){ # three scenarios
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = obs_label,
        Simulated_Endemicity = sim_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Unknown"
      ))
    }
  }
  return(mse_results)
}

# Function for known endemicity (compare only to the correct scenario)
compute_MSE_known <- function(obs_data, sim_data, endemic_label) {
  scenarios <- c("Enhanced", "Reference", "Minimal")
  mse_results <- data.frame()
  for (i in 1:3) {
    sim_scenario <- sim_data[[i]] * 100
    df_sim <- data.frame(Year = floor(time), Prevalence = sim_scenario)
    annual_mean <- aggregate(Prevalence ~ Year, data = df_sim, mean)
    obs_data$Year <- round(obs_data$Survey_year)
    merged_df <- merge(obs_data, annual_mean, by = "Year")
    if(nrow(merged_df) > 0) {
      mse <- mean((merged_df$PREV_cr - merged_df$Prevalence)^2)
      mse_results <- rbind(mse_results, data.frame(
        Observed_Endemicity = endemic_label,
        Simulated_Endemicity = endemic_label,
        Scenario = scenarios[i],
        MSE = mse,
        Baseline_Status = "Known"
      ))
    }
  }
  return(mse_results)
}

# Simulations as before
endemicity_simulations <- list(
  hypoendemic = list(list1[[2]], list11[[2]], list111[[2]]),
  mesoendemic = list(list2[[2]], list22[[2]], list222[[2]]),
  hyperendemic = list(list3[[2]], list33[[2]], list333[[2]]),
  holoendemic = list(list4[[2]], list44[[2]], list444[[2]])
)

# Observed unknown baseline
observed_datasets <- list(
  hypoendemic = MaritimeNoEND1,
  mesoendemic = MaritimeNoEND2,
  hyperendemic = MaritimeNoEND3)

# Observed known baseline
known_datasets <- list(
  hypoendemic = MaritimeNo1[MaritimeNo1$Survey_year >= 1992, ],
  hyperendemic = MaritimeNo3[MaritimeNo3$Survey_year >= 1992, ]
)

# Calculate for unknown
final_MSE_comparison <- data.frame()
for (obs_name in names(observed_datasets)){
  obs_data <- observed_datasets[[obs_name]]
  if(nrow(obs_data) == 0) next
  for (sim_name in names(endemicity_simulations)){
    sim_data <- endemicity_simulations[[sim_name]]
    mse_res <- compute_MSE(obs_data, sim_data, obs_name, sim_name)
    final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
  }
}

# Calculate for known
for (known_name in names(known_datasets)){
  obs_data <- known_datasets[[known_name]]
  if(nrow(obs_data) == 0) next
  sim_data <- endemicity_simulations[[known_name]]
  mse_res <- compute_MSE_known(obs_data, sim_data, known_name)
  final_MSE_comparison <- rbind(final_MSE_comparison, mse_res)
}

# Show results
print(final_MSE_comparison)







